# -*- coding: utf-8 -*-
"""
AAV qPCR Titer Calculator 
- Mirrors the yellow cells in your spreadsheet
- Computes ds/ss molecules from DNA mass
- Builds 8-point 10× standards and fits Ct = m * log10(copies) + b
- Converts sample Ct replicates to vg/µL and vg/mL using your dilution scheme
"""

import math
import numpy as np
import pandas as pd

# Optional plotting (safe to disable if not installed)
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

AVOGADRO = 6.022e23  # mol^-1 (Avogadro's constant)

# =========================
# 0) USER CONFIG — NOW USING COMMAND LINE ARGUMENTS
# =========================

import argparse

parser = argparse.ArgumentParser(
    description="AAV qPCR Titer Calculator — command line version"
)

# --- DNA Standard block ---
parser.add_argument("--dna_name", type=str, required=True)
parser.add_argument("--dna_size", type=float, required=True)
parser.add_argument("--bp_weight", type=float, default=650.0, help="average molar mass per bp for dsDNA (g/mol/bp)")
parser.add_argument("--std_stock", type=float, required=True,
                    help="Standard DNA concentration in ng/µL")
parser.add_argument("--v_qpcr", type=float, required=True,
                    help="Volume loaded into qPCR (µL)")

# --- Dilution Factor block ---
parser.add_argument("--pk_total", type=float, required=True)
parser.add_argument("--virus_input", type=float, required=True)
parser.add_argument("--protease", type=float, required=True)
parser.add_argument("--template", type=float, required=True)

# --- Standard dilutions (variable length) ---
parser.add_argument("--dilutions", nargs="+", required=True,
                    help="Standard dilution series, e.g. 1e-1 1e-2 1e-3")

# --- Standard Ct values ---
parser.add_argument("--std_cts", nargs="+", type=float, required=True,
                    help="Ct means OR single replicate per dilution point.")

# --- Sample Ct input (variable number of samples)
# Format: --sample SampleName Ct1 Ct2 Ct3  (repeatable)
parser.add_argument("--sample", action="append", nargs="+",
                    help="Provide sample name followed by its Ct values.")

parser.add_argument("--out_csv", help="optional output CSV filename") 
args = parser.parse_args()

# Assign variables back to your original names
DNA_NAME = args.dna_name
DNA_SIZE_BP = args.dna_size
STD_STOCK_NG_PER_UL = args.std_stock
V_TEMPLATE_QPCR_ML = args.v_qpcr / 1000.0   # convert µL → mL

PK_TOTAL_UL = args.pk_total
VIRUS_INPUT_UL = args.virus_input
PROTEASE_DILUTION = args.protease
TEMPLATE_PER_WELL_UL = args.template

# Convert dilutions to float list
STANDARD_DILUTIONS = [float(d) for d in args.dilutions]

# Standard Ct must match
STD_CTS_TRIPLICATE = [[ct] for ct in args.std_cts]

# Samples dict
SAMPLES_CTS = {}
if args.sample:
    for block in args.sample:
        name = block[0]
        cts = [float(x) for x in block[1:]]
        SAMPLES_CTS[name] = cts

STD_FROM_CSV = None
SAMPLES_FROM_CSV = None
MAKE_PLOT = True


# =========================
# Helper functions
# =========================

def backcalc_factor(pk_total_ul, virus_input_ul, protease_dilution, template_per_well_ul):
    """
    Convert copies measured per qPCR well back to vg/µL in the original stock.
    Factor = (PK_total / virus_input) * (protease_dilution) * (1 / template_per_well)
    """
    return (pk_total_ul / virus_input_ul) * protease_dilution * (1.0 / template_per_well_ul)

def standard_mass_g_used(std_stock_ng_per_ul, v_template_qpcr_ml):
    """
    Mass of standard DNA loaded into one qPCR well (g).
    std_stock_ng_per_ul (ng/µL) × volume (µL) × 1e-9.
    Note: v_template_qpcr_ml is given in mL, so we convert to µL (×1000).
    """
    vol_uL = v_template_qpcr_ml * 1000.0
    return std_stock_ng_per_ul * 1e-9 * vol_uL

def molecules_from_mass(mass_g, bp_len, bp_weight=650.0, ss=True):
    """
    Convert DNA mass (g) to number of molecules (copies).
    moles = mass / (bp_len * bp_weight), copies = moles * Avogadro.
    If ss=True, multiply by 2 to convert ds-molecules to ss copies
    (AAV genome is ssDNA; most plasmid standards are dsDNA).
    """
    mol = mass_g / (bp_len * bp_weight)
    copies_ds = mol * AVOGADRO
    return copies_ds * 2.0 if ss else copies_ds

def fit_standard_curve(log10_copies, ct_means):
    """
    Fit Ct = m * log10(copies) + b.
    Returns slope m, intercept b, R^2, and amplification efficiency.
    Efficiency = 10^(-1/m) - 1
    """
    x = np.asarray(log10_copies, dtype=float)
    y = np.asarray(ct_means, dtype=float)
    m, b = np.polyfit(x, y, 1)
    y_pred = m * x + b
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan
    efficiency = 10 ** (-1.0 / m) - 1.0
    return m, b, r2, efficiency

def copies_from_ct(ct, m, b):
    """Invert the standard curve to get copies per well from Ct."""
    return 10 ** ((ct - b) / m)

def read_standards_csv(path):
    """
    Read standards from CSV with columns like:
      Dilution,Ct1,Ct2,Ct3
    'Dilution' is relative to the undiluted standard (e.g., 1e-8 ... 1e-1).
    """
    df = pd.read_csv(path)
    ct_cols = [c for c in df.columns if c.lower().startswith("ct")]
    ct_means = df[ct_cols].astype(float).mean(axis=1).tolist()
    dilutions = df["Dilution"].astype(float).tolist()
    return dilutions, ct_means

def read_samples_csv(path):
    """
    Read samples from CSV with columns like:
      Sample,Ct1,Ct2,Ct3,...
    Returns dict: {sample_name: [replicate Cts]}
    """
    df = pd.read_csv(path)
    ct_cols = [c for c in df.columns if c.lower().startswith("ct")]
    samples = {}
    for _, r in df.iterrows():
        name = str(r["Sample"])
        vals = [v for v in r[ct_cols].astype(float).tolist() if not math.isnan(v)]
        if vals:
            samples[name] = vals
    return samples

# =========================
# 1) DNA Standard block: mass -> ds/ss molecules
# =========================

mw_g_per_mol = DNA_SIZE_BP * args.bp_weight
mass_g_qpcr = standard_mass_g_used(STD_STOCK_NG_PER_UL, V_TEMPLATE_QPCR_ML)
ds_molecules = molecules_from_mass(mass_g_qpcr, DNA_SIZE_BP, args.bp_weight, ss=False)
ss_molecules = molecules_from_mass(mass_g_qpcr, DNA_SIZE_BP, args.bp_weight, ss=True)

print("=== DNA Standard block ===")
print(f"Name: {DNA_NAME}")
print(f"Size (bp): {DNA_SIZE_BP}")
print(f"Molecular weight (g/mol): {mw_g_per_mol:.0f}")
print(f"AAV genome used in qPCR (g): {mass_g_qpcr:.3e}")
print(f"ds AAV genomic DNA molecules in qPCR: {ds_molecules:.3e}")
print(f"ss AAV genomic DNA molecules in qPCR: {ss_molecules:.3e}")

# =========================
# 2) Dilution Factor (copies/well -> vg/µL)
# =========================

factor = backcalc_factor(PK_TOTAL_UL, VIRUS_INPUT_UL, PROTEASE_DILUTION, TEMPLATE_PER_WELL_UL)
print("\n=== Dilution Factor Calculation ===")
print(f"PK total / virus input = {PK_TOTAL_UL}/{VIRUS_INPUT_UL}")
print(f"Protease treated virus dilution = 1:{PROTEASE_DILUTION:.0f}")
print(f"Template per well (µL) = {TEMPLATE_PER_WELL_UL}")
print(f"Dilution Factor (copies/well -> vg/µL) = {factor:.0f}")

# =========================
# 3) Standards: copies/well & Ct, fit the curve
# =========================

# Compute true copies/well for each dilution point
std_copies = [ss_molecules * d for d in STANDARD_DILUTIONS]
std_log10 = [math.log10(c) for c in std_copies]

if STD_FROM_CSV:
    dilutions_csv, ct_means_csv = read_standards_csv(STD_FROM_CSV)
    STANDARD_DILUTIONS = dilutions_csv
    std_ct_means = ct_means_csv
else:
    # No Ct provided
    if not STD_CTS_TRIPLICATE:
        raise ValueError("ERROR: STD_CTS_TRIPLICATE is empty. Please provide Ct values.")
    
    # Check matching lengths
    if len(STANDARD_DILUTIONS) != len(STD_CTS_TRIPLICATE):
        raise ValueError(
            f"ERROR: Dilution points = {len(STANDARD_DILUTIONS)}, "
            f"but Ct rows = {len(STD_CTS_TRIPLICATE)}. These must match."
        )
    
    # Compute Ct means from user-provided replicates
    std_ct_means = [float(np.mean(rep)) for rep in STD_CTS_TRIPLICATE]


std_table = pd.DataFrame({
    "Dilution": STANDARD_DILUTIONS if not STD_FROM_CSV else dilutions_csv,
    "copies_per_well": std_copies,
    "log10_copies": std_log10,
    "Ct_mean": std_ct_means
})

# Fit only if we have >= 3 standard points with Ct
fit_mask = ~np.isnan(std_table["Ct_mean"].values)
if fit_mask.sum() >= 3:
    m, b, r2, eff = fit_standard_curve(std_table.loc[fit_mask, "log10_copies"].values,
                                       std_table.loc[fit_mask, "Ct_mean"].values)
    print("\n=== Standard curve ===")
    print(f"Ct = {m:.4f} * log10(copies) + {b:.4f}")
    print(f"R² = {r2:.5f}")
    print(f"Efficiency = {(eff*100):.2f}%")
else:
    m = b = r2 = eff = np.nan
    print("\n[WARN] Not enough standard points to fit the curve.")

# Optional plot
if MAKE_PLOT and HAS_MPL and fit_mask.sum() >= 3:
    x = std_table.loc[fit_mask, "log10_copies"].values
    y = std_table.loc[fit_mask, "Ct_mean"].values
    xp = np.linspace(x.min() - 0.5, x.max() + 0.5, 100)
    yp = m * xp + b
    plt.figure(figsize=(5.6, 4.2))
    plt.scatter(x, y)
    plt.plot(xp, yp)
    plt.gca().invert_yaxis()  # visual convention: lower Ct toward top
    plt.xlabel("log10(copies per well)")
    plt.ylabel("Ct")
    plt.title("Standards: Ct vs. log10(copies)")
    plt.tight_layout()
    try:
        plt.show()
    except Exception:
        pass

print("\n=== Standards table (preview) ===")
with pd.option_context('display.precision', 3):
    print(std_table)


# =========================
# 4) Samples: Ct -> copies -> vg/µL, vg/mL
# =========================

if SAMPLES_FROM_CSV:
    samples = read_samples_csv(SAMPLES_FROM_CSV)
else:
    samples = SAMPLES_CTS

rows = []
if not np.isnan(m):
    for name, cts in samples.items():
        arr = np.array(cts, dtype=float)
        arr = arr[~np.isnan(arr)]
        if arr.size == 0:
            continue
        ct_mean = float(arr.mean())
        copies = copies_from_ct(ct_mean, m, b)
        vg_per_uL = copies * factor
        vg_per_mL = vg_per_uL * 1e3
        rows.append({
            "Sample": name,
            "Ct_mean": ct_mean,
            "copies_per_well": copies,
            "vg_per_uL": vg_per_uL,
            "vg_per_mL": vg_per_mL
        })

if rows:
    out_df = pd.DataFrame(rows)
    print("\n=== Sample titers ===")
    with pd.option_context('display.precision', 3):
        print(out_df)

    # Save to CSV if requested
    if args.out_csv:
        out_df.to_csv(args.out_csv, index=False)
        print(f"\n[INFO] Sample titers saved to CSV: {args.out_csv}")

else:
    print("\n[INFO] No sample rows computed yet (no samples or no fitted curve).")
  
