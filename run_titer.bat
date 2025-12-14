@echo off
REM === Auto-run AAV qPCR Titer Calculator ===
REM Change directory to the folder where this script is located
cd /d "%~dp0"

REM Run titer.exe with parameters
titer.exe --dna_name 20251031 --dna_size 5882 --std_stock 1.74 --v_qpcr 2 --pk_total 227 --virus_input 2 --protease 300 --template 2 --dilutions 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 --std_cts 12.904 17.079 20.851 24.323 27.759 29.646 --sample AAV-prep-1 21.06 21.32 --sample AAV-prep-2 24.94 24.86 --out_csv titer_results.csv

echo.
echo ✅ Analysis complete! Results have been saved as titer_results.csv
pause
