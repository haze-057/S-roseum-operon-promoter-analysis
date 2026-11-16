# Run from PowerShell: .\scripts\run_likehood_analysis.ps1
param(
    [string]$Python = $env:PYTHON
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$PreviousLocation = Get-Location

try {
    Set-Location -LiteralPath $ProjectRoot

    if (-not $Python) {
        $VenvPython = Join-Path $ProjectRoot '.venv\Scripts\python.exe'
        if (Test-Path -LiteralPath $VenvPython -PathType Leaf) {
            $Python = $VenvPython
        } else {
            $Python = 'python'
        }
    }
    if (-not (Get-Command $Python -ErrorAction SilentlyContinue)) {
        throw "Python not found: $Python. Use -Python C:\path\to\python.exe."
    }

    & $Python -c 'import pandas, numpy, Bio, matplotlib'
    if ($LASTEXITCODE -ne 0) { throw 'Python dependency check failed.' }

    foreach ($InputPath in @(
        'data/S_coelicolor/sequence.gb',
        'data/S_coelicolor/gene_pairs.CSV'
    )) {
        if (-not (Test-Path -LiteralPath $InputPath -PathType Leaf)) {
            throw "Required input missing: $InputPath"
        }
        if ((Get-Item -LiteralPath $InputPath).Length -eq 0) {
            throw "Required input empty: $InputPath"
        }
    }

    New-Item -ItemType Directory -Force -Path 'results/operon_prediction_ll' | Out-Null
    # Set $env:MPLBACKEND = 'Agg' to run without plot windows.
    Write-Host 'Calculating inter-CDS distance likelihoods'
    & $Python 'scripts/likehood_analysis/ll_calc.py'
    if ($LASTEXITCODE -ne 0) { throw 'Likelihood analysis failed.' }

    Write-Host "Analysis complete: $ProjectRoot/results/operon_prediction_ll/likehood.json"
} finally {
    Set-Location -LiteralPath $PreviousLocation
}
