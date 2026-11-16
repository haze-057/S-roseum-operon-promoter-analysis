# Run from PowerShell: .\scripts\run_promoter_analysis.ps1
# Requires WSL with Bash, MEME, FIMO, and the Linux Python dependencies installed.
# A Windows virtual environment cannot be used by the Linux runner.
param(
    [string]$Distribution,
    [string]$LinuxPython
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw 'WSL is required for the promoter pipeline. Install WSL and MEME Suite first.'
}

$WslArguments = @()
if ($Distribution) {
    $WslArguments += @('--distribution', $Distribution)
}

$LinuxRoot = & wsl.exe @WslArguments --exec wslpath -a -u $ProjectRoot
if ($LASTEXITCODE -ne 0 -or -not $LinuxRoot) {
    throw 'Could not convert the project path to a WSL path.'
}
$LinuxRoot = ($LinuxRoot -join "`n").Trim()

# --cd preserves paths containing spaces without constructing shell commands.
$PipelineArguments = @('--cd', $LinuxRoot, '--exec', 'env')
if ($LinuxPython) {
    $PipelineArguments += "PYTHON=$LinuxPython"
}
$PipelineArguments += @('bash', 'scripts/run_promoter_analysis.sh')

Write-Host 'Running promoter analysis in WSL'
& wsl.exe @WslArguments @PipelineArguments
if ($LASTEXITCODE -ne 0) {
    throw 'Promoter analysis failed. Check the WSL dependencies and inputs reported above.'
}
