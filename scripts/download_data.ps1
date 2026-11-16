param([string]$Python = $env:PYTHON)

$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent $PSScriptRoot
$DataDir = Join-Path $ProjectRoot "data"

if (-not $Python) {
    $VenvPython = Join-Path $ProjectRoot '.venv\Scripts\python.exe'
    if (Test-Path -LiteralPath $VenvPython -PathType Leaf) {
        $Python = $VenvPython
    } else {
        $Python = 'python'
    }
}
if (-not (Get-Command $Python -ErrorAction SilentlyContinue)) {
    throw "Python not found: $Python"
}

$CoelicolorDir = Join-Path $DataDir "S_coelicolor"
$LividansDir = Join-Path $DataDir "S_lividans"
$RoseumDir = Join-Path $DataDir "S_roseum"

New-Item -ItemType Directory -Force -Path $CoelicolorDir, $LividansDir, $RoseumDir | Out-Null

function Download-File {
    param(
        [string]$Url,
        [string]$Output
    )

    Write-Host "Downloading $(Split-Path $Output -Leaf)"

    Invoke-WebRequest -Uri $Url -OutFile $Output
}

# Streptomyces coelicolor

Download-File `
    "https://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-11-578/MediaObjects/12864_2010_9973_MOESM2_ESM.CSV" `
    (Join-Path $CoelicolorDir "gene_pairs.CSV")

Download-File `
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AL645882.2&rettype=gbwithparts&retmode=text" `
    (Join-Path $CoelicolorDir "sequence.gb")


# Streptomyces lividans

Download-File `
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP009124.1&rettype=fasta&retmode=text" `
    (Join-Path $LividansDir "genome.fasta")

Download-File `
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/6742748/bin/Table_4.XLSX" `
    (Join-Path $LividansDir "TSS_table.xlsx")


# Streptosporangium roseum

$BaseUrl = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/024/865/GCF_000024865.1_ASM2486v1"
$Prefix = "GCF_000024865.1_ASM2486v1"

foreach ($Filename in @(
    "genomic.fna",
    "cds_from_genomic.fna",
    "feature_table.txt"
)) {
    $GzipPath = Join-Path $RoseumDir "$Filename.gz"
    $OutputPath = Join-Path $RoseumDir $Filename

    Download-File "$BaseUrl/${Prefix}_${Filename}.gz" $GzipPath

    $InputStream = [System.IO.File]::OpenRead($GzipPath)
    $OutputStream = [System.IO.File]::Create($OutputPath)

    try {
        $GzipStream = [System.IO.Compression.GZipStream]::new(
            $InputStream,
            [System.IO.Compression.CompressionMode]::Decompress
        )

        try {
            $GzipStream.CopyTo($OutputStream)
        }
        finally {
            $GzipStream.Dispose()
        }
    }
    finally {
        $InputStream.Dispose()
        $OutputStream.Dispose()
    }

    Remove-Item $GzipPath
}

& $Python (Join-Path $ProjectRoot 'scripts/prepare_cds.py') `
    (Join-Path $RoseumDir 'feature_table.txt') `
    (Join-Path $RoseumDir 'CDS_with_protein.csv')
if ($LASTEXITCODE -ne 0) { throw 'CDS table preparation failed.' }

Write-Host "All source datasets downloaded and CDS table prepared."
