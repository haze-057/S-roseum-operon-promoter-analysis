# About

This project is designed to identify putative operons in the *Streptosporangium roseum* genome. It includes modules for log-likelihood-based operon prediction, promoter identification, and the identification of promoter–operon pairs. Each module can be executed independently using a dedicated script.

## Requirements

The analysis requires Python 3 and the external bioinformatics software MEME Suite, including MEME and FIMO.

Python dependencies are listed in `requirements.txt` and can be installed using:

```bash
python3 -m pip install -r requirements.txt
```

On Windows, use:

```powershell
py -m pip install -r requirements.txt
```

MEME Suite must be installed separately. Installation instructions are available in the [official MEME Suite documentation](https://meme-suite.org/meme/doc/install.html).

For Windows users, MEME Suite requires a compatible environment, such as WSL2 or Docker.

All analyses should be executed from the project root directory.

## Running on Linux/macOS

First, clone the GitHub repository:

```bash
git clone https://github.com/haze-057/S-roseum-operon-promoter-analysis.git
```

Navigate to the project root directory:

```bash
cd S.-roseum-genome-analysis
```

### 1. Download the data

Run the following command from the project root:

```bash
bash download_data.sh
```

### 2. Operon prediction

To perform log-likelihood-based operon prediction using data from *Streptomyces coelicolor*, run:

```bash
bash scripts/run_likelihood_analysis.sh
```

It is recommended to execute this script before proceeding with the subsequent analyses.

Generated figures are available in `results/operon_prediction_ll/`.

### 3. Promoter prediction

To identify putative promoters using MEME and sequences derived from *Streptomyces lividans*, run:

```bash
bash scripts/run_promoter_analysis.sh
```

FIMO and MEME results are available in `results/promoter_prediction/`.

### 4. Promoter–operon pair identification

To identify promoter–operon pairs and generate a CSV file containing the results, run:

```bash
python3 scripts/operon_prediction/operon_prediction.py
```

The resulting table is saved as `results/operons_with_promoters.csv`.

## Running on Windows (PowerShell)

First, clone the GitHub repository:

```powershell
git clone https://github.com/haze-057/S-roseum-operon-promoter-analysis.git
cd S.-roseum-genome-analysis
```

### 1. Download the data

Run the following command from the project root:

```powershell
.\scripts\download_data.ps1
```

### 2. Operon prediction

To perform log-likelihood-based operon prediction using data from *Streptomyces coelicolor*, run:

```powershell
.\scripts\run_likelihood_analysis.ps1
```

It is recommended to execute this script before proceeding with the subsequent analyses.

Generated figures are available in `results/operon_prediction_ll/`.

### 3. Promoter prediction

To identify putative promoters using MEME and sequences derived from *Streptomyces lividans*, run:

```powershell
.\scripts\run_promoter_analysis.ps1
```

FIMO and MEME results are available in `results/promoter_prediction/`.

### 4. Promoter–operon pair identification

To identify promoter–operon pairs and generate a CSV file containing the results, run:

```powershell
python scripts/operon_prediction/operon_prediction.py
```

The resulting table is saved as `results/operons_with_promoters.csv`.

## Additional output files

Additional intermediate files generated during the analyses are available in the `results/` directory.
