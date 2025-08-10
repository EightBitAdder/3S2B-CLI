# 3S2B CLI - Molecular Fragment Database Analysis Tool

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A comprehensive command-line interface and Python library for Fragment Pattern Intensity Explained (FPIE) scoring and molecular fragment database analysis, implementing the 3S2B algorithm.

## 📖 About

The 3S2B CLI provides tools for:
- **Fragment Analysis**: Systematic molecular fragmentation with mass calculation
- **FPIE Scoring**: Fragment Pattern Intensity Explained scoring against mass spectral data
- **Database Management**: Molecular fragment database creation and querying
- **Mass Spectrometry**: Peak matching and intensity analysis

For detailed methodology, see our paper: [3S2B Algorithm](https://chemrxiv.org/engage/chemrxiv/article-details/685aa01d1a8f9bdab548c2f4)

## 🚀 Installation

### Prerequisites

- Python 3.10 or higher
- RDKit (for molecular handling)
- IsoSpecPy (for isotope calculations)

### Install from Source

1. **Clone the repository:**
   ```bash
   git clone https://github.com/EightBitAdder/3S2B-CLI.git
   cd 3S2B-CLI
   ```

2. **Install the package:**
   ```bash
   pip install .
   ```

3. **Or install in development mode:**
   ```bash
   pip install -e ".[dev]"
   ```

### Install Dependencies

If you encounter issues with RDKit installation:

```bash
# Using conda (recommended for RDKit)
conda install -c conda-forge rdkit

# Or using pip (may require additional setup)
pip install rdkit
```

## 🎯 Quick Start

### 1. Create Database (First Time Setup)

```bash
# Create database from SDF file
3s2b make-db --sdf path/to/molecules.sdf

# Or let it auto-detect SDF location
export THREE_S_TWO_B_SDF_PATH="path/to/SWGDRUG313.SDF"
3s2b make-db
```

### 2. Fragment Analysis

```bash
# Generate fragments from SMILES
3s2b fragment "CCO" --max-mult 2

# Save fragments to file
3s2b fragment "c1ccccc1" --save benzene_fragments.csv
```

### 3. FPIE Scoring

```bash
# Compare MS data against fragment list
3s2b compare spectrum.csv fragments.csv 0.1

# Compare against entire database
3s2b compare-all spectrum.csv 0.1 --top 20

# With custom weight function and plotting
3s2b compare data.txt masses.txt 0.05 \
    --weight-function EXPONENTIAL-DECAY \
    --plot --save-plot results.png
```

### 4. Database Operations

```bash
# View database index
3s2b index --limit 50

# Search by SMILES, name, or ID
3s2b search "CCO"
3s2b search "CL123"

# Search by mass
3s2b mass-search 100.0 --tolerance 0.1

# Add molecule to database
3s2b add "CC(=O)O" --name "acetic acid"

# Database information
3s2b db-info
```

## Using the Tool

After running the cli.bat script, five menu options will be displayed to the console. They are as follows (*Note: type 3s2b " command " --help in the console for more information . . .*):

1. **View Index Table**: Pressing **i** (followed by **enter**) will display the MFD Index Table. This is the *master table* which pairs each of the SWGDrug molecules with their associated Fragment Table.

2. **Search and Fetch**: Pressing **f**, followed by an entry in the MFD Index Table (followed by **enter**) will display the Fragment Table of the corresponding entry. In particular, the user may select an entry either by *SMILES*; or, by *CRAFTS Lab Entry Number*. To select an entry by *CRAFTS Lab Entry Number*, simply type **CL** followed by the associated entry number in the MFD Index Table.

3. **Search and Fetch by Mass**: Pressing **m** followed by a mass value (*significant to two digits*) (followed by **enter**) will display a table containing all of the fragments in the MFD with the given mass value.

4. **Compare**: Pressing **c**, followed by two file paths, as well as a tolerance about which to calculate the FPIE score (followed by **enter**) will cause the two files to be compared against one another. Then the FPIE score will be displayed. *Note*: The first file path must refer to a .csv or .txt file containing mass-spectral data (m/z and intensity observations). The second file path must refer to a .csv or .txt file containing a list of masses.

    The **c** command also boasts the following options (*Note: type 3s2b c --help in the console for more information . . .*):

    - *-sr* selects the file reader to be used for the given mass-spectral data file. Currently, only the default **MS-DATA** file reader is available. Support for SDF files will be available soon.
    
    - *-lr* selects the file reader to be used for the given mass-list data file. Currently, only the default **MASS-LIST** file reader is available.
    
    - *-dr* specifies the delimiter used in both the mass-secptral and mass-list data files. If a delimiter is not specified, then the file reader will attempt to sniff the appropriate delimiter; otherwise, the file reader will default to using: *" \s+ "*.
    
    - *-wf* specifies the weight function to be used in the FPIE calculation. Currently only two weight functions are supported: **CONSTANT** (default when the given mass-list data file is not annotated with multiplicities) and **EXPONENTIAL-DECAY** (default when the given mass-list data file is annotated with multiplicities).
    
    - *--plot* specifies whether an annoted FPIE stem-plot should/should not be generated for the given comparison.

5. **Compare All**: Pressing **a**, followed a file path, as well as a tolerance about which to calculate the FPIE scores (followed by **enter**) will cause the file to be compared against the entire MFD. Then a table containing all of the calculated FPIE scores will be displayed. *Note*: The file path must refer to a .csv or .txt file containing mass-spectral data (m/z and intensity observations).

    The **a** command also boasts the following options (*Note: type 3s2b a --help in the console for more information . . .*):

    - *-sr* selects the file reader to be used for the given mass-spectral data file. Currently, only the default **MS-DATA** file reader is available. Support for SDF files will be available soon.
    
    - *-dr* specifies the delimiter used in both the mass-secptral and mass-list data files. If a delimiter is not specified, then the file reader will attempt to sniff the appropriate delimiter; otherwise, the file reader will default to using: *" \s+ "*.
    
    - *-wf* specifies the weight function to be used in the FPIE calculation. Currently only two weight functions are supported: **CONSTANT** (default when the given mass-list data file is not annotated with multiplicities) and **EXPONENTIAL-DECAY** (default when the given mass-list data file is annotated with multiplicities).

6. **Fragment**: Pressing **fr**, followed by a *SMILES* string (followed by **enter**) will display the Fragment Table of the associated *SMILES* string.

7. **Add**: Pressing **add**, followed by a *SMILES* string (followed by **enter**) inserts the given *SMILES* string into the MDF.

## Additional Notes

We are actively working on providing support for Linux and MacOS. Stay tuned for updates!