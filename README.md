# 3S2B CLI

Welcome to the 3S2B CLI. This is a command-line interface tool for the 3S2B algorithm and the MFD. Please visit https://chemrxiv.org/engage/chemrxiv/article-details/685aa01d1a8f9bdab548c2f4 to read our paper.

## Getting Started

To start using the tool, simply follow the steps listed below:

### Step 1: Set-up
*NOTE*: Python version 3.12.7 is required in order to use the 3S2B CLI.

1. Clone the repository to your local machine:
    ```bash
    git clone https://github.com/EightBitAdder/3S2B-CLI.git
    cd 3S2B-CLI
    ```

    Or else, simply download the repository as a .zip file onto your computer; making sure to extract the contents of the file before proceeding.

2. Run the bootstrap.py script in order to configure the project environment:
    ```bash
    python bootstrap.py
    ```

    and follow the on-screen instrcutions.

### Step 2: Use the Tool
Simply type:
```bash
3s2b
```

or

```bash
3s2b --help
```

into the command prompt where you will be greeted with a welcome message, as well as a list of all of the available commands.

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