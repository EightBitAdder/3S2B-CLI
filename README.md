# 3S2B CLI

Welcome to the 3S2B CLI. This is a command-line interface tool for the 3S2B algorithm and the MFD.

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

2. Run the setup.bat script in order to configure your environment on Windows:
    ```bash
    setup.bat
    ```

    *(A .sh script for Linux and MacOS is actively being working on. Stay tuned...)*

### Step 2: Use the Tool
Run the cli.bat script in order to start using the tool.
```bash
cli.bat
```

## Using the Tool

After running the cli.bat script, five menu options will be displayed to the console. They are as follows:

1. **View Index Table**: Pressing **i** (followed by **enter**) will display the MFD Index Table. This is the *master table* which pairs each of the SWGDrug molecules with their associated Fragment Table.

2. **Search and Fetch**: Pressing **f**, followed an entry in the MFD Index Table (followed by **enter**) will display the Fragment Table of the corresponding entry. In particular, the user may select an entry either by *SMILES*; or, by *CRAFTS Lab Entry Number*. To select an entry by *CRAFTS Lab Entry Number*, simply type **CL** followed by the associated entry number in the MFD Index Table.

3. **Compare**: Pressing **c**, followed by two file paths, as well as a tolerance about which to calculate the FPIE score (followed by **enter**) will cause the two files to be compared against one another. Then the FPIE score together with a stem plot will be displayed. *Note*: The first file path must refer to a .csv or .txt file containing mass-spectral data (m/z and intensity observations), separated by a whitespace (*This is currently in the process of being updated to allow for arbitrary separator types. Stay tuned...*). The second file path must refer to a .csv or .txt file containing a list of masses.

4. **Compare All**: Pressing **a**, followed a file path, as well as a tolerance about which to calculate the FPIE scores (followed by **enter**) will cause the file to be compared against the entire MFD. Then a table containing all of the calculated FPIE scores will be displayed. *Note*: The file path must refer to a .csv or .txt file containing mass-spectral data (m/z and intensity observations), separated by a whitespace (*This is currently in the process of being updated to allow for arbitrary separator types. Stay tuned...*).

5. **Fragment**: Pressing **fr**, followed by a *SMILES* string (followed by **enter**) will display the Fragment Table of the associated *SMILES* string.

## Additional Notes

We are actively working on providing support for Linux and MacOS. Stay tuned for updates!