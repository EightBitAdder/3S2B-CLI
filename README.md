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
Run the main.py script in order to start using the tool.
```bash
python main.py
```

## Using the Tool

After running the main.py script, five menu options will be displayed to the console. They are as follows:

1. **View Index Table**: Pressing **V** (or **v**), followed by **enter** on your keyboard will display the MFD Index Table. This is the *master table* which pairs each of the SWGDrug molecules with their associated Fragment Table (*see below*).

2. **Search and Fetch**: Pressing **F** (or **f**), followed by **enter** on your keyboard will prompt the user to select an entry in the MFD Index Table to view. In particular, the user may select an entry either by *SMILES*; or, by *CRAFTS Lab Entry Number*. To select an entry by *CRAFTS Lab Entry Number*, simply type **CL** followed by the associated entry number in the MFD Index Table.

3. **Download**: Pressing **D** (or **d**), followed by **enter** on your keyboard will prompt the user to select an entry in the MFD Index Table to download. In particular, the user may select an entry either by *SMILES*; or, by *CRAFTS Lab Entry Number*. To select an entry by *CRAFTS Lab Entry Number*, simply type **CL** followed by the associated entry number in the MDF Index Table.

4. **Compare**: Pressing **C** (or **c**), followed by enter on your keyboard will prompt the user to enter in two file paths, separated by the following string: **" "**. The first file path must refer to a .csv or .txt file containing mass-spectral data (m/z and intensity observations), separated by a whitespace (*This is currently in the process of being updated to allow for arbitrary separator types. Stay tuned...*). The second file path must refer to a .csv or .txt file containing a list of masses. The user will also be prompted to enter a tolerance about which to calculate the FPIE score. Upon hitting enter, the two files will be compared against one another. Then an FPIE score together with a stem plot will be produced.

5. **Compare All**: Pressing **A** (or **a**), followed by enter on your keyboard will prompt the user to enter in one file path. The file path must refer to a .csv or .txt file containing mass-spectral data (m/z and intensity observations), separated by a whitespace (*This is currently in the process of being updated to allow for arbitrary separator types. Stay tuned...*). The user will also be prompted to enter a tolerance about which to calculate the FPIE score. Upon hitting enter, the file will be compared against the entire MFD. Then a table containing all of the calculated FPIE scores will be produced.

## Additional Notes

We are actively working on providing support for Linux and MacOS. Stay tuned for updates!