# Gaussian Prepper

Gaussian Prepper is a tool for producing Gaussian16 input files without all the hassle of exact whitespace-sensitive syntax surgery required to prepare proper input files.

## Getting Started

Clone the repository into your local space (or download all .py files into an empty directory) and configure your input in the `config.py` file. When you are ready, run `gaussian_prepper.py` to generate input files in the selected output directory

## Required Input Files

You will need:

1. **A .cif file** containing necessary structural information about your molecule (EXCEPT geometry!)
2. **An .xyz file** containing the geometry of your molecule
3. **For QM/MM geometry calculations**: charge distribution information provided as a raw Gaussian .log file from a prior charge calculation job

## Supported Job Types

Currently you can create inputs for:

1. **Geometry optimization jobs**
2. **TD-DFT jobs**
3. **Charge calculation jobs** (output of which can be used for QM/MM)
4. **QM/MM jobs**
5. **Interaction energy jobs (counterpoise)**

## Configuration

You can specify the following parameters in `config.py`:

- **Functional and basis set**: Configure `FUNCTIONAL` and `BASE` lists
- **Electronic states**: Specify states like S0, S1, T1 etc. in the `STATE` parameter
- **Custom pseudopotentials**: Use `GENECP_FILE` for custom basis sets and ECP pseudopotentials (a file has to be provided containing the basis/ECP specification)
- **QM/MM settings**: Configure cluster radius and other ONIOM parameters

## Features

- **Hydrogen bond standardization**: Option to standardize H-X bonds (C-H, O-H, N-H) for your uploaded geometry
- **Custom basis sets**: Specify your own basis set and/or ECP pseudopotential
- **Multiple job creation**: Jobs are created for every specified electronic state using the corresponding functional/base/name from the lists (all combinations are considered)
- **Automatic shell script generation**: Creates submission scripts with optimized computational parameters for different job types, which can also be manually specified

## File Structure

- `gaussian_prepper.py` - Main execution file
- `config.py` - Configuration file for all input parameters
- `lib/utils.py` - Core calculation and structure handling functions
- `lib/file_operations.py` - Input file generation functions

## Usage Notes

The pseudopotential GenECP files will be accessed in order, whenever a GenECP base is specified. The tool automatically handles the complex syntax requirements of Gaussian input files and generates properly formatted inputs for your specified calculations.

For any additional help or bug reports, please contact: **pa.laski@uw.edu.pl**
