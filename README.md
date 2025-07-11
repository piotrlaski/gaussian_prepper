# Gaussian Prepper

Gaussian Prepper is a tool for producing Gaussian16 input files without all the hassle of exact whitespace-sensitive syntax surgery required to prepare proper input files. Useful for doing broad searches of the basis set/funtional theory level which is most suited for your molecular system. Pair it with QM Distance Calculator to obtain reports on selected geometry parameters among all of the resulting calculations.

It is the public version of the project "CM_Batch", which is referenced in my thesis (chapter 2.4.3.1 - p.63-64). 

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

- **Functional and basis set**: Configure `FUNCTIONAL` and `BASE` lists - remeber to put a single functional keyword for a single base keyword, as they will be taken one-by-one
- **Electronic states**: Specify the list of states like S0, S1, T1 etc. in the `STATE` parameter
- **Custom pseudopotentials**: Use `GENECP_FILE` path parameter for custom basis sets and ECP pseudopotentials (a file has to be provided containing the basis/ECP specification, which is the same as you would put at the end of a Gaussian input file, i.e. all the characters past the molecular geometry)
- **QM/MM settings**: Configure cluster radius and other ONIOM parameters
- **Counterpoise settings**: Specify the symmetries of the desired counterpoise twins. The central molecule (x,y,z) is always present, another will be created for each provided symmetry operation.

## Features

- **Hydrogen bond standardization**: Option to standardize H-X bonds (C-H, O-H, N-H) to the neutron-standarized values for your uploaded geometry
- **Custom basis sets**: Specify your own basis set and/or ECP pseudopotential
- **Multiple job creation**: Jobs are created for every specified electronic state using the corresponding functional/base/name from the lists (all combinations are considered)
- **Automatic shell script generation**: Creates submission scripts with optimized computational parameters for different job types, which can also be manually specified. This is mostly for the SLURM queuing system installed on WCSS network.

The script also enables the creation of complex QM/MM input files using the ONIOM method. This functionality requires an additional log file from a prior Hirshfeld charge calculation. The script then:
1.	Constructs a cluster of a user-specified radius (default 15 Å) around the central molecule.
2.	Assigns the appropriate partial charges to all atoms within the cluster.
3.	Partitions the cluster into high-level (quantum) and low-level (molecular mechanics) regions.
4.	Automatically adjusts the total charge and multiplicity declarations to match the requested electronic state.


## File Structure

- `gaussian_prepper.py` - Main execution file
- `config.py` - Configuration file for all input parameters
- `lib/utils.py` - Core calculation and structure handling functions
- `lib/file_operations.py` - Input file generation functions

## Usage Notes

The pseudopotential GenECP files will be accessed in order, whenever a GenECP base is specified. 

For any additional help or bug reports, please contact: **pa.laski@uw.edu.pl**
