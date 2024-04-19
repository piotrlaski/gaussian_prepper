Gaussian Prepper is a tool for producing Gaussian16 input files without all the hassle of exact whitespace-sensitive syntax surgery required to prepare proper input files.

Clone the repository into your local space (or download all .py files into an empty directory) and configure your input in the config.py file. When you are ready, run the gaussian_prepper.py to generate input files in the same directory.

You will need a .cif file containing necessary structural information about your molecule, as well as an .xyz containing the geometry of your molecule. For QM/MM geometry calculation, charge distribution information is required, which you can provide in form of a raw gaussian .log file from a prior charge calculation job.

Currently you can create inputs for:
1. Geometry optimization jobs
2. TD:DFT jobs
3. Charge calculation jobs (output of which can be used for QM/MM)
4. QM/MM jobs
5. Interaction energy jobs (counterpoise)

For any additional help or bug reports, please contact:
pa.laski@uw.edu.pl
