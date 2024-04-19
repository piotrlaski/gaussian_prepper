
# full paths .cif containing symmetry data and .xyz containing inquired molecular geometry
CIF_PATH = r'C:\Users\piotr\Documents\VS_Code\cluster\example_files\Rh(4-Br-SA)(CO)2__Q1_21Dlk1__CCDC.cif'
XYZ_PATH = r'C:\Users\piotr\Documents\VS_Code\cluster\example_files\Rh_4-Br-SA_CO_2.xyz'

# hydrogen extension for CH, OH, NH bonds (True/False)
EXTEND_HYDROGENS = True

# which Gaussian outputs to create (True/False)
CREATE_CHARGE_CALC = True
CREATE_ISOLATED_OPT = True
CREATE_TDDFT_CALC = True

# interaction energy calculation input creation (True/False)
CREATE_COUNTERPOISE = True
# provide symmetry identifiers eg. (['-x+3', '-y+1', '-z+2'], ['x+0', 'y+0', 'z+0']) of the 2 molecules (tuple(list[str,str,str], list [str,str,str]))
COUNTERPOISE_SYMMETRY = (['-x+3', '-y+1', '-z+2'], ['x+0', 'y+0', 'z+0'])

# QMMM calculation input creation (True/False)
CREATE_ONIOM = True
# post-Gaussian charge calculation *.log file full path
CHARGE_LOG_PATH = r'C:\Users\piotr\Documents\VS_Code\cluster\example_files\chrg.log'
# molecular cluster radius in Angstroms (integer)
RADIUS = 15

# Name of the Gaussian job
NAME = 'cam_6311_lan'

# Gaussian alculation parameters
FUNCTIONAL = 'CAM-B3LYP'
BASE = 'GenECP'
STATE = 'T1'

#Number of states for TD calculations (also when calculating higher excited states)
NSTATES = 20
# Additional keywords to be put at the end of header in inputs
ADDITIONAL = 'EmpiricalDispersion=GD3BJ'

# only for GenECP custom base/pseudopotential - the appendix you would normally put after molecule geometry in a Gaussian input (text file full path)
GENECP_FILE = r'C:\Users\piotr\Documents\VS_Code\cluster\example_files\gsen_ecp.txt'