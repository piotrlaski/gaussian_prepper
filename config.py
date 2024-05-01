# full paths .cif containing symmetry data and .xyz containing inquired molecular geometry
CIF_PATH = r'C:\Users\piotr\Documents\working_dirs_lapek\gaussian_prepper\example_files\Rh(4-Br-SA)(CO)2__Q1_21Dlk1__CCDC.cif'
XYZ_PATH = r'C:\Users\piotr\Documents\working_dirs_lapek\gaussian_prepper\example_files\Rh_4-Br-SA_CO_2.xyz'
OUTPUT_DIRECTORY = r'C:\Users\piotr\Desktop\last_att_Rh\pure_try\testin'

# hydrogen extension for CH, OH, NH bonds (True/False)
EXTEND_HYDROGENS = True

CREATE_SHELL_SCRIPT = True
#Slurm computing params, CORES/MEMORY/TIME (-n, -m, -t):
CHARGE_PARAMS = (4, 30, 20)
ISOLATED_OPT_PARAMS = (12, 50, 100)
TDDFT_PARAMS = (12, 50, 100)
COUNTERPOISE_PARAMS = (4, 30, 20)
ONIOM_PARAMS = (24, 80, 150)

# which Gaussian outputs to create (True/False)
CREATE_CHARGE_CALC = False
CREATE_ISOLATED_OPT = False
CREATE_TDDFT_CALC = True

# interaction energy calculation input creation (True/False)
CREATE_COUNTERPOISE = False
# provide symmetry identifiers eg. (['-x+3', '-y+1', '-z+2'], ['x+0', 'y+0', 'z+0']) of the 2 molecules (tuple(list[str,str,str], list [str,str,str]))
COUNTERPOISE_SYMMETRY = (['-x+3', '-y+1', '-z+2'], ['x+0', 'y+0', 'z+0'])

# QMMM calculation input creation (True/False)
CREATE_ONIOM = False
# post-Gaussian charge calculation *.log file full path
CHARGE_LOG_PATH = [r'C:\Users\piotr\Desktop\last_att_Rh\pure_try\cam_6311_lan_S0_chrg.log',
                   r'C:\Users\piotr\Desktop\last_att_Rh\pure_try\cam_6311_lan_S0_chrgprim.log',
                   r'C:\Users\piotr\Desktop\last_att_Rh\pure_try\cam_6311_lan_S0_chrgbeta.log']
# molecular cluster radius in Angstroms (integer)
RADIUS = 10

# Name of the Gaussian job (it should ifentify the functional/base, different states will be named automatically)
NAME = ['b3lyp_6311g' ,'cam_genecp', 'pbe_lanl2dz']
FUNCTIONAL = ['B3LYP','CAM-B3LYP', 'PBE1PBE']
BASE = ['6-31G**' ,'GenECP', 'Lanl2dz']
# only for GenECP custom base/pseudopotential - the appendix you would normally put after molecule geometry in a Gaussian input (text file full path)
GENECP_FILE = [r'C:\Users\piotr\Documents\working_dirs_lapek\gaussian_prepper\example_files\gsen_ecp.txt']
STATE = ['T1','S0', 'S1']
# Additional keywords to be put at the end of header in inputs
ADDITIONAL = 'IOp(3/174=1000000) IOp(3/175=2067400) IOp(3/177=0370800) IOp(3/178=5474300)'
#Number of states for TD calculations (also when calculating higher excited states)
NSTATES = '20, singlets'    #you can also inject TD keywords like 'singlets' or '50/50' etc. here, e.g. NSTATES = '20, singlets'