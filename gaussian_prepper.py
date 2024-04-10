import numpy as np
import itertools
import re

def apply_symmetry(old_fracs, sym_op):
    x, y, z = old_fracs
    new_x_frac, new_y_frac, new_z_frac = eval(sym_op[0]), eval(sym_op[1]), eval(sym_op[2])
    return ([new_x_frac, new_y_frac, new_z_frac])

def read_xyz_file(xyz_file):
    molecule_ort = []
    with open (xyz_file, '+r') as f:
        while line := f.readline():
            if len(line.split()) == 4:  #[Atom x y z] line
                atomic_line = line.split()
                atomic_line.append(0)  # adding dud charge
                molecule_ort.append(atomic_line)
    return (molecule_ort)

def read_charge_log_file(log_file):
    with open (log_file, '+r') as f:
        xyz_regex = r'(?<=\sSymbolic\sZ-matrix:\n\sCharge\s=\s\s.\sMultiplicity\s=\s.\n)[\s\S]*(?=Stoichiometry)'
        xyz = re.findall(xyz_regex, f.read())[0].split('\n')
        xyz = [atom.split() for atom in xyz][:-2]   #removing 2 empty \n lines at the end
        f.seek(0)
        chrg_regex = r"(?<=Q-CM5\s\s\s\n)[\s\S]*(?=\n[\s]{7}Tot)"
        chrg = re.findall(chrg_regex, f.read())[0].split('\n')
        chrg = [atom.split()[2] for atom in chrg]   #taking only the charge from matched lines
    for i in range(len(chrg)):
        xyz[i].append(chrg[i])
    return (xyz)

def read_cif(cif_file):
    with open (cif_file, '+r') as f:
        sym_ops_regex = r'(?<=_space_group_symop_operation_xyz\n)[\s\S]*(?=\n\n_cell_length_a)'
        sym_ops = re.findall(sym_ops_regex, f.read())[0].split(r"'")[1:-1]
        sym_ops = list(set(sym_ops))
        sym_ops.remove('\n ')
        sym_ops = [i.split(', ') for i in sym_ops]
        f.seek(0)
        cell_params_regex = r'(?<=cell_length_a)[\s\S]*(?=\n_cell_volume)' 
        cell_params = re.findall(cell_params_regex, f.read())[0]
        cell_params = re.findall(r'(?<=\s)\d+[\.]*\d*', cell_params)
        cell_params = [float(i) for i in cell_params]
    return (sym_ops, cell_params)

def create_QMMM_input(qmmm_file, crystal, state = 'S0', base = 'lanl2dz', functional = 'PBE1PBE', nstates = 20, MM_theory = 'UFF', additional = '', genecp = ''):
    # prepare necessary headline strings
    if state[1] == '0' or state == 'T1' or state == 't1':
        headline_1 = f'#N ONIOM({functional}/{base}:{MM_theory})=EmbedCharge Opt Freq Int=UltraFine SCF=XQC {additional}\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_1 = f'#N ONIOM(TD(triplets,nstates={nstates},Root={state[1]}) {functional}/{base}:{MM_theory})=EmbedCharge Opt Freq Int=UltraFine SCF=XQC {additional}\n'
    else:
        headline_1 = f'#N ONIOM(TD(nstates={nstates},Root={state[1]}) {functional}/{base}:{MM_theory})=EmbedCharge Opt Freq Int=UltraFine SCF=XQC {additional}\n'
    if state[0] == 's' or state[0] == 'S':
        headline_2 = '0 1 0 1 0 1\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_2 = '0 1 0 3 0 3\n'
    # find the central molecule for QM level
    for molecule in crystal.get_molecules():
        if molecule.sym_id == ['x+0', 'y+0', 'z+0']:
            central_molecule = molecule
            break
    with open(qmmm_file, '+w') as f:
        f.write(headline_1)
        f.write('\n')
        f.write(f'{state} QMMM Optimization\n')
        f.write('\n')
        f.write(headline_2)
        for atom in central_molecule.get_atoms():
            f.write('{:<17}0{:>16}{:>15}{:>15} H\n'.format(f'{atom.name}--{atom.charge:.5f}', f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))
        for molecule in crystal.get_molecules():
            if molecule is central_molecule:
                continue
            for atom in molecule.get_atoms():
                f.write('{:<16}-1{:>16}{:>15}{:>15} L\n'.format(f'{atom.name}--{atom.charge:.5f}', f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))      
        f.write('\n')
        f.write(genecp)

def create_charge_calc_input(charge_calc_file, molecule, state = 'S0', base = 'lanl2dz', functional = 'PBE1PBE', additional = '', genecp = ''):

    headline_1 = f'#N {functional}/{base} pop=(full,Hirshfeld) Int=UltraFine {additional}\n'

    if state[0] == 's' or state[0] == 'S':
        headline_2 = '0 1\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_2 = '0 3\n'

    with open(charge_calc_file, '+w') as f:
        f.write(headline_1)
        f.write('\n')
        f.write('Charge Calculation\n')
        f.write('\n')
        f.write(headline_2)
        for atom in molecule.get_atoms():
            f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))   
        f.write('\n')
        f.write(genecp)

def create_iso_opt_input(isolated_opt_file, molecule, state = 'S0', base = 'lanl2dz', functional = 'PBE1PBE', nstates = 20, additional = '', genecp = ''):
    
    if state[1] == '0' or state == 'T1' or state == 't1':
        headline_1 = f'#N {functional}/{base} Opt Freq Int=UltraFine SCF=XQC {additional}\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_1 = f'#N {functional}/{base} TD=(triplets, NStates={nstates},Root={state[1]}) Opt Freq Int=UltraFine SCF=XQC {additional}\n'
    else:
        headline_1 = f'#N {functional}/{base} TD=(NStates={nstates},Root={state[1]}) Opt Freq Int=UltraFine SCF=XQC {additional}\n'
    if state[0] == 's' or state[0] == 'S':
        headline_2 = '0 1\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_2 = '0 3\n'

    with open(isolated_opt_file, '+w') as f:
        f.write(headline_1)
        f.write('\n')
        f.write(f'{state} Isolated Molecule Optimization\n')
        f.write('\n')
        f.write(headline_2)
        for atom in molecule.get_atoms():
            f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))   
        f.write('\n')
        f.write(genecp)

def create_tddft_input(tddft_file, molecule, state = 'S0', base = 'lanl2dz', functional = 'PBE1PBE', nstates = 20, additional = '', genecp = ''):

    headline_1 = f'#N {functional}/{base} TD=(nstates={nstates}, 50-50) pop=(full,Hirshfeld) Int=UltraFine {additional}\n'

    if state[0] == 's' or state[0] == 'S':
        headline_2 = '0 1\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_2 = '0 3\n'

    with open(tddft_file, '+w') as f:
        f.write(headline_1)
        f.write('\n')
        f.write('TDDFT Calculation\n')
        f.write('\n')
        f.write(headline_2)
        for atom in molecule.get_atoms():
            f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))   
        f.write('\n')
        f.write(genecp)

def create_counterpoise_input(counterpoise_file, crystal, state = 'S0', base = 'lanl2dz', functional = 'PBE1PBE', nstates = 20, additional = '', genecp = ''):
    
    headline_1 = f'#N {functional}/{base} Counterpoise=2 {additional}\n'
    headline_2 = '0 1 0 1 0 1\n'

    with open(counterpoise_file, '+w') as f:
        f.write(headline_1)
        f.write('\n')
        f.write('Counterpoise Calculation\n')
        f.write('\n')
        f.write(headline_2)
        for i, molecule in enumerate(crystal.get_molecules()):
            for atom in molecule.get_atoms():
                f.write('{:<18}{:>12}{:>12}{:>12}\n'.format(f'{atom.name}(Fragment={i + 1})', f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))   
        f.write('\n')
        f.write(genecp)

def save_xyz(new_xyz_file, data_object):
    with open (new_xyz_file, '+w') as f:
        if type(data_object) is Crystal:
            for molecule in data_object.get_molecules():
                for atom in molecule.get_atoms():
                    f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))
        elif type(data_object) is Molecule:
            for atom in data_object.get_atoms():
                    f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))

def dist(v1, v2):
    return np.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2)

def unit_cell_vectors(cell_parameters):
    # Convert angles to radians
    alpha, beta, gamma = np.radians(cell_parameters[3:])
    # Calculate cell volume
    a, b, c = cell_parameters[:3]
    volume = a * b * c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    # Calculate vectors
    a_vector = np.array([a, 0, 0])
    b_vector = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    c_vector = np.array([c * np.cos(beta),
                         (c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma))) / np.sin(gamma),
                         volume / (a * b * np.sin(gamma))])
    return a_vector, b_vector, c_vector

def ort2frac(xyz_coords, A):
    return (list(np.matmul(xyz_coords, np.linalg.inv(A))))

def frac2ort(frac_coords, A):
    return (list(np.matmul(frac_coords, A)))

class Crystal():
    def __init__(self, cell_params, sym_ops) -> None:
        self.cell_params = cell_params
        self.sym_ops = sym_ops
        self.cell_vec_a, self.cell_vec_b, self.cell_vec_c = unit_cell_vectors(cell_params)
        self.ort_matrix = np.array([self.cell_vec_a, self.cell_vec_b, self.cell_vec_c])
        self.__molecules = []
        self.__unique_atom_set = set()
    def add_molecule(self, molecule):
        self.__molecules.append(molecule)
    def get_molecules(self):
        return self.__molecules
    def del_molecule(self, molecule):
        self.__molecules.remove(molecule)
    def spawn_sym_mate(self, main_molecule, sym_op):
        new_molecule = Molecule(crystal = self, sym_id = sym_op)
        for main_atom in main_molecule.get_atoms():
            main_atom_fracs = [main_atom.x_frac, main_atom.y_frac, main_atom.z_frac]
            new_fracs = apply_symmetry(main_atom_fracs, sym_op)
            new_atom_unique_fracs = tuple(round(i, 2) for i in new_fracs)
            if new_atom_unique_fracs not in self.__unique_atom_set:
                Atom(molecule = new_molecule, name = main_atom.name, fracs = new_fracs, charge = main_atom.charge)
    def build_infinite_crystal(self, main_molecule, p = 3):
        for sym_op in self.sym_ops:
            for mod_x, mod_y, mod_z in itertools.product(range(-p, p+1), repeat=3):
                new_sym_id = [f'{sym_op[0]}+{str(mod_x)}', f'{sym_op[1]}+{str(mod_y)}', f'{sym_op[2]}+{str(mod_z)}']
                self.spawn_sym_mate(main_molecule, sym_op = new_sym_id)
    def cut_out_cluster(self, radius):
        deletion_list = []
        for molecule in self.get_molecules():
            for atom in molecule.get_atoms():
                if dist([0, 0, 0], [atom.x_ort, atom.y_ort, atom.z_ort]) < radius:
                    break
            else:
                deletion_list.append(molecule)
        self.__molecules = list(set(self.get_molecules()).difference(deletion_list))
    def add_unique_atom(self, atom):
        self.__unique_atom_set.add(atom)
    def get_unique_atom_set(self):
        return self.__unique_atom_set
    def clear_unique_atom_set(self):
        self.__unique_atom_set = set()



class Molecule(Crystal):
    def __init__(self, crystal, sym_id, add_to_crystal = True) -> None:
        super().__init__(cell_params = crystal.cell_params, sym_ops = crystal.sym_ops)
        self.sym_id = sym_id
        self.__atoms = []
        self.parent_crystal = crystal
        if add_to_crystal:
            crystal.add_molecule(self)
    def add_atom(self, atom):
        self.__atoms.append(atom)
    def get_atoms(self):
        return self.__atoms
    def del_atom(self, atom):
        self.__atoms.remove(atom)
    def find_nearest_non_H_atom(self, h_atom): #returns the h_atom if no atoms in 1.2 range
        lowest_dist = 1.2
        nearest_atom = h_atom
        for atom in self.__atoms:
            atom_dist = dist([h_atom.x_ort, h_atom.y_ort, h_atom.z_ort], [atom.x_ort, atom.y_ort, atom.z_ort])
            if atom_dist < lowest_dist and atom.name != 'H':
                nearest_atom = atom
                lowest_dist = atom_dist
        return nearest_atom       
    def extend_hydrogen_bonds(self, C_H_BOND = 1.089, N_H_BOND = 1.015, O_H_BOND = 0.993):
        for atom in self.__atoms:
            if atom.name == 'H':
                nearest_atom = self.find_nearest_non_H_atom(atom)
                if nearest_atom.name == 'C':
                    new_dist = C_H_BOND
                elif nearest_atom.name == 'N':
                    new_dist = N_H_BOND
                elif nearest_atom.name == 'O':
                    new_dist = O_H_BOND
                else: continue
                h_orts = [atom.x_ort, atom.y_ort, atom.z_ort]
                nearest_orts = [nearest_atom.x_ort, nearest_atom.y_ort, nearest_atom.z_ort]
                x_to_h = np.subtract(h_orts, nearest_orts)
                norm_x_to_h = x_to_h / np.sqrt(np.sum(x_to_h ** 2))
                x_to_h_extension = norm_x_to_h * new_dist
                h_ext = np.add(nearest_orts, x_to_h_extension)
                atom.update_ort(h_ext)

class Atom(Molecule):
    def __init__(self, molecule, name = None, orts = None, fracs = None, charge = None, add_to_molecule = True, add_to_unique_atom_set = True) -> None:
        super().__init__(molecule.parent_crystal, sym_id = molecule.sym_id, add_to_crystal = False)
        self.name = name
        self.parent_molecule = molecule
        if orts is not None and fracs is None:
            self.x_ort, self.y_ort, self.z_ort = orts
            fracs = ort2frac(orts, self.ort_matrix)
            self.x_frac, self.y_frac, self.z_frac = fracs
        elif orts is  None and fracs is not None:
            self.x_frac, self.y_frac, self.z_frac = fracs
            orts = frac2ort(fracs, self.ort_matrix)
            self.x_ort, self.y_ort, self.z_ort = orts       
        else:
            self.x_ort, self.y_ort, self.z_ort = None, None, None
            self.x_frac, self.y_frac, self.z_frac = None, None, None
        self.charge = charge
        if add_to_molecule:
            molecule.add_atom(self)
        if add_to_unique_atom_set:
            molecule.parent_crystal.add_unique_atom((round(self.x_frac, 2), round(self.y_frac, 2), round(self.z_frac, 2)))
    def update_ort(self, new_orts):
        self.x_ort, self.y_ort, self.z_ort = new_orts
        new_fracs = ort2frac(new_orts, self.ort_matrix)
        self.x_frac, self.y_frac, self.z_frac = new_fracs
    def update_frac(self, new_fracs):
        self.x_frac, self.y_frac, self.z_frac = new_fracs
        new_orts = frac2ort(new_fracs, self.ort_matrix)
        self.x_ort, self.y_ort, self.z_ort = new_orts

sym_ops, cell_params = read_cif(r'C:\Users\piotr\Documents\VS_Code\cluster\Rh(4-Br-SA)(CO)2__Q1_21Dlk1__CCDC.cif')
xyz_chrgd_ort = read_charge_log_file(r'C:\Users\piotr\Documents\VS_Code\cluster\chrg.log')
molecule_ort = read_xyz_file(r'C:\Users\piotr\Documents\VS_Code\cluster\nonextend.xyz')


CIF_PATH = r'C:\Users\piotr\Documents\VS_Code\cluster\Rh(4-Br-SA)(CO)2__Q1_21Dlk1__CCDC.cif'
XYZ_PATH = r'C:\Users\piotr\Documents\VS_Code\cluster\Rh_4-Br-SA_CO_2.xyz'

EXTEND_HYDROGENS = True

CREATE_CHARGE_CALC = True
CREATE_ISOLATED_OPT = True
CREATE_TDDFT_CALC = True

CREATE_COUNTERPOISE = True
COUNTERPOISE_SYMMETRY = (['-x+3', '-y+1', '-z+2'], ['x+0', 'y+0', 'z+0'])

CREATE_ONIOM = False
CHARGE_LOG_PATH = r'C:\Users\piotr\Documents\VS_Code\cluster\chrg.log'

NAME = 'cam_6311_lan'

FUNCTIONAL = 'CAM-B3LYP'
BASE = 'GenECP'
STATE = 'T1'
NSTATES = 20
RADIUS = 15
ADDITIONAL = 'EmpiricalDispersion=GD3BJ'

GENECP_FILE = r'C:\Users\piotr\Documents\VS_Code\cluster\gen_ecp.txt'

with open(GENECP_FILE, '+r') as f:
    GENECP = f.read()

if CREATE_ONIOM:
    loaded_data = xyz_chrgd_ort
else:
    loaded_data = molecule_ort

xtal = Crystal(cell_params, sym_ops)
main_molecule = Molecule(crystal = xtal,
                        sym_id = ['x+0', 'y+0', 'z+0'],
                        add_to_crystal = True)


for atomic_line in loaded_data:
    Atom(molecule = main_molecule,
         name = atomic_line[0],
         orts = [float(coord) for coord in atomic_line[1:4]],
         charge = float(atomic_line[4]),
         add_to_unique_atom_set = True)

if EXTEND_HYDROGENS:
    main_molecule.extend_hydrogen_bonds()

if CREATE_COUNTERPOISE:
    xtal.del_molecule(main_molecule)
    xtal.clear_unique_atom_set()
    xtal.spawn_sym_mate(main_molecule, COUNTERPOISE_SYMMETRY[0])
    xtal.spawn_sym_mate(main_molecule, COUNTERPOISE_SYMMETRY[1])
    save_xyz(f'counterpoise_dimer.xyz', xtal)
    create_counterpoise_input(f'{NAME}_counterpoise.inp', xtal, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = GENECP)

if CREATE_ONIOM:
    xtal.build_infinite_crystal(main_molecule)
    xtal.cut_out_cluster(radius = RADIUS)
    save_xyz(f'{NAME}_{STATE}_cluster.xyz', xtal)
    create_QMMM_input(f'{NAME}_{STATE}_qmmm.inp',xtal, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = GENECP)

if CREATE_CHARGE_CALC:
    create_charge_calc_input(f'{NAME}_{STATE}_chrg.inp', main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, additional = ADDITIONAL, genecp = GENECP)

if CREATE_ISOLATED_OPT:
    create_iso_opt_input(f'{NAME}_{STATE}_isolated.inp', main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = GENECP)

if CREATE_TDDFT_CALC:
    create_tddft_input(f'{NAME}_{STATE}_tddft.inp', main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = GENECP)

