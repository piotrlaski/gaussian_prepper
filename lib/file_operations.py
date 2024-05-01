from lib.utils import *

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

def create_QMMM_input(outputdir, name, crystal, state, base, functional, molecule = None, nstates = 20, MM_theory = 'UFF', additional = '', genecp = ''):
    qmmm_file = os.path.join(outputdir, f'{name}_{state}_qmmm.inp')
    # prepare necessary headline strings
    if state[1] == '0' or state == 'T1' or state == 't1':
        headline_1 = f'#N ONIOM({functional}/{base}:{MM_theory})=EmbedCharge Opt Int=UltraFine SCF=XQC {additional}\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_1 = f'#N ONIOM(TD(triplets,nstates={nstates},Root={state[1]}) {functional}/{base}:{MM_theory})=EmbedCharge Opt Int=UltraFine SCF=XQC {additional}\n'
    else:
        headline_1 = f'#N ONIOM(TD(nstates={nstates},Root={state[1]}) {functional}/{base}:{MM_theory})=EmbedCharge Opt Int=UltraFine SCF=XQC {additional}\n'
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

def create_charge_calc_input(outputdir, name, molecule, state, base, functional, nstates = None, crystal = None, MM_theory = None, additional = '', genecp = ''):
    charge_calc_file = os.path.join(outputdir, f'{name}_{state}_charge.inp')
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

def create_iso_opt_input(outputdir, name, molecule, state, base, functional, nstates = None, crystal = None, MM_theory = None, additional = '', genecp = ''):
    isolated_opt_file = os.path.join(outputdir, f'{name}_{state}_opt.inp')
    if state[1] == '0' or state == 'T1' or state == 't1':
        headline_1 = f'#N {functional}/{base} Opt Int=UltraFine SCF=XQC {additional}\n'
    elif state[0] == 't' or state[0] == 'T':
        headline_1 = f'#N {functional}/{base} TD=(triplets, NStates={nstates},Root={state[1]}) Opt Int=UltraFine SCF=XQC {additional}\n'
    else:
        headline_1 = f'#N {functional}/{base} TD=(NStates={nstates},Root={state[1]}) Opt Int=UltraFine SCF=XQC {additional}\n'
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

def create_tddft_input(outputdir, name, molecule, state, base, functional, nstates = None, crystal = None, MM_theory = None, additional = '', genecp = ''):
    tddft_file = os.path.join(outputdir, f'{name}_{state}_tddft.inp')
    headline_1 = f'#N {functional}/{base} TD=(nstates={nstates}) pop=(full,Hirshfeld) Int=UltraFine {additional}\n'
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

def create_counterpoise_input(outputdir, name, crystal, state, base, functional, nstates = None, molecule = None, MM_theory = None, additional = '', genecp = ''):
    counterpoise_file = os.path.join(outputdir, f'{name}_{state}_counterpoise.inp')
    headline_1 = f'#N {functional}/{base} Counterpoise=2 {additional}\n'
    headline_2 = '0 1 0 1 0 1\n'

    with open(counterpoise_file, '+w') as f:
        f.write(headline_1)
        f.write('\n')
        f.write('Counterpoise Calculation\n')
        f.write('\n')
        f.write(headline_2)
        for i, mate_molecule in enumerate(crystal.get_molecules()):
            for atom in mate_molecule.get_atoms():
                f.write('{:<18}{:>12}{:>12}{:>12}\n'.format(f'{atom.name}(Fragment={i + 1})', f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))   
        f.write('\n')
        f.write(genecp)

def multiple_inputs(input_creation_function, outputdir, name, state, base, functional, crystal = None, molecule = None, nstates = 20, MM_theory = 'UFF', additional = '', genecp = []):
    if genecp:
        it_genecp = iter(genecp)
    for cur_base, cur_functional, cur_name in zip(base, functional, name):
        if cur_base.lower() == 'genecp':
            cur_genecp = next(it_genecp)
        else:
            cur_genecp = ''
        for cur_state in state:
            input_creation_function(outputdir = outputdir, name = cur_name, crystal = crystal, molecule = molecule, state = cur_state, base = cur_base, functional = cur_functional, nstates = nstates, MM_theory = MM_theory, additional = additional, genecp = cur_genecp)

def save_xyz(new_xyz_file, data_object):
    with open (new_xyz_file, '+w') as f:
        if type(data_object) is Crystal:
            for molecule in data_object.get_molecules():
                for atom in molecule.get_atoms():
                    f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))
        elif type(data_object) is Molecule:
            for atom in data_object.get_atoms():
                    f.write('{:<7}{:>12}{:>12}{:>12}\n'.format(atom.name, f'{atom.x_ort:.5f}', f'{atom.y_ort:.5f}', f'{atom.z_ort:.5f}'))