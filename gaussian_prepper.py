from config import *
from lib.utils import *
from lib.file_operations import *

if __name__ == '__main__':

    if os.path.isfile(CIF_PATH) and os.path.isfile(CHARGE_LOG_PATH):
        loaded_data = read_charge_log_file(CHARGE_LOG_PATH)
    elif os.path.isfile(CIF_PATH) and os.path.isfile(XYZ_PATH):
        loaded_data = read_xyz_file(XYZ_PATH)
    else:
        print ('Invalid input files!')
        exit()

    sym_ops, cell_params = read_cif(CIF_PATH)

    if os.path.isfile(GENECP_FILE):
        with open(GENECP_FILE, '+r') as f:
            genecp = f.read()
    else:
        genecp = ''

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
        create_counterpoise_input(os.path.join(OUTPUT_DIRECTORY, f'{NAME}_counterpoise.inp'), xtal, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)

    if CREATE_ONIOM:
        xtal.build_infinite_crystal(main_molecule)
        xtal.cut_out_cluster(main_molecule=main_molecule, radius = RADIUS)
        save_xyz(os.path.join(OUTPUT_DIRECTORY, f'{NAME}_{STATE}_cluster.xyz'), xtal)
        create_QMMM_input(os.path.join(OUTPUT_DIRECTORY, f'{NAME}_{STATE}_qmmm.inp'),xtal, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)

    if CREATE_CHARGE_CALC:
        create_charge_calc_input(os.path.join(OUTPUT_DIRECTORY, f'{NAME}_{STATE}_chrg.inp'), main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, additional = ADDITIONAL, genecp = genecp)

    if CREATE_ISOLATED_OPT:
        create_iso_opt_input(os.path.join(OUTPUT_DIRECTORY, f'{NAME}_{STATE}_isolated.inp'), main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)

    if CREATE_TDDFT_CALC:
        create_tddft_input(os.path.join(OUTPUT_DIRECTORY, f'{NAME}_{STATE}_tddft.inp'), main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)
