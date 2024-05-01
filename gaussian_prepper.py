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

    genecp = []
    for genecp_file in GENECP_FILE:
        with open(genecp_file, '+r') as f:
            genecp.append(f.read())


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
        save_xyz(os.path.join(OUTPUT_DIRECTORY, f'counterpoise_dimer.xyz'), xtal)
        multiple_inputs(input_creation_function=create_counterpoise_input, outputdir=OUTPUT_DIRECTORY, name=NAME, crystal = xtal, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)

    if CREATE_ONIOM:
        xtal.build_infinite_crystal(main_molecule)
        xtal.cut_out_cluster(main_molecule=main_molecule, radius = RADIUS)
        save_xyz(os.path.join(OUTPUT_DIRECTORY, f'cluster.xyz'), xtal)
        multiple_inputs(input_creation_function = create_QMMM_input, outputdir=OUTPUT_DIRECTORY, name=NAME, crystal=xtal, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)
    if CREATE_CHARGE_CALC:
        multiple_inputs(input_creation_function=create_charge_calc_input, outputdir=OUTPUT_DIRECTORY, name=NAME, molecule = main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, additional = ADDITIONAL, genecp = genecp)

    if CREATE_ISOLATED_OPT:
        multiple_inputs(input_creation_function=create_iso_opt_input, outputdir=OUTPUT_DIRECTORY, name=NAME, molecule=main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)

    if CREATE_TDDFT_CALC:
        multiple_inputs(input_creation_function=create_tddft_input, outputdir=OUTPUT_DIRECTORY, name=NAME, molecule=main_molecule, functional = FUNCTIONAL, base = BASE, state = STATE, nstates = NSTATES, additional = ADDITIONAL, genecp = genecp)
