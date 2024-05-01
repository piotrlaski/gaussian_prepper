import numpy as np
import itertools
import os
import re
import time
import functools

def timing_decorator(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Function '{func.__name__}' took {elapsed_time:.4f} seconds to complete.")
        return result
    return wrapper


def apply_symmetry(old_fracs, sym_op):
    x, y, z = old_fracs
    new_x_frac, new_y_frac, new_z_frac = eval(sym_op[0]), eval(sym_op[1]), eval(sym_op[2])
    return ([new_x_frac, new_y_frac, new_z_frac])

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
    @timing_decorator
    def build_infinite_crystal(self, main_molecule, p = 3):
        for sym_op in self.sym_ops:
            for mod_x, mod_y, mod_z in itertools.product(range(-p, p+1), repeat=3):
                new_sym_id = [f'{sym_op[0]}+{str(mod_x)}', f'{sym_op[1]}+{str(mod_y)}', f'{sym_op[2]}+{str(mod_z)}']
                self.spawn_sym_mate(main_molecule, sym_op = new_sym_id)
    @timing_decorator
    def cut_out_cluster(self, main_molecule, radius):
        deletion_list = []
        for molecule in self.get_molecules():
            for atom in molecule.get_atoms():
                if any([dist([main_atom.x_ort, main_atom.y_ort, main_atom.z_ort], [atom.x_ort, atom.y_ort, atom.z_ort]) < radius for main_atom in main_molecule.get_atoms()]):
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

def spawn_crystal(cell_params, sym_ops, loaded_data):
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
    return (xtal, main_molecule)

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