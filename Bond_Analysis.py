#!/usr/bin/env python3

import os
import sys
import json
import numpy as np
from ase.io import read

# Bond cutoff distances for specific bond types

bond_cutoffs = {
    ('C', 'C'): 1.54,
    ('C', 'H'): 1.09,
    ('C', 'O'): 1.50,  # Increased slightly
    ('C', 'N'): 1.50,  # Increased slightly
    ('C', 'O(double)'): 1.30,  # Increased slightly
    ('O', 'H'): 1.00,  # Increased slightly
    ('N', 'H'): 1.10,  # Increased slightly
}

def get_xyz_files(directory="."):
    """Fetch all .xyz files from a directory."""
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        sys.exit(1)
    
    xyz_files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.xyz')]
    if not xyz_files:
        print(f"No .xyz files found in the directory: {directory}")
        sys.exit(1)
    
    return xyz_files
# Helper functions
def minimum_image_distance(atom1, atom2, cell):
    delta = atom1 - atom2
    fractional_delta = np.dot(delta, np.linalg.inv(cell))  # Convert to fractional coordinates
    fractional_delta -= np.round(fractional_delta)  # Apply periodic boundary conditions
    cartesian_delta = np.dot(fractional_delta, cell)  # Convert back to Cartesian coordinates
    return np.linalg.norm(cartesian_delta)

def calculate_angle(atom1, atom2, atom3, cell):
    vec1 = atom1 - atom2
    vec2 = atom3 - atom2
    vec1 = vec1 - np.dot(np.round(np.dot(vec1, np.linalg.inv(cell))), cell)
    vec2 = vec2 - np.dot(np.round(np.dot(vec2, np.linalg.inv(cell))), cell)
    cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    cos_theta = max(-1.0, min(1.0, cos_theta))  # Clamp cos_theta to prevent numerical issues
    return np.degrees(np.arccos(cos_theta))


#################
def find_intra_bonds(atom_positions, atom_types, cell, bond_cutoffs):
    """
    Identify intramolecular NH, CO(double), and CO bonds for Paracetamol.
    """
    intra_bonds = {"NH": [], "CO(double)": [], "CO": []}
    num_atoms = len(atom_types)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            bond_type = (atom_types[i], atom_types[j])
            reverse_bond_type = bond_type[::-1]

            if bond_type in bond_cutoffs or reverse_bond_type in bond_cutoffs:
                distance = minimum_image_distance(atom_positions[i], atom_positions[j], cell)
                cutoff = bond_cutoffs.get(bond_type, bond_cutoffs.get(reverse_bond_type))

                if distance <= cutoff:
                    # Debugging: Print bond candidates
                    print(f"Evaluating bond: {bond_type} (distance: {distance:.3f} Å, cutoff: {cutoff:.3f} Å)")

                    # NH bond
                    if bond_type == ('N', 'H') or reverse_bond_type == ('N', 'H'):
                        print(f"  -> Possible NH bond detected between {i} (N) and {j} (H)")
                        # Check if nitrogen is bonded to a carbon
                        nitrogen_idx = i if atom_types[i] == 'N' else j
                        hydrogen_idx = j if atom_types[i] == 'N' else i

                        bonded_to_carbon = any(
                            minimum_image_distance(atom_positions[nitrogen_idx], atom_positions[k], cell)
                            <= bond_cutoffs.get(('N', 'C'), bond_cutoffs.get(('C', 'N')))
                            for k, atom_type in enumerate(atom_types)
                            if atom_type == 'C' and k != nitrogen_idx and k != hydrogen_idx
                        )

                        if bonded_to_carbon:
                            intra_bonds["NH"].append({
                                'Nitrogen': nitrogen_idx,
                                'Hydrogen': hydrogen_idx,
                                'Bond_Length (Å)': distance
                            })
                            print(f"    -> NH bond confirmed between {nitrogen_idx} (N) and {hydrogen_idx} (H)")

                    # CO(double) bond
                    if bond_type == ('C', 'O') or reverse_bond_type == ('C', 'O'):
                        print(f"  -> Possible CO(double) bond detected between {i} (C) and {j} (O)")
                        carbon_idx = i if atom_types[i] == 'C' else j
                        oxygen_idx = j if atom_types[i] == 'C' else i

                        # Check if carbon is bonded to NH and another carbon with three hydrogens
                        bonded_to_nh = any(
                            minimum_image_distance(atom_positions[carbon_idx], atom_positions[k], cell)
                            <= bond_cutoffs.get(('C', 'N'))
                            and any(
                                minimum_image_distance(atom_positions[k], atom_positions[h], cell)
                                <= bond_cutoffs.get(('N', 'H'))
                                for h, atom_type_h in enumerate(atom_types) if atom_type_h == 'H'
                            )
                            for k, atom_type in enumerate(atom_types) if atom_type == 'N'
                        )

                        bonded_to_carbon_with_hydrogens = any(
                            minimum_image_distance(atom_positions[carbon_idx], atom_positions[k], cell)
                            <= bond_cutoffs.get(('C', 'C'))
                            and sum(
                                1 for h, atom_type_h in enumerate(atom_types)
                                if atom_type_h == 'H' and
                                minimum_image_distance(atom_positions[k], atom_positions[h], cell)
                                <= bond_cutoffs.get(('C', 'H'))
                            ) == 3
                            for k, atom_type in enumerate(atom_types) if atom_type == 'C'
                        )

                        if bonded_to_nh and bonded_to_carbon_with_hydrogens:
                            intra_bonds["CO(double)"].append({
                                'Carbon': carbon_idx,
                                'Oxygen': oxygen_idx,
                                'Bond_Length (Å)': distance
                            })
                            print(f"    -> CO(double) bond confirmed between {carbon_idx} (C) and {oxygen_idx} (O)")

                    # CO bond
                    if bond_type == ('C', 'O') or reverse_bond_type == ('C', 'O'):
                        intra_bonds["CO"].append({
                            'Carbon': i if atom_types[i] == 'C' else j,
                            'Oxygen': j if atom_types[i] == 'C' else i,
                            'Bond_Length (Å)': distance
                        })
                        print(f"    -> CO bond confirmed between {i} (C) and {j} (O)")

    return intra_bonds


#################


def calculate_bond_lengths(atom_positions, atom_types, cell, bond_cutoffs):
    """
    Calculate bond lengths and classify them into categories.
    """
    bond_lengths = {"NH": [], "COH": [], "CO(double)": [], "C-C": []}
    num_atoms = len(atom_types)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):  # Only calculate for pairs (i, j)
            bond_type = (atom_types[i], atom_types[j])
            reverse_bond_type = bond_type[::-1]

            if bond_type in bond_cutoffs or reverse_bond_type in bond_cutoffs:
                distance = minimum_image_distance(atom_positions[i], atom_positions[j], cell)
                cutoff = bond_cutoffs.get(bond_type, bond_cutoffs.get(reverse_bond_type))

                if distance <= cutoff:
                    # Classify bonds
                    if bond_type == ('N', 'H') or reverse_bond_type == ('N', 'H'):
                        bond_lengths["NH"].append({
                            'Atom1': i,
                            'Atom2': j,
                            'Bond_Type': bond_type,
                            'Bond_Length (Å)': distance
                        })
                    elif bond_type == ('C', 'O(double)') or reverse_bond_type == ('C', 'O(double)'):
                        bond_lengths["CO(double)"].append({
                            'Atom1': i,
                            'Atom2': j,
                            'Bond_Type': bond_type,
                            'Bond_Length (Å)': distance
                        })
                    elif bond_type == ('C', 'O') or reverse_bond_type == ('C', 'O'):
                        bond_lengths["COH"].append({
                            'Atom1': i,
                            'Atom2': j,
                            'Bond_Type': bond_type,
                            'Bond_Length (Å)': distance
                        })
                    elif bond_type == ('C', 'C') or reverse_bond_type == ('C', 'C'):
                        bond_lengths["C-C"].append({
                            'Atom1': i,
                            'Atom2': j,
                            'Bond_Type': bond_type,
                            'Bond_Length (Å)': distance
                        })

    return bond_lengths


def find_inter_hbonds(atom_positions, atom_types, donor_pair, acceptor, cell, max_distance=2.5, min_angle=150):
    hydrogen_bonds = []
    num_atoms = len(atom_types)
    for i, atom_type in enumerate(atom_types):
        if atom_type == donor_pair[0]:
            for j, atom_type_h in enumerate(atom_types):
                if atom_type_h == donor_pair[1]:
                    for k, acceptor_type in enumerate(atom_types):
                        if acceptor_type == acceptor:
                            h_o_distance = minimum_image_distance(atom_positions[j], atom_positions[k], cell)
                            if h_o_distance <= max_distance:
                                angle = calculate_angle(atom_positions[i], atom_positions[j], atom_positions[k], cell)
                                if angle >= min_angle:
                                    hydrogen_bonds.append({
                                        'Donor': i,
                                        'Hydrogen': j,
                                        'Acceptor': k,
                                        'H---O Distance (Å)': h_o_distance,
                                        'Donor-H---O Angle (°)': angle
                                    })
    return hydrogen_bonds

def create_master_selections(reference_xyz, bond_cutoffs, output_file):
    """
    Create master selections using a reference `.xyz` file.
    """
    cell_file = f"{reference_xyz[:-4]}.txt"

    if not os.path.exists(cell_file):
        print(f"Error: Cell file {cell_file} not found for {reference_xyz}.")
        exit()

    # Load the cell matrix
    try:
        with open(cell_file, 'r') as f:
            cell = np.loadtxt(f).reshape(3, 3)
    except Exception as e:
        print(f"Error reading the cell matrix for {reference_xyz}: {e}")
        exit()

    # Load the molecular structure
    structure = read(reference_xyz)
    atom_positions = structure.get_positions()
    atom_types = structure.get_chemical_symbols()

    # Extract bond lengths and hydrogen bonds
    bond_lengths = calculate_bond_lengths(atom_positions, atom_types, cell, bond_cutoffs)
    nh_hbonds = find_inter_hbonds(atom_positions, atom_types, donor_pair=('N', 'H'), acceptor='O', cell=cell)
    oh_hbonds = find_inter_hbonds(atom_positions, atom_types, donor_pair=('O', 'H'), acceptor='O', cell=cell)

    # Create master selections interactively
    master_selections = {"bonds": {}, "angles": {"NH--O": [], "OH--O": []}}

    print("\n--- Intramolecular Bonds ---")
    for category, bonds in bond_lengths.items():
        print(f"\n{category} Bonds:")
        for idx, bond in enumerate(bonds):
            print(f"  {idx}: Bond {bond['Bond_Type']} between Atom {bond['Atom1']} and Atom {bond['Atom2']} "
                  f"Length: {bond['Bond_Length (Å)']:.6f}")
        selected_indices = input(f"Enter the indices of {category} bonds to track (comma-separated): ")
        selected_indices = [int(idx.strip()) for idx in selected_indices.split(',') if idx.strip().isdigit()]
        master_selections["bonds"][category] = selected_indices

    print("\n--- Intermolecular Hydrogen Bonds: NH--O ---")
    for idx, hbond in enumerate(nh_hbonds):
        print(f"  {idx}: Donor {hbond['Donor']} - Hydrogen {hbond['Hydrogen']} - Acceptor {hbond['Acceptor']} "
              f"Angle: {hbond['Donor-H---O Angle (°)']:.2f}, Distance: {hbond['H---O Distance (Å)']:.6f}")
    nh_indices = input("Enter the indices of NH--O bonds to track (comma-separated): ")
    master_selections["angles"]["NH--O"] = [int(idx.strip()) for idx in nh_indices.split(',') if idx.strip().isdigit()]

    print("\n--- Intermolecular Hydrogen Bonds: OH--O ---")
    for idx, hbond in enumerate(oh_hbonds):
        print(f"  {idx}: Donor {hbond['Donor']} - Hydrogen {hbond['Hydrogen']} - Acceptor {hbond['Acceptor']} "
              f"Angle: {hbond['Donor-H---O Angle (°)']:.2f}, Distance: {hbond['H---O Distance (Å)']:.6f}")
    oh_indices = input("Enter the indices of OH--O bonds to track (comma-separated): ")
    master_selections["angles"]["OH--O"] = [int(idx.strip()) for idx in oh_indices.split(',') if idx.strip().isdigit()]

    # Save to output file
    with open(output_file, 'w') as f:
        json.dump(master_selections, f)
    print(f"Master selections saved to {output_file}.")

def apply_fixed_indices(master_indices, bond_lengths, nh_hbonds, oh_hbonds):
    """
    Apply fixed indices from a reference selection to the current file.
    """
    selections = {"bonds": {}, "angles": {"NH--O": [], "OH--O": []}}

    # Apply fixed indices for bonds
    for category, indices in master_indices["bonds"].items():
        selections["bonds"][category] = [
            bond_lengths[category][idx] for idx in indices if idx < len(bond_lengths[category])
        ]

    # Apply fixed indices for angles
    for category, indices in master_indices["angles"].items():
        current_hbonds = nh_hbonds if category == "NH--O" else oh_hbonds
        selections["angles"][category] = [
            current_hbonds[idx] for idx in indices if idx < len(current_hbonds)
        ]

    return selections

# Main processing script
try:
    
    # Directory containing .xyz and .txt files
    directory = ""     #Please Enter Directory

    # Get all .xyz and .txt files
    xyz_files = get_xyz_files(directory)
    txt_files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.txt')]

    # Ensure there's a corresponding .txt file for each .xyz file
    if len(txt_files) < len(xyz_files):
        print("Error: Not enough .txt files for the .xyz files.")
        sys.exit(1)

    # Load or create master selections
    master_selection_file = "master_selections.json"
    if not os.path.exists(master_selection_file):
        reference_xyz = input("Enter the reference `.xyz` file name: ")
        create_master_selections(reference_xyz, bond_cutoffs, master_selection_file)

    with open(master_selection_file, 'r') as f:
        master_selections = json.load(f)
    print(f"Loaded master selections from {master_selection_file}.")

    # Process each XYZ file with its corresponding cell matrix file
    for filename in sorted(os.listdir(directory)):
        if filename.endswith('.xyz'):
            base_name = filename[:-4]  # Remove ".xyz" to find matching .txt file
            cell_file = os.path.join(directory, f"{base_name}.txt")

            if not os.path.exists(cell_file):
                print(f"Missing cell matrix file for {filename}. Skipping.")
                continue

            # Load the cell matrix
            try:
                with open(cell_file, 'r') as f:
                    cell = np.loadtxt(f).reshape(3, 3)
            except Exception as e:
                print(f"Error reading the cell matrix for {filename}: {e}")
                continue

            # Load the molecular structure
            structure = read(os.path.join(directory, filename))
            atom_positions = structure.get_positions()
            atom_types = structure.get_chemical_symbols()

            # Extract bond lengths and hydrogen bonds
            bond_lengths = calculate_bond_lengths(atom_positions, atom_types, cell, bond_cutoffs)
            nh_hbonds = find_inter_hbonds(atom_positions, atom_types, donor_pair=('N', 'H'), acceptor='O', cell=cell)
            oh_hbonds = find_inter_hbonds(atom_positions, atom_types, donor_pair=('O', 'H'), acceptor='O', cell=cell)
            intra_bonds = find_intra_bonds(atom_positions, atom_types, cell, bond_cutoffs)

            print(f"\nIntra Bonds for {filename}:")
            for category, bonds in intra_bonds.items():
                print(f"{category}: {len(bonds)} bonds detected.")
                for bond in bonds:
                    print(f"  {bond}")

            # Apply fixed indices from the master selections
            selections = apply_fixed_indices(master_selections, bond_lengths, nh_hbonds, oh_hbonds)

            # Append tracked bond lengths
            for category, bonds in selections["bonds"].items():
                bond_file = f"{category}_tracked_bond_length.txt"
                with open(bond_file, 'a') as f:
                    for bond in bonds:
                        f.write(f"{base_name}: {bond['Bond_Length (Å)']:.6f}\n")

            # Append tracked angles
            for category, angles in selections["angles"].items():
                angle_file = f"{category}_tracked_angle.txt"
                distance_file = f"{category}_tracked_distance.txt"
                with open(angle_file, 'a') as f:
                    for angle in angles:
                        f.write(f"{base_name}: {angle['Donor-H---O Angle (°)']:.2f}\n")
                with open(distance_file, 'a') as f:
                    for angle in angles:
                        f.write(f"{base_name}: {angle['H---O Distance (Å)']:.6f}\n")

            print(f"Processed {filename} with its corresponding cell matrix {cell_file}.")

except Exception as e:
    print(f"Error processing files: {e}")