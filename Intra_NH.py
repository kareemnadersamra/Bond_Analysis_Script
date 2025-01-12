#!/usr/bin/env python3

import os
import sys
import numpy as np
from ase.io import read

# Helper functions
def minimum_image_distance(atom1, atom2, cell):
    """
    Calculate the shortest distance between two atoms using the minimum image convention.
    """
    delta = atom1 - atom2
    fractional_delta = np.dot(delta, np.linalg.inv(cell))
    fractional_delta -= np.round(fractional_delta)  # Apply periodic boundary conditions
    cartesian_delta = np.dot(fractional_delta, cell)
    return np.linalg.norm(cartesian_delta)

def calculate_angle(atom1, atom2, atom3, cell):
    """
    Calculate the angle (in degrees) formed by three atoms, considering periodicity.
    """
    vec1 = atom1 - atom2
    vec2 = atom3 - atom2

    # Apply minimum image convention to vectors
    fractional_vec1 = np.dot(vec1, np.linalg.inv(cell))
    fractional_vec1 -= np.round(fractional_vec1)
    cartesian_vec1 = np.dot(fractional_vec1, cell)

    fractional_vec2 = np.dot(vec2, np.linalg.inv(cell))
    fractional_vec2 -= np.round(fractional_vec2)
    cartesian_vec2 = np.dot(fractional_vec2, cell)

    cos_theta = np.dot(cartesian_vec1, cartesian_vec2) / (
        np.linalg.norm(cartesian_vec1) * np.linalg.norm(cartesian_vec2)
    )
    cos_theta = max(-1.0, min(1.0, cos_theta))  # Clamp to avoid numerical issues
    return np.degrees(np.arccos(cos_theta))

def find_inter_hbonds(atom_positions, atom_types, cell, max_distance=3.5, min_distance=1.5, min_angle=120):
    """
    Identify realistic intermolecular hydrogen bonds based on length and angle criteria.
    """
    hydrogen_bonds = []
    num_atoms = len(atom_types)

    for i, donor_type in enumerate(atom_types):
        if donor_type not in ['O', 'N']:  # Donor must be O or N
            continue

        # Find hydrogens covalently bonded to the donor
        for j, hydrogen_type in enumerate(atom_types):
            if hydrogen_type != 'H':  # Skip non-hydrogen atoms
                continue

            # Verify donor-hydrogen bond length (covalent bond)
            d_h_distance = minimum_image_distance(atom_positions[i], atom_positions[j], cell)
            if d_h_distance > 1.1:  # Covalent bond threshold
                continue

            # Check for acceptor atoms
            for k, acceptor_type in enumerate(atom_types):
                if acceptor_type != 'O' or k == i or k == j:  # Acceptor must be oxygen, different from donor/hydrogen
                    continue

                # Calculate H---O distance
                h_o_distance = minimum_image_distance(atom_positions[j], atom_positions[k], cell)
                if not (min_distance <= h_o_distance <= max_distance):  # Realistic hydrogen bond length
                    continue

                # Calculate D-H---O angle
                angle = calculate_angle(atom_positions[i], atom_positions[j], atom_positions[k], cell)
                if angle < min_angle:  # Realistic hydrogen bond angle
                    continue

                # Add valid hydrogen bond to the list
                hydrogen_bonds.append({
                    'Donor': i,
                    'Hydrogen': j,
                    'Acceptor': k,
                    'D-H Distance (Å)': d_h_distance,
                    'H---O Distance (Å)': h_o_distance,
                    'D-H---O Angle (°)': angle
                })

    return hydrogen_bonds

def process_xyz_files(directory):
    """
    Process all .xyz files in the directory and allow user to select specific hydrogen bonds to track.
    """
    xyz_files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.xyz')]
    if not xyz_files:
        print(f"No .xyz files found in the directory: {directory}")
        sys.exit(1)

    # Process the first .xyz file to select indices
    first_file = xyz_files[0]
    base_name = first_file[:-4]
    cell_file = os.path.join(directory, f"{base_name}.txt")

    if not os.path.exists(cell_file):
        print(f"Missing cell matrix file for {first_file}. Skipping.")
        sys.exit(1)

    with open(cell_file, 'r') as f:
        cell = np.loadtxt(f).reshape(3, 3)

    structure = read(first_file)
    atom_positions = structure.get_positions()
    atom_types = structure.get_chemical_symbols()

    # Find hydrogen bonds in the first file
    hydrogen_bonds = find_inter_hbonds(atom_positions, atom_types, cell)

    print(f"\n--- Processing First File: {first_file} ---")
    print(f"  Detected {len(hydrogen_bonds)} hydrogen bonds:")
    for idx, hbond in enumerate(hydrogen_bonds):
        print(f"    {idx}: Donor {hbond['Donor']} - Hydrogen {hbond['Hydrogen']} - Acceptor {hbond['Acceptor']} "
              f"H---O Distance: {hbond['H---O Distance (Å)']:.3f}, D-H---O Angle: {hbond['D-H---O Angle (°)']:.2f}")

    # User selects indices to track
    selected_indices = input("Enter the indices of hydrogen bonds to track (comma-separated): ")
    selected_indices = [int(idx.strip()) for idx in selected_indices.split(',') if idx.strip().isdigit()]

    # Process all .xyz files
    for filename in sorted(xyz_files):
        base_name = filename[:-4]
        cell_file = os.path.join(directory, f"{base_name}.txt")

        if not os.path.exists(cell_file):
            print(f"Missing cell matrix file for {filename}. Skipping.")
            continue

        with open(cell_file, 'r') as f:
            cell = np.loadtxt(f).reshape(3, 3)

        structure = read(filename)
        atom_positions = structure.get_positions()
        atom_types = structure.get_chemical_symbols()

        # Find hydrogen bonds
        hydrogen_bonds = find_inter_hbonds(atom_positions, atom_types, cell)

        print(f"\nProcessing {filename}:")
        print(f"  Detected {len(hydrogen_bonds)} hydrogen bonds.")

        # Write only the selected bonds to files
        angle_file = "tracked_angle.txt"
        distance_file = "tracked_distance.txt"

        with open(angle_file, 'a') as af, open(distance_file, 'a') as df:
            for idx in selected_indices:
                if idx < len(hydrogen_bonds):
                    hbond = hydrogen_bonds[idx]
                    af.write(f"{base_name}: {hbond['D-H---O Angle (°)']:.2f}\n")
                    df.write(f"{base_name}: {hbond['H---O Distance (Å)']:.3f}\n")

if __name__ == "__main__":
    directory = input("Enter the directory path containing .xyz and .txt files: ").strip()

    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        sys.exit(1)

    process_xyz_files(directory)
