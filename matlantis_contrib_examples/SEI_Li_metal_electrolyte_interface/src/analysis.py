"""Trajectory analysis utilities for the Li-metal / carbonate-electrolyte interface example."""

from collections import Counter

import numpy as np
from ase import Atoms
from ase.formula import Formula
from pfcc_extras.structure.connectivity import get_connected_components, get_connectivity_matrix
from scipy.sparse import coo_matrix


def distance_from_surface(positions, cell, fixed_idx, top_idx, bottom_idx):
    """Signed distance of every atom from the nearer Li surface (frames x atoms).

    The fixed central layers never move, so their mean height is a stable reference in every frame.
    The surfaces are the median heights of the outermost Li layers. The distance is positive in the
    electrolyte and negative inside the slab.
    """
    lz = cell[2, 2]
    z_center = positions[:, fixed_idx, 2].mean(axis=1)
    z_rel = (positions[:, :, 2] - z_center[:, None] + lz / 2) % lz - lz / 2
    top_plane = np.median(z_rel[:, top_idx], axis=1)
    bottom_plane = np.median(z_rel[:, bottom_idx], axis=1)
    dist = np.where(z_rel >= 0, z_rel - top_plane[:, None], bottom_plane[:, None] - z_rel)
    return dist, top_plane, bottom_plane


def hill(symbols):
    """Chemical formula in Hill notation, e.g. ['O', 'C'] -> 'CO'."""
    return Formula.from_list(list(symbols)).format("hill")


def bond_analysis(frames, atom_indices, mol_of_atom, n_molecules, bond_scale=1.2):
    """Covalent fragments in every frame, using the connectivity of pfcc_extras.

    atom_indices: atoms included in the bond network (Li is usually excluded).
    mol_of_atom: index of the original molecule of every atom of the system.
    bond_scale: bonded if distance < bond_scale x (sum of covalent radii), the `mult` of pfcc_extras.

    Returns
        n_fragments: (frames x molecules) number of fragments of each original molecule
        intra_labels: per frame, fragment label of each atom in atom_indices (bonds within one molecule only)
        census: per frame, Counter of the formulas of all bonded species (bonds between any atoms)
    """
    atom_indices = np.asarray(atom_indices)
    symbols = np.array(frames[0].get_chemical_symbols())[atom_indices]
    mol = mol_of_atom[atom_indices]

    n_fragments = np.zeros((len(frames), n_molecules), dtype=int)
    intra_labels, census = [], []
    for f, atoms in enumerate(frames):
        sub = atoms[atom_indices]
        bonds = get_connectivity_matrix(sub, mult=bond_scale).tocoo()
        _, all_labels = get_connected_components(sub, matrix=bonds)

        same = mol[bonds.row] == mol[bonds.col]  # keep only bonds within one original molecule
        intra = coo_matrix((bonds.data[same], (bonds.row[same], bonds.col[same])), shape=bonds.shape)
        _, labels = get_connected_components(sub, matrix=intra)
        pairs = np.unique(np.stack([mol, labels], axis=1), axis=0)
        n_fragments[f] = np.bincount(pairs[:, 0], minlength=n_molecules)
        intra_labels.append(labels)

        groups = {}
        for label, symbol in zip(all_labels, symbols):
            groups.setdefault(label, []).append(symbol)
        census.append(Counter(hill(g) for g in groups.values()))
    return n_fragments, intra_labels, census


def fragment_formula(symbols, labels):
    """Formulas of the fragments of one molecule, e.g. 'C2H4O2 + CO' or '6 F + P'."""
    groups = {}
    for label, symbol in zip(labels, symbols):
        groups.setdefault(label, []).append(symbol)
    counted = Counter(hill(g) for g in groups.values())
    ordered = sorted(counted.items(), key=lambda item: (-len(item[0]), item[0]))
    return " + ".join(name if n == 1 else f"{n} {name}" for name, n in ordered)


def first_persistent(mask, width):
    """Index of the first frame from which mask stays True for `width` frames (None if never)."""
    if width <= 0:
        raise ValueError("width must be positive")
    for f in range(len(mask) - width + 1):
        if mask[f:f + width].all():
            return f
    return None


def whole_molecule(atoms, start, end):
    """Atoms start:end made whole across the periodic boundary (for display)."""
    sub = atoms[start:end]
    vectors = sub.get_distances(0, range(len(sub)), mic=True, vector=True)
    return Atoms(sub.get_chemical_symbols(), positions=sub.positions[0] + vectors)


def region_counts(dist, element_masks, cutoff):
    """Atoms of each element within `cutoff` of the surface, per frame (frames x elements)."""
    inside = dist <= cutoff
    return np.array([[np.count_nonzero(row & em) for em in element_masks] for row in inside])


def to_percent(counts, masses, basis):
    """Convert element counts to atomic % (basis='atomic') or weight % (basis='weight')."""
    weights = counts * masses if basis == "weight" else counts
    return 100 * weights / weights.sum(axis=-1, keepdims=True)
