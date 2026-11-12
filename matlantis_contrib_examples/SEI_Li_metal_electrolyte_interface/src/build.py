"""Structure building helpers for the Li-metal / carbonate-electrolyte interface example.

Molecules from SMILES (`smiles_to_atoms`), connectivity checks (`get_connected_components`) and the Packmol packing
around the slab (`LiquidGenerator`) are provided by pfcc_extras; this module contains only what pfcc_extras does not.
"""

import numpy as np
from ase import Atoms


def build_pf6(d_pf=1.62):
    """Return an octahedral PF6 with P-F distance d_pf (A); RDKit cannot embed hexacoordinated P."""
    directions = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    return Atoms("PF6", positions=np.vstack([[0.0, 0.0, 0.0], d_pf * directions]))


def assign_layers(z, tol=0.5):
    """Return the layer index (0 = bottom) of each atom, grouping z coordinates closer than tol (A)."""
    order = np.argsort(z)
    layer = np.empty(len(z), dtype=int)
    current = 0
    layer[order[0]] = 0
    for prev, idx in zip(order[:-1], order[1:]):
        if z[idx] - z[prev] > tol:
            current += 1
        layer[idx] = current
    return layer


def molecule_ranges(composition, molecules, start=0):
    """Atom index range [start, end) of every molecule, in the order Packmol writes them (composition order).

    composition: list of {"name": ..., "number": ...}; molecules: {name: Atoms}.
    """
    table = []
    for entry in composition:
        n_atoms = len(molecules[entry["name"]])
        for k in range(entry["number"]):
            table.append({"label": f"{entry['name']}_{k:03d}", "species": entry["name"],
                          "start": start, "end": start + n_atoms})
            start += n_atoms
    return table
