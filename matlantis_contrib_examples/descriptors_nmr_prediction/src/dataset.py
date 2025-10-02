###
#This module was provided from Preferred Networks Inc.
#Copyright Preferred Networks Inc. as contributors to Matlantis contrib project
###

import h5py
import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np

class HDF5Dataset(Dataset):
    def __init__(self,
        hdf5_path,
        data_scaler = None,
        keys = None,
    ):
        """
        Args:
            hdf5_path (str): Path to the HDF5 file.
        """
        self.hdf5_path = hdf5_path
        self.data_scaler = data_scaler
        
        if keys:
            self.molecule_ids = keys
        else:
            with h5py.File(self.hdf5_path, "r") as f:
                self.molecule_ids = list(f.keys())

    def __len__(self):
        return len(self.molecule_ids)

    def __getitem__(self, idx):
        with h5py.File(self.hdf5_path, "r") as f:
            mol_id = self.molecule_ids[idx]
            mol_group = f[mol_id]

            molecule_data = []
            for atom_id in mol_group.keys():
                atom_group = mol_group[atom_id]

                # per-atom features
                original_value = np.reshape(atom_group["value"][()], (1, 1))   # (1, 1)                
                value = torch.tensor(original_value, dtype=torch.float32).squeeze()

                h_s = torch.tensor(atom_group["scalar"][()], dtype=torch.float32)

                molecule_data.append({
                    "value": value,
                    "scalar": h_s,
                })

        value_all = torch.stack([mol['value'] for mol in molecule_data])
        h_s_all = torch.stack([mol['scalar'] for mol in molecule_data])

        return value_all, h_s_all
