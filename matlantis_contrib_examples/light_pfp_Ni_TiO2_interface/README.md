# LightPFP model training and valiation for Ni / anatase TiO2 interface

## How to use the notebooks
* Please run `model_training_Ni_TiO2.ipynb` and `pfp_md.ipynb` at first, the two notebooks can be run in parallel.
  * `model_training_Ni_TiO2.ipynb`: Dataset collection and model training with active learning method.
  * `pfp_md.ipynb`: Generate several MD trajectories with PFP for the comparision and valiation of LightPFP performance.
* After both notebooks are finished, please run `validate.ipynb`.
  * `validate.ipynb`: Compare "density", "radial distribution function" and "diffusion behavior" between the LightPFP trajectories and PFP trajectories.
* At last, please run "interface_energy.ipynb"
  * `interface_energy.ipynb`: Calculate the interface energies with PFP and LightPFP

## Appendix notebook
* `appendix_make_interface.ipynb`: How to generate interface structures with pymatgen? These structures are severed as initial structures of validation MD tasks.
* `appendix_view_training_structure.ipynb`: How to visualize the atomic structure from the H5DF training dataset?
* `appendix_view_trajectory.ipynb`: How to visualize the MD trajectories?

## Results summary
* The results are summerized in `outputs/results.pdf`