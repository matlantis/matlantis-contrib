### Notion
- These scripts consume large storage (> 10 GB) and memory (> 16GB on learning).

### How to use
Please execute the script in order.
The summary of each script are shown as below:
01_extract_from_qm9_dataset: the structure and magnetic shielding value was extracted from raw data
02_get_descriptors: the PFP descriptors was obtained from the structure
03_prepare_dataset: the dataset for prediction model was generated
04_pycaret_screening: brief model screening with pycaret
05_condition_screening: preprocess and other condition screening with selected model
06_analyze_condition: the analysis of condition screening
07_hyperparameter_tuning: the tuning with selected model and conditions to obtain best prediction model

### Results
- Finally we can achieve the SOTA result for 1H NMR chemical shift prediction for QM9NMR dataset.
  - valid MAE: 0.1156 ppm

### Dependency
The script was confirmed with the libraries below:
- python 3.11.16  
- lightgbm 4.6.0
- numpy 1.26.4
- optuna4.4.0
- pandas 2.3.2
- scikit-learn 1.7.1
- xgboost 3.0.2
- ase 3.25.0
- pfp-api-client 1.23.1
- pfcc-extras 0.12.1
- rdkit 2024.9.4

### Citation
- [QM9NMR dataset](https://moldis-group.github.io/qm9nmr/)
