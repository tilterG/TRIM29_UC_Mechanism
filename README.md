TRIM29-UC Mechanism Analysis
This repository contains the reproducible R code for the manuscript "TRIM29 Drives Ulcerative Colitis by Disrupting Lipid Metabolism via Lysosomal Dysfunction: A Multi-Omics and Experimental Study".

Quick Reproducibility Guide
1. Clone & Setup
git clone https://github.com/yourusername/TRIM29_UC_Mechanism.git
cd TRIM29_UC_Mechanism
2. Restore R Environment
# Install renv if needed
install.packages("renv")
renv::restore()  # Installs exact package versions
3. Run Full Analysis
source("run_analysis.R")
4. Find Results
All outputs are saved to results/:

Figures: results/figures/

Tables: results/tables/

Data objects: results/*.RData

Code Structure

TRIM29_UC_Mechanism/

├── run_analysis.R              # Master script

├── renv.lock                   # Package versions

├── code/                       # Core analysis

│   ├── 01_WGCNA_analysis.R

│   ├── 02_mendelian_randomization.R

│   ├── 03_machine_learning.R

│   └── 04_xgboost_selection.R

Citation
If using this code, please cite our manuscript. For questions, contact [your.name@institution.edu].

License
MIT License. See LICENSE file.
