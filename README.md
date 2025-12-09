TRIM29-UC Mechanism Analysis
This repository contains the reproducible R code for the manuscript "TRIM29 Drives Ulcerative Colitis by Disrupting Lipid Metabolism via Lysosomal Dysfunction: A Multi-Omics and Experimental Study" which is under review.

**Quick Reproducibility Guide**

1. Clone & Setup
git clone https://github.com/yourusername/TRIM29_UC_Mechanism.git
cd TRIM29_UC_Mechanism

2. Restore R Environment

**Install renv if needed**

install.packages("renv")

renv::restore()  # Installs exact package versions

3. Run Full Analysis
source("run_analysis.R")

4. Find Results

**All outputs are saved to results/:**

Figures: results/figures/

Tables: results/tables/

Data objects: results/*.RData

**Code Structure**

TRIM29_UC_Mechanism/

├── run_analysis.R              # Master scrips

├── renv.lock                   # Package versions

├── code/                       # Core analysis

│   ├── 01_WGCNA_analysis.R

│   ├── 02_mendelian_randomization.R

│   ├── 03_machine_learning.R

│   └── 04_xgboost_selection.R

**Citation**
If using this code, please cite our manuscript after. 
For questions, contact [gxz479233@163.com].

**License**
MIT License. See LICENSE file.

**Markdown**
- This code repository is provided **solely for the purpose of peer review and academic verification** of the associated manuscript.
- The authors make no warranties regarding the accuracy or completeness of the code.
- **This code is not intended for clinical use or direct medical decision-making.**
- Users are responsible for complying with all applicable data use agreements for any input data they employ.
