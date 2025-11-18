# CoDaTools - Source Code

This directory contains MATLAB functions for **Compositional Data Analysis (CoDa)**. The tools are organized into three main categories:

---

## 📂 src Folder Structure

```
src/
│
├── check/        # Functions for data validation and consistency checks
│   ├── check_closure.m   # Checks if data satisfies closure property.
│   ├── check_clr.m       # Verifies if data satisfies the centered log-ratio (clr) property, i.e., the sum of clr-transformed components equals zero.
│   ├── check_comp.m      # Checks if data are valid compositional vectors (strictly positive).
│   ├── check_pw.m        # Checks if data are valid pairwise log-ratio vectors derived from compositional data.
│   └── check_sbp.m       # Validates sequential binary partition.
│
├── metric/       # Functions for computing distances and norms
│   ├── Ait_innerp.m      # Aitchison inner product
│   ├── Ait_norm.m        # Aitchison norm
│   ├── CoDa_distance.m   # Distance measures for CoDa
│   ├── L1CoDa.m          # L1 norm induced in the CoDa space
│   ├── L1clr_norm.m      # L1 norm in clr space
│   ├── L1plr_norm.m      # L1 norm of log-pairwises
│   ├── LinfCoDa_norm.m   # L∞ norm induced in the CoDa space
│   └── LpCoDa_norm.m     # Lp norm induced in the CoDa space
│
└── transform/    # Functions for CoDa transformations
    ├── alr_coord.m       # Additive log-ratio coordinates
    ├── clr_scores.m      # Centered log-ratio scores
    ├── mlr_scores.m      # Median centered log-ratio scores
    ├── nalr_coord.m      # nested additive log-ratio coordinates
    ├── olr_coord.m       # Orthonormal log-ratio coordinates
    ├── pw_matrix.m       # Pairwise matrix for CoDa
    ├── pv_basis.m        # Pivot basis
    └── sbp_basis.m       # Basis generated from a Sequential Binary Partition
```

---

## ✅ Purpose
The functions in this repository implement common operations for **Compositional Data Analysis**, including:
- **Validation**: Ensure data meets compositional constraints.
- **Metrics**: Compute distances and norms in CoDa space.
- **Transformations**: Apply log-ratio transformations (alr, clr, olr, etc.).

---

## 🔧 Installation
### Option 1: Install Toolbox (.mltbx)
1. Open MATLAB.
2. Go to **Home → Add-Ons → Install from File**.
3. Select `CoDaTools.mltbx`.
4. MATLAB will install and add it to your path automatically.

### Option 2: Manual Setup
```matlab
addpath(genpath('path_to_CoDaTools/src'));
```

---

## ✅ Quick Usage
```matlab
data = [0.2, 0.3, 0.5];
check_closure(data);      % Check closure property
clr_scores(data);         % Compute CLR scores
```

---

## 👤 Author
```matlab
Name: Jordi Saperas-Riera  
Email: jordi.saperas@udg.edu  
Institution: Universitat de Girona
```
---

## ⚠️ Notes
- Ensure input compositions are strictly positive. Each row represents an individual composition, and each column represents a part.


## 📚 References

### Books
- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman & Hall.
- Pawlowsky-Glahn, V., & Buccianti, A. (Eds.). (2011). *Compositional Data Analysis: Theory and Applications*. Wiley.
- Pawlowsky-Glahn, V., Egozcue, J. J., & Tolosana-Delgado, R. (2015). *Modeling and Analysis of Compositional Data*. Wiley.

### Papers
- Aitchison, J., Barceló-Vidal, C., Martín-Fernández, J., Pawlowsky-Glahn, V. (2000). Logratio analysis and compositional distance. *Mathematical Geology*, 32, 271–275.
- Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., Barceló-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology*, 35, 279–300.
- Barceló-Vidal, C., Martín-Fernández, J.A. (2016). The mathematics of compositional analysis. *Austrian Journal of Statistics*, 45, 57–71.
- Martín-Fernández, J.A. (2019). Comments on: compositional data: the sample space and its structure. *TEST*, 28, 653–657.
- Saperas-Riera, J. (2025). *Avenços en els fonaments matemàtics de l'anàlisi composicional de dades: convexitat i normes Lp. Aplicació a la regressió lineal LASSO amb covariable composicional* (Doctoral dissertation). Universitat de Girona. http://hdl.handle.net/10803/693962
- Saperas-Riera, J., Mateu-Figueras, G., & Martín-Fernández, J. A. (2023). Lasso regression method for a compositional covariate regularised by the norm L1 pairwise logratio. *Journal of Geochemical Exploration*, 255, 107327. https://doi.org/10.1016/j.gexplo.2023.107327
- Saperas-Riera, J., Mateu-Figueras, G., & Martín-Fernández, J. A. (2024). Lp-norm for compositional data: Exploring the CoDa L1-norm in penalised regression. *Mathematics*, 12(9), 1388. https://doi.org/10.3390/math12091388
