# pvca
Principal Variance Component Analysis using lme4 implemented in Python

Notes: 

installation has to be done through conda-forge due to dependencies on pymer4, follow the instructions on pymer4. r-matrix has some incompatibility breaks that aren't accounted for, so you need to specify `r-matrix=1.6_1`. 
```
conda create --name pvca -c ejolly -c conda-forge -c defaults pymer4 tqdm scipy pandas r-matrix=1.6_1
```