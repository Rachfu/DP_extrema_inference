# DP-Extrema-inference

This repository is the official implementation of paper "**Differentially Private Inference for Extrema of
Parameters in Clinical Studies**".

# Requirements

To install requirements:
```
packages_to_install <- readLines("requirements.txt")
install.packages(packages_to_install)
```

# Simulation
```
cd simulation
```
*fun_data.R*: generate data for simulation, including Gaussian case, linear regression case, etc.
*fun_model.R*: models including our proposed method.
*utils.R*: functions for general use, such as generating Laplacian noises.
*run.R*: run this file to do simulation, including Gaussian case, linear regression case, etc.

# Real data application
```
cd real_application
```
we apply our proposed method to AIDS Clinical Trials Group Study 175 (ACTG175). ACTG 175 was a randomized clinical trial to compare monotherapy with zidovudine or didanosine with combination therapy with zidovudine and didanosine or zidovudine and zalcitabine in adults infected with the human immunodeficiency virus type I and whose CD4 T cell counts were between 200 and 500 per cubic millimeter. We explore the treatment effect of zidovudine-incorporated therapy with a binary treatment indicator and consider CD8 T cell count at 20+/-5 weeks (cd820) as the response vector and 15 pre-treatment covariates in regression adjustment, including age, weight in kilograms (wtkg), race, Karnofsky score on a scale of 0-100 (karnof) and etc.

*fun_data.R*: two scenarios to explore the corresponding best subgroup effect with the linear regression models respectively. The first scenario consists of four subpopulations with disciplines of race (white or non-white) and age (below or above 34), while the second scneraio consists of four subpopulations with disciplines of symptom (symptomatic indicator equals 0 meaning asymptomatic and vice versa) and offtrt (indicator of off-treatment before 96\%pm5 weeks equals 0 and vice versa).
*fun_model.R*: models including our proposed method. Parallel computation is implemented in this real application to improve efficiency.
*utils.R*: functions for general use, such as generating Laplacian noises.
*run.R*: run this file to get experiment results with real data application, including the two scenarios of subgroup splitting.

