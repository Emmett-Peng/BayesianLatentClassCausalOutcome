# BayesianLatentClassCausalOutcome

R code and simulation results for estimating causal effects on latent class membership with a fully Bayesian joint latent class model and posterior predictive standardization (Bayesian g-computation).

**Citation.** Peng Y, Lu Z, Anderson GM, Liu K. *Bayesian Latent Class Approach for Causal Estimation with Multiple Outcomes.* Manuscript submitted to *Statistical Methods in Medical Research*.

## Files

| File | Description |
|---|---|
| `simulation_functions.R` | All functions: data generation, the JAGS model, label switching, posterior predictive standardization, the comparator (two-step inverse probability weighted latent class estimator of Lanza et al.) and the summaries of MSE and coverage. Also the simulation design (scenario grid and seeds). Defines functions only. |
| `simulation_code.R` | Runs the simulation study: proposed model and comparator on 16 scenarios x 100 datasets, and writes the MSE and coverage of the 95% interval for each ATE_c. |
| `lanza_bootstrap_B1000.R` | Re-runs the comparator with a 1000-replicate bootstrap on exactly the same datasets. The comparator results in the paper come from this script. Does not need JAGS. |
| `run_lanza_bootstrap.sbatch` | Slurm array script for `lanza_bootstrap_B1000.R` (one task per dataset). |

## Simulation design

Each dataset has n = 2000 subjects, K = 4 latent classes and 6 continuous outcomes (two memory and four executive function measures), with three baseline confounders (education, sex, age), a binary exposure with prevalence 0.30 and a common outcome SD of 0.5. The four classes are high on both domains (HH), low on both (LL), high memory with low executive function (HL) and low memory with high executive function (LH).

The 16 scenarios cross three factors:

| Factor | Levels |
|---|---|
| Confounding strength (`weak`, `strong`) | scale 0.6 and 1.6 on the confounder coefficients of the exposure and class models |
| Class separation (`low`, `high`) | scale delta = 0.6 and 1.2 on the distance between the class means |
| Class imbalance | `balanced` (0.25, 0.25, 0.25, 0.25); `one_rare` (0.30, 0.30, 0.30, 0.10); `two_rare` (0.40, 0.40, 0.10, 0.10); `severe_one` (0.30, 0.35, 0.30, 0.05) |

The data-generating values are in `simulate_one()` in `simulation_functions.R` and in Tables S1 and S2 of the supplementary material.

Scenario ids are numbered 1 to 16 with confounding varying fastest, then imbalance, then separation (ids 1 to 8 are `low` separation and ids 9 to 16 are `high`). The seed of dataset `d` in scenario `r` is `seeds[(r - 1) * 100 + d]`, where `seeds` comes from `make_design()`. The seeds were originally drawn for a 24-scenario development grid that also contained a third, smaller separation level (delta = 0.3) which is not reported. `make_design()` draws the same seed vector and drops the unused first block, so every dataset uses the same seed as in the paper. The class and exposure intercepts are calibrated numerically (`optim`), so the generated datasets can differ slightly across R versions and platforms.

## Proposed model

The joint model (exposure, covariates, latent class and outcomes) is specified in JAGS (`JAGS_MODEL` in `simulation_functions.R`) and fitted with three chains of 4000 adaptation, 6000 burn-in and 2000 retained iterations (6000 posterior draws). Label switching is resolved after sampling by matching each draw's class means to the data-generating class means with the Hungarian algorithm (`clue::solve_LSAP`), which is possible in a simulation because the true class means are known. The causal effects ATE_c are obtained by posterior predictive standardization with 1000 synthetic covariate vectors per posterior draw. The 95% credible interval is the 2.5th and 97.5th percentile of the posterior draws of ATE_c.

## Requirements

R (4.5.0 was used on the cluster) and, for `simulation_code.R`, [JAGS](https://mcmc-jags.sourceforge.io/).

R packages: `rjags`, `coda`, `truncnorm`, `clue`, `flexmix`, `WeightIt`, `nnet`, `future.apply`.

## How to run

Run the scripts from the folder that contains them. Both accept environment variables (listed at the top of each script) for running a subset.

```bash
# small test: scenarios 1 and 9, 2 datasets each, 20 bootstrap replicates
SCENARIOS=1,9 N_DATASETS=2 B_BOOT=20 Rscript simulation_code.R

# full study (long: run on a cluster, with WORKERS parallel workers)
WORKERS=8 Rscript simulation_code.R

# comparator with B = 1000 (one Slurm array task per dataset), then aggregate
sbatch run_lanza_bootstrap.sbatch
AGGREGATE_ONLY=1 Rscript lanza_bootstrap_B1000.R
```

Every dataset is saved as its own checkpoint, so interrupted runs can be restarted and finished datasets are skipped.

## Notes on reproducibility

- The simulated datasets and the comparator bootstrap are seeded and reproducible. In `simulation_code.R`, the JAGS chains are not seeded explicitly, so the posterior summaries can differ slightly from the paper because of Monte Carlo variation and software versions.
- In `simulation_code.R` the comparator uses 50 bootstrap replicates. The comparator intervals in the paper use 1000 replicates (`lanza_bootstrap_B1000.R`).
