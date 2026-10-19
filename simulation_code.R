################################################################################
# simulation_code.R
#
# Simulation study: proposed Bayesian joint latent class causal model versus the
# two-step inverse probability weighted latent class estimator (comparator).
#
# For each of the 16 scenarios (confounding x class imbalance x class separation)
# and each of 100 datasets, the script simulates data, fits the proposed model
# (JAGS), computes the posterior distribution of ATE_c by posterior predictive
# standardization, fits the comparator with a bootstrap, and records the MSE
# and 95% interval coverage of both methods.
#
# Usage
#   Rscript simulation_code.R                       # all 16 scenarios x 100 datasets
#   SCENARIOS=1,9 N_DATASETS=2 Rscript simulation_code.R   # small test
#   AGGREGATE_ONLY=1 Rscript simulation_code.R      # only (re)build the summary files
#
# Environment variables (all optional)
#   SCENARIOS       scenario ids to run, e.g. "1,3,5-8"            (default: all 16)
#   N_DATASETS      datasets per scenario                          (default: 100)
#   B_BOOT          bootstrap replicates of the comparator         (default: 50)
#   M_SYNTH         synthetic covariate vectors per posterior draw (default: 1000)
#   WORKERS         parallel workers over datasets                 (default: 1)
#   OUT_DIR         output directory                               (default: output)
#   AGGREGATE_ONLY  1 to skip the simulation and only aggregate    (default: 0)
#   FORCE           1 to recompute datasets that already have a checkpoint
#
# Each (scenario, dataset) pair is saved as its own checkpoint, so an interrupted
# run can simply be restarted. The comparator intervals reported in the paper come
# from lanza_bootstrap_B1000.R (B = 1000), run on the same datasets.
#
# Outputs (in OUT_DIR)
#   final_simulation_results.csv   MSE and coverage by scenario, method and class
#   scenario_index.csv             scenario settings, achieved exposure and class proportions
#   checkpoints/                   one .rds per dataset
################################################################################

source("simulation_functions.R")

## settings -------------------------------------------------------------------
B_BOOT         <- as.integer(Sys.getenv("B_BOOT", "50"))
M_SYNTH        <- as.integer(Sys.getenv("M_SYNTH", "1000"))
N_DATASETS     <- as.integer(Sys.getenv("N_DATASETS", SIM$D_per_scenario))
OUT_DIR        <- Sys.getenv("OUT_DIR", "output")
AGGREGATE_ONLY <- Sys.getenv("AGGREGATE_ONLY", "0") == "1"
FORCE          <- Sys.getenv("FORCE", "0") == "1"
WORKERS        <- get_workers()

design    <- make_design()
scenarios <- design$scenarios
seeds     <- design$seeds
scn_ids   <- parse_ids(Sys.getenv("SCENARIOS", ""), scenarios$scenario_id)
ds_ids    <- seq_len(N_DATASETS)

ckpt_dir  <- file.path(OUT_DIR, "checkpoints")
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
ckpt_file <- function(r, d) file.path(ckpt_dir, sprintf("rep_scn%02d_d%03d.rds", r, d))

message("simulation_code.R: scenarios = ", paste(scn_ids, collapse = ","),
        "; datasets = ", N_DATASETS, "; B_BOOT = ", B_BOOT, "; workers = ", WORKERS)

## simulation -----------------------------------------------------------------
if (WORKERS > 1L) plan(multisession, workers = WORKERS) else plan(sequential)

if (!AGGREGATE_ONLY) {
  for (r in scn_ids) {
    scn  <- scenarios[r, ]
    base <- (r - 1L) * SIM$D_per_scenario
    message("[", format(Sys.time(), "%H:%M:%S"), "] scenario ", r, "/", nrow(scenarios),
            "  ", scenario_tag(scn))

    todo <- if (FORCE) ds_ids else ds_ids[!file.exists(vapply(ds_ids, function(d) ckpt_file(r, d), ""))]
    if (!length(todo)) next

    invisible(future_lapply(todo, function(d) {
      res <- tryCatch(
        run_one_in_scenario(scn, seed = seeds[base + d], M_synth = M_SYNTH, B_boot = B_BOOT),
        error = function(e) list(error = conditionMessage(e)))
      saveRDS(res, ckpt_file(r, d))
      NULL
    }, future.seed = TRUE))
  }
}
plan(sequential)

## aggregation ------------------------------------------------------------------
message("Aggregating results ...")
summary_list <- list(); index_list <- list()
for (r in scn_ids) {
  scn <- scenarios[r, ]
  files <- ckpt_file(r, ds_ids)
  files <- files[file.exists(files)]
  resL  <- lapply(files, readRDS)
  failed <- vapply(resL, function(x) !is.null(x$error), logical(1))
  if (any(failed)) message("  scenario ", r, ": ", sum(failed), " dataset(s) failed and are excluded")
  resL <- resL[!failed]
  if (!length(resL)) { message("  scenario ", r, ": no results"); next }

  s <- summarize_scenario(scn, resL)
  s$n_datasets <- length(resL)
  summary_list[[length(summary_list) + 1L]] <- s
  index_list[[length(index_list) + 1L]]     <- scenario_index_row(scn, resL)
}

write.csv(do.call(rbind, summary_list), file.path(OUT_DIR, "final_simulation_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, index_list),   file.path(OUT_DIR, "scenario_index.csv"),           row.names = FALSE)
message("Wrote final_simulation_results.csv and scenario_index.csv in ", OUT_DIR)
message("Done.")
