# Makefile -----------------------------------------------------------------
#
# Convenience targets. Everything is plain Rscript underneath, so nothing here
# is required -- `Rscript scripts/03_leaderboard_reliability.R` works on its own.
#
#   make check     validate the G-theory engine against the paper's equations
#   make summary   design summary and coverage tables (seconds, no fits needed)
#   make analysis  every analysis script, from cached fits (minutes)
#   make fit       re-estimate all models (hours -- read the README first)
#   make all       fit, then analyse
#   make clean     remove generated figures, tables, and results (keeps fits)
#   make distclean also remove cached model fits

R := Rscript

.PHONY: all check summary analysis fit ablations methods clean distclean help

help:
	@sed -n '3,14p' Makefile | sed 's/^# \{0,1\}//'

check:
	$(R) tests/test_gtheory.R
	$(R) tests/test_ranks.R

summary:
	$(R) scripts/00_design_summary.R

# Re-estimation. Hours. See the timing table in README.md.
fit:
	REFIT=TRUE $(R) scripts/01_fit_models.R

analysis: summary
	$(R) scripts/02_benchmark_reliability.R
	$(R) scripts/03_leaderboard_reliability.R
	$(R) scripts/04_rank_analysis.R

# Frequentist ablation arms only. Add RUN_LOBO_BAYES=TRUE REFIT=TRUE for the
# Bayesian leave-one-benchmark-out refits (~57 h).
ablations:
	$(R) scripts/05_ablations.R

# Add RUN_DISCO=TRUE to include the nonparametric decomposition (~17 h).
methods:
	$(R) scripts/06_method_comparison.R

all: check fit analysis ablations methods

clean:
	rm -f outputs/figures/*.pdf outputs/tables/*.csv outputs/tables/*.tex outputs/results/*.csv

distclean: clean
	rm -f outputs/fits/*.rds
