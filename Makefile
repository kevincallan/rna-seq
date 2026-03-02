# =============================================================================
# RNA-seq Pipeline Makefile
# =============================================================================
# Convenience targets for running the pipeline.
#
# Usage:
#   make run                    # Full pipeline with default config
#   make run SUBSET=day3_wt_vs_tet1   # Run with specific subset
#   make qc                     # Steps 0-4 only (QC)
#   make map                    # Steps 0-5 (through mapping)
#   make count                  # Steps 0-8 (through counting)
#   make report                 # Generate report (step 11 only)
#   make compare                # Compare methods (step 10 only)
#   make test                   # Run pytest
#   make clean                  # Remove results/ work/ logs/
# =============================================================================

CONFIG   ?= config/config.yaml
METHODS  ?=
SUBSET   ?=
PYTHON   ?= ./py

# Build optional flags
METHODS_FLAG = $(if $(METHODS),--methods $(METHODS),)
SUBSET_FLAG  = $(if $(SUBSET),--subset $(SUBSET),)

.PHONY: run qc map count de report compare test clean help

help:
	@echo "Available targets:"
	@echo "  run       Full pipeline"
	@echo "  qc        Steps 0-4 (validate + samples + trim + QC)"
	@echo "  map       Steps 0-5 (through STAR mapping)"
	@echo "  count     Steps 0-8 (through matrix filtering)"
	@echo "  de        Steps 0-9 (through DESeq2)"
	@echo "  report    Step 11 only (generate report)"
	@echo "  compare   Step 10 only (compare methods)"
	@echo "  test      Run pytest"
	@echo "  clean     Remove results/ work/ logs/"

run:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run $(METHODS_FLAG) $(SUBSET_FLAG)

qc:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run --steps 0 1 2 3 4 $(METHODS_FLAG) $(SUBSET_FLAG)

map:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run --steps 0 1 2 3 4 5 $(METHODS_FLAG) $(SUBSET_FLAG)

count:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run --steps 0 1 2 3 4 5 7 8 $(METHODS_FLAG) $(SUBSET_FLAG)

de:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run --steps 0 1 2 3 4 5 7 8 9 $(METHODS_FLAG) $(SUBSET_FLAG)

report:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run --steps 11 $(METHODS_FLAG)

compare:
	$(PYTHON) scripts/run_pipeline.py --config $(CONFIG) run --steps 10 $(METHODS_FLAG)

test:
	$(PYTHON) -m pytest tests/ -v

clean:
	rm -rf results/ work/ logs/
	@echo "Cleaned results/, work/, logs/"
