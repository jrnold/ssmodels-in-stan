STAN = stan/includes/ssm.stan
PKG_DIR = StanStateSpace
STAN_INCLUDE_DIR = $(PKG_DIR)/inst/stan/include/
SSM_STAN = stan/ssm.stan
SSM_STAN_CLEAN = $(STAN_INCLUDE_DIR)/$(notdir $(SSM_STAN))
STAN_FUNCTION_DOC = doc/stanfunctions.Rmd

build: doc clean

doc: $(STAN_FUNCTION_DOC)

clean: $(SSM_STAN_CLEAN)

$(STAN_FUNCTION_DOC): scripts/standoc.py $(SSM_STAN)
	python $^ $@

$(SSM_STAN_CLEAN): scripts/stanclean.py $(SSM_STAN)
	python $^ $@
