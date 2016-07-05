CMDSTAN_VERSION = 2.10.0
CMDSTAN_URL = https://github.com/stan-dev/cmdstan/releases/download/v$(VERSION)/cmdstan-$(VERSION).tar.gz
CMDSTAN_DIR = cmdstan

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

test:
	make -C stan/tests test

cmdstan-build:
	wget $(CMDSTAN_URL)
	tar -xf $(notdir $(CMDSTAN_URL))
	mv $(patsubst %.tar.gz,%,$(notdir $(CMDSTAN_URL))) $(CMDSTAN_DIR)
	rm -f $(notdir $(CMDSTAN_URL))
	cd $(CDMSTAN_DIR) && make build
