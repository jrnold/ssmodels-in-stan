CMDSTAN_VERSION = 2.10.0
CMDSTAN_URL = https://github.com/stan-dev/cmdstan/releases/download/v$(CMDSTAN_VERSION)/cmdstan-$(CMDSTAN_VERSION).tar.gz
CMDSTAN_DIR = lib/cmdstan-$(CMDSTAN_VERSION)

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

test: build
	make -C stan/tests test

cmdstan-build:
	mkdir -p $(notdir $(CMDSTAN_DIR))
	wget -nc $(CMDSTAN_URL)
	tar -xf $(notdir $(CMDSTAN_URL))
	mv $(patsubst %.tar.gz,%,$(notdir $(CMDSTAN_URL))) $(CMDSTAN_DIR)
	rm -f $(notdir $(CMDSTAN_URL))
	make -C $(CMDSTAN_DIR) build

# Watchman triggers to automatically build stuff
watch:
	watchman watch .
	watchman -- trigger . buildssm 'stan/*.stan' -- make build
	watchman -- trigger . tests $(SSM_STAN_CLEAN) 'stan/*.stan' 'stan/tests/src/*.stan' 'stan/tests/Makefile' 'stan/tests/*.R' -- make test

unwatch:
	watchman watch-del .
