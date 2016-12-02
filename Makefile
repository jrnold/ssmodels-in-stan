CMDSTAN_VERSION = 2.12.0
CMDSTAN_URL = https://github.com/stan-dev/cmdstan/releases/download/v$(CMDSTAN_VERSION)/cmdstan-$(CMDSTAN_VERSION).tar.gz
CMDSTAN_DIR = lib/cmdstan-$(CMDSTAN_VERSION)


STAN = stan/includes/ssm.stan
PKG_DIR = StanStateSpace
STAN_INCLUDE_DIR = $(PKG_DIR)/inst/stan/include/
SSM_STAN = stan/ssm.stan

SSM_STAN_CLEAN = $(STAN_INCLUDE_DIR)/$(notdir $(SSM_STAN))
STAN_FUNCTION_DOC = doc/stanfunctions.Rmd

TEST_DIR = tests
TESTS = $(patsubst tests/src/%.stan,tests/build/%,$(wildcard tests/src/*.stan))


build: clean doc

doc: $(STAN_FUNCTION_DOC)

clean: $(SSM_STAN_CLEAN)

$(STAN_FUNCTION_DOC): ./bin/standoc.py $(SSM_STAN)
	python $^ $@

$(SSM_STAN_CLEAN): ./bin/stanclean.py $(SSM_STAN)
	python $^ $@

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
	#watchman -- trigger . tests $(SSM_STAN_CLEAN) 'stan/*.stan' 'stan/tests/src/*.stan' 'stan/tests/Makefile' 'stan/tests/*.R' -- make test

unwatch:
	watchman watch-del .

tests: $(TESTS)
  -Rscript -e 'testthat::test_dir("$(TEST_DIR)")'

test-%: $(TEST_DIR)/build/test_%

$(TEST_DIR)/build/%.stan: $(TEST_DIR)/src/%.stan $(SSM_STAN)
	mkdir -p build
	# hacky #include alternative
	# use _ as alterate separator
	# sed r command includes the contents of a file
	sed -e  '\_#include *ssm.stan_r$(SSM_STAN)' < $< > $@

$(TEST_DIR)/build/%: $(TEST_DIR)/build/%.stan $(TEST_DIR)/src/%.stan $(SMM_STAN)
	-make -C $(CMDSTAN_DIR) $(abspath $@)
