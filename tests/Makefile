CMDSTAN_DIR = ../../lib/cmdstan-2.12.0
SSM_STAN = ../../StanStateSpace/inst/stan/include/ssm.stan
STAN_TESTS = $(patsubst src/%.stan,build/%,$(wildcard src/*.stan))

all: test

test: $(STAN_TESTS)
	Rscript -e 'testthat::test_dir(".")'

test-%: /build/test_%

build/%.stan: src/%.stan $(SSM_STAN)
	mkdir -p build
	# hacky #include alternative
	# use _ as alterate separator
	# sed r command includes the contents of a file
	sed -e  '\_#include *ssm.stan_r$(SSM_STAN)' < $< > $@

build/%: build/%.stan src/%.stan $(SMM_STAN)
	make -C $(CMDSTAN_DIR) $(abspath $@)
