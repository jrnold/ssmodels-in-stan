STAN = "stan/includes/ssm.stan"
PKGDIR =

build: doc clean

doc: doc/stanfunctions.Rmd

clean: StanStateSpace/inst/stan/ssm.stan

doc/stanfunctions.Rmd: scripts/standoc.py stan/includes/ssm.stan
	python $^ $@

StanStateSpace/inst/stan/ssm.stan: scripts/stanclean.py stan/includes/ssm.stan
	python $^ $@

	
