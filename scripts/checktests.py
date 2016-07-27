"""
Check the coverage of tests for Stan functions.  If determines which functions
are covered from

- names of the test-*.R
- functions listed in tets.yaml: since some models test multiple functions, especially the ones that test the extractor functions.

"""
import sys
import os
from os import path
import re

import yaml

from standoc import parse_stan_function_defs

def check_functions(filename, testdir):
    tests = []
    with open(path.join(testdir, 'tests.yaml'), 'r') as f:
        known_tests = yaml.load(f)
    for k, v in known_tests.items():
        for tst in v:
            tests.append(tst)
    for f in os.listdir(testdir):
        m = re.match(r'^test-(.*)\.R$', f)
        if m:
            tests.append(m.group(1))
    with open(filename, 'r') as f:
        functions = parse_stan_function_defs(f.read())
    missing = []
    for func in functions:
        if func not in tests:
            missing.append(func)
    if len(missing) == 0:
        print("All functions have tests!!!\n")
    else:
        print("Missing tests for %d functions" % len(missing))
        print("Functions missing tests are:")
        print('\n'.join("- %s" % x for x in missing))

def main():
    filename, testdir = sys.argv[1:3]
    check_functions(filename, testdir)

if __name__ == "__main__":
    main()
