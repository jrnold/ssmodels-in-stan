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
import fnmatch

import yaml

from standoc import parse_stan_function_defs

def check_functions(filename, testdir):
    covered = set()
    for fname in os.listdir(testdir):
        if fnmatch.fnmatch(fname, 'test-*.R'):
            with open(path.join(testdir, fname), 'r') as f:
                for line in f:
                    m = re.match(r'#function:\s+(\w+)', line)
                    if m:
                        covered.add(m.group(1))
    with open(filename, 'r') as f:
        functions = set(parse_stan_function_defs(f.read()))
    missing = functions - covered
    if len(missing) == 0:
        print("All functions have tests!!!\n")
    else:
        print("Missing tests for %d functions:" % len(missing))
        print('\n'.join("- %s" % x for x in missing))
    unknown = covered - functions
    if len(unknown) > 0:
        print("\nThere are some unknown functions in tests")
        print('\n'.join("- %s" % x for x in unknown))

def main():
    filename, testdir = sys.argv[1:3]
    check_functions(filename, testdir)

if __name__ == "__main__":
    main()
