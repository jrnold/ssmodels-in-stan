import sys
import re
outdir = "tests"
with open("test-ssm.R", "r") as f:
    lines = f.readlines()

re_testthat = re.compile(r"\s*test_that\(")
re_end = re.compile(r"^}\)")

current = []
in_test = False
unused_lines = []
for line in lines:
    if in_test:
        if re_end.match(line):
            current.append(line)
            print("end of %s" % function_name)
            # write it out
            with open('newtests/%s.R' % function_name, 'w') as f:
                f.write(''.join(current))
            current = []
            in_test = False
            function_name = None
        else:
            current.append(line)
    else:
        if re_testthat.match(line):
            try:
                function_name = re.search(r"stan\s*functions?\s*(\w+)", line, re.I).group(1)
            except AttributeError as e:
                print(line)
                raise e
            if function_name == "s":
                print(line)
                sys.exit(1)
            print(function_name)
            current.append(line)
            in_test = True
        else:
            unused_lines.append(line)

with open('unused.txt', 'w') as f:
    f.write(''.join(unused_lines))
