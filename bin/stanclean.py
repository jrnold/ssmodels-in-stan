#!/usr/bin/env %python
import re
import sys
import argparse

""" Regex matching the start of a non-function doc block /** """
RE_BLOCK_START = re.compile(r"^\s*/\*\*")
""" Regex matching the end of a doc block """
RE_BLOCK_END = re.compile(r"^\s*\*/")

def remove_line_comments(x):
    return re.sub('(#|//).*$', '', x)

def clean(f):
    """Remove comments from a Stan file"""
    text = f.readlines()
    cleantext = []
    in_comment = False
    for linenum, line in enumerate(text):
        if in_comment:
            if RE_BLOCK_END.match(line):
                in_comment = False
        else:
            if RE_BLOCK_START.match(line):
                in_comment = True
            else:
                line = remove_line_comments(line).rstrip()
                if line != "":
                    cleantext.append(line)
    return '\n'.join(cleantext) + '\n'

def remove_comments(src, dst):
    docs = clean(src)
    dst.write(docs)

def main():
    parser = argparse.ArgumentParser(description = "Remove comments from a Stan file")
    parser.add_argument('infile', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('outfile', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
    args = parser.parse_args()
    remove_comments(args.infile, args.outfile)

if __name__ == "__main__":
    main()
