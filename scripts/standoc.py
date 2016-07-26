#!/usr/bin/env %python
import re
import sys
import argparse

""" Regex matching the start of function doc block /** function_name """
RE_FUNCTION_BLOCK_START = re.compile(r"^\s*/\*\*+\s*([A-Za-z][A-Za-z0-9_]*)\b")
""" Regex matching the start of a non-function doc block /** """
RE_BLOCK_START = re.compile(r"^\s*/\*\*")
""" Regex matching the end of a doc block """
RE_BLOCK_END = re.compile(r"^\s*\*/")

_CHAPTER = "# "
_SECTION = "## "
_SUBSECTION = "### "
_FUNCTION = "### "

def parse_stan_function_defs(text):
    text = re.sub('\s+', ' ', text, re.DOTALL + re.M)
    basic_types = ('int', 'real', 'matrix', 'vector', 'row_vector')
    brackets = (r'\[(?:\s|,)*\]')
    identifier = r'[A-Za-z][A-Za-z0-9_]*'
    return_types =  r'(?:void|(?:{basic_types})\s*(?:{brackets})?)'.\
        format(basic_types = '|'.join(basic_types), brackets = brackets)
    arg_type = r'(?:{basic_types})\s*(?:{brackets})?'.\
        format(basic_types = '|'.join(basic_types), brackets = brackets)
    arg = r'{type}\s+{name}'.format(type = arg_type, name = identifier)
    arglist = r'({arg})(?:\s*,\s*({arg}))*'.format(arg = arg)
    function_def = r''.join((r'^\s*(?P<return>{return_type})',
                             r'\s+(?P<func>{function_name})',
                             r'\s*\(\s*(?P<arglist>{arglist})?\s*\)')).\
        format(return_type = return_types, function_name = identifier,
               arglist = arglist)
    functions = {}
    for fun_def in re.finditer(function_def, text, re.M):
        newfun = {'return_type': fun_def.group('return'),
                  'args': []}
        for arg in re.match(arglist, fun_def.group('arglist')).groups():
            if arg:
                argtype, argname = re.split('\s+', arg)
                newfun['args'].append({'type': argtype, 'name': argname})
        functions[fun_def.group('func')] = newfun
    return functions

class CodeBlock:
    """ CodeBlock class

    The code block class contains blocks of code
    """
    before = ""
    after = ""

    def __init__(self, lines = []):
        self.lines = []
        if not isinstance(lines, list):
            lines = list(lines)
        for ln in lines:
            self.append(ln)

    def append(self, line):
        line = line.rstrip()
        self.lines.append(line)

    def to_string(self):
        return self.before + '\n'.join(self.lines) + self.after

    def is_empty(self):
        return ''.join(self.lines).strip() == ''

class DocBlock:
    """ CodeBlock class

    The code block class contains blocks of comments
    """
    _re_leading_comment = re.compile(r'\s*\*')
    before = "/**\n"
    after = "\n*/"

    def __init__(self, lines = []):
        self.lines = []
        if not isinstance(lines, list):
            lines = list(lines)
        for ln in lines:
            self.append(ln)

    def append(self, line):
        line = self._re_leading_comment.sub('', line).rstrip()
        self.lines.append(line)

    def to_string(self):
        return self.before + '\n'.join(self.lines) + self.after

    def is_empty(self):
        return ''.join(self.lines).strip() == ''

class FunctionDocBlock(DocBlock):
    """ CodeBlock class

    The code block class contains blocks of comments
    """
    _re_leading_comment = re.compile(r'\s*\*')

    def __init__(self, funcname, lines = []):
        self.lines = []
        self.params = []
        self.name = funcname
        self.return_type = None
        if not isinstance(lines, list):
            lines = list(lines)
        for ln in lines:
            self.append(ln)

    def append(self, line):
        line = self._re_leading_comment.sub('', line).rstrip()
        is_tag = process_tag(line)
        if is_tag:
            if is_tag[0] == 'param':
                self.params.append({'type': is_tag[2][0],
                                    'name': is_tag[2][1],
                                    'description': is_tag[2][2]})
            elif is_tag[0] == 'return':
                self.return_type = is_tag[2]
            else:
                print("tag not recognized:", line)
                self.lines.append(line)
        else:
            self.lines.append(line)

    def to_string(self):
        template = """/**
---
name: {funname}
param:
{params}
return: {rtrn}
---

{body}
*/
"""
        if len(self.params) > 0:
            param_list = '\n'.join(["- name: {name}\n  description: {description}".format(**x)
                                    for x in self.params])
        else:
            param_list = ""
        if self.return_type:
            try:
                rtrn = "{1}".format(*self.return_type)
            except:
                print("WARNING: %s has bad return type" % self.name)
                print(self.return_type)
                rtrn = ""
        else:
            print("WARNING: %s has no return type" % self.name)
            rtrn = ""
        msg = template.format(funname = self.name,
                               params = param_list,
                               rtrn = rtrn,
                               body = '\n'.join(self.lines))
        return msg

    def is_empty(self):
        return False

class Document(object):
    # The Document class is a list of CodeBlocks and DocBlocks

    def __init__(self, data = []):
        if not isinstance(data, list):
            data = list(data)
        self.data = data

    def append(self, block):
        self.data.append(block)

    def number_blocks(self):
        """ Number of blocks in the Document """
        return len(self.data)

    def to_string(self):
        txt = '\n'.join([x.to_string() for x in self.data if not x.is_empty()])
        return txt


def process_tag(x):
    _TAGS = {
        'param' : lambda x: re.split(r'\s+', x.strip(), 2),
        'return' : lambda x: re.split(r'\s+', x.strip(), 1),
        'section': lambda x: x.strip(),
        'subsection': lambda x: x.strip()
    }
    m = re.match('\s*@([A-Za-z][A-Za-z0-9_]*)\s*(.*)', x)
    if m:
        tag, text = m.groups()
        if tag in _TAGS:
            parsed = _TAGS[tag](text)
        else:
            parsed = None
        return (tag, text, parsed)
    else:
        return None

def is_function_docblock_start(x):
    """ Check if line starts a function block
    """
    return RE_FUNCTION_BLOCK_START.search(x)

def is_docblock_start(x):
    """ Check if comment starts a doc block """
    return RE_BLOCK_START.search(x)

def is_docblock_end(x):
    """ Check if comment ends a doc block """
    return RE_BLOCK_END.search(x)

def parse(f):
    """Parse a Stan File"""
    text = f.readlines()
    current_block = None
    doc = Document()
    for linenum, line in enumerate(text):
        m_function_block = is_function_docblock_start(line)
        if current_block is None:
            if m_function_block:
                current_block = FunctionDocBlock(m_function_block.group(1))
            elif is_docblock_start(line):
                current_block = DocBlock()
            else:
                current_block = CodeBlock([line])

        elif (isinstance(current_block, DocBlock) or
            isinstance(current_block, FunctionDocBlock)):
            if is_docblock_end(line):
                doc.append(current_block)
                current_block = CodeBlock()
            else:
                current_block.append(line)

        elif isinstance(current_block, CodeBlock):
            if m_function_block:
                doc.append(current_block)
                current_block = FunctionDocBlock(m_function_block.group(1))
            elif is_docblock_start(line):
                doc.append(current_block)
                current_block = DocBlock()
            else:
                current_block.append(line)
        else:
            print("bad line: %s" % line)
    doc.append(current_block)
    return doc

def create_docfile(src, dst):
    docs = parse(src).to_string()
    dst.write(docs)

def main():
    parser = argparse.ArgumentParser(description = "Parse a Stan file and a output markdown formatted file")
    parser.add_argument('infile', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('outfile', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
    args = parser.parse_args()
    create_docfile(args.infile, args.outfile)

if __name__ == "__main__":
    main()
