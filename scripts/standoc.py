#!/usr/bin/env %python
import re
import sys
import yaml

import argparse

""" Regex matching the start of function doc block /** function_name """
RE_FUNCTION_BLOCK_START = re.compile(r"^\s*/\*\*\s*function")
""" Regex matching the start of a non-function doc block /** """
RE_BLOCK_START = re.compile(r"^\s*/\*\*")
""" Regex matching the end of a doc block """
RE_BLOCK_END = re.compile(r"^\s*\*/")

_CHAPTER = "# "
_SECTION = "## "
_SUBSECTION = "### "
_FUNCTION = "### "

def text_or_list(x):
    if isinstance(x, list):
        return ''.join(x)
    else:
        return str(x)

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
    template = '```stan\n{text}\n```'

    def __init__(self, text):
        self.text = text_or_list(text)
        self.functions = parse_stan_function_defs(self.text)

    def format(self):
        return self.template.format(text = self.text)

    def is_empty(self):
        return self.text.strip() == ''

class DocBlock:
    """ CodeBlock class

    The code block class contains blocks of comments
    """
    _re_leading_comment = re.compile(r'\s*\*')
    template = "{text}"

    def __init__(self, text):
        self.text = text_or_list(text)

    # def append(self, line):
    #     line = self._re_leading_comment.sub('', line).rstrip()
    #     is_tag = process_tag(line)
    #     if is_tag:
    #         if is_tag[0] == 'chapter':
    #             self.lines.append(_CHAPTER + " " + is_tag[1])
    #         elif is_tag[0] == 'section':
    #             self.lines.append(_SECTION + " " + is_tag[1])
    #         elif is_tag[0] == 'subsection':
    #             self.lines.append(_SUBSECTION + " " + is_tag[1])
    #         else:
    #             print("tag not recognized:", line)
    #             self.lines.append(line)
    #     else:
    #         self.lines.append(line)

    def format(self):
        return self.template.format(text = self.text)

    def is_empty(self):
        return self.text.strip() == ''

class FunctionDocBlock(DocBlock):
    """ CodeBlock class

    The code block class contains blocks of comments
    """
    _re_leading_comment = re.compile(r'\s*\*')
    template = """
### {funname}

**Parameters:**

{params}

**Return Value:** {rtrn}

{body}
"""

    def __init__(self, function, text):
        self.text = text_or_list(text)
        data = yaml.loads(self.text)

    def format(self):
        return self.text

    # def format(self):
    #     if len(self.params) > 0:
    #         param_list = '\n'.join(["- `{type}`  **{name}** {description}".format(**x)
    #                                 for x in self.params])
    #     else:
    #         param_list = ""
    #     if self.return_type:
    #         try:
    #             rtrn = "`{0}` {1}".format(*self.return_type)
    #         except:
    #             print("WARNING: %s has bad return type" % self.name)
    #             print(self.return_type)
    #             rtrn = ""
    #     else:
    #         print("WARNING: %s has no return type" % self.name)
    #         rtrn = ""
    #     msg = template.format(funname = self.name,
    #                            params = param_list,
    #                            rtrn = rtrn,
    #                            body = '\n'.join(self.lines))
    #     return msg

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

    def format(self):
        txt = '\n'.join([x.format() for x in self.data if not x.is_empty()])
        if not txt.endswith('\n'):
            txt += '\n'
        return txt

def process_tag(x):
    m = re.match('\s*@([A-Za-z][A-Za-z0-9_]*)\s*(.*)', x)
    if m:
        tag, text = m.groups()
        return (tag, text.strip())
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
    sink = []
    for linenum, line in enumerate(text):
        print("parsing line %d" % linenum)
        if current_block is None:
            if is_function_docblock_start(line):
                current_block = 'function'
            elif is_docblock_start(line):
                current_block = 'doc'
            else:
                current_block = 'code'
                sink.append(line)

        elif current_block in ('function', 'doc'):
            if is_docblock_end(line):
                if current_block == 'function':
                    doc.append(FunctionDocBlock(sink))
                else:
                    doc.append(DocBlock(sink))
                current_block = None
                sink = []
            else:
                sink.append(line)

        elif current_block in ('code',):
            if is_function_docblock_start(line):
                doc.append(CodeBlock(sink))
                sink = []
                current_block = 'function'
            elif is_docblock_start(line):
                doc.append(CodeBlock(sink))
                sink = []
                current_block = 'doc'
            else:
                sink.append(line)
        else:
            print("bad line: %s" % line)
    if current_block == 'function':
        doc.append(FunctionDocBlock(sink))
    elif current_block == 'doc':
        doc.append(DocBlock(sink))
    elif current_block == 'code':
        doc.append(CodeBlock(sink))
    else:
        print("Somethings wrong, document ended in state %s" % current_block)
    return doc

def create_docfile(src, dst):
    docs = parse(src).format()
    dst.write(docs)

def main():
    parser = argparse.ArgumentParser(description = "Parse a Stan file and a output markdown formatted file")
    parser.add_argument('infile', nargs = '?', type = argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('outfile', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
    args = parser.parse_args()
    create_docfile(args.infile, args.outfile)

if __name__ == "__main__":
    main()
