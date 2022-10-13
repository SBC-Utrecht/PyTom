'''
Created on Jun 29, 2011

@author: yuxiangchen
'''

import sys
from pytom.tools.files import checkFileExists, checkDirExists

def parse_script_options(args, helper):
    """
    parse script options
    @param args: array of input arguments (typically from command line)
    @type args: L{list}
    @param helper: script helper
    @type helper: L{pytom.tools.script_helper.ScriptHelper}
    @return: result of parsing
    @rtype: L{list}
    """
    import getopt

    if '--help' in args or '-h' in args:
        print(helper)
        sys.exit()
    
    res = [] # the result of parsing
    opt_str = ""
    long_opt = []
    for opt in helper.options:
        res.append(None)
        for name in opt.option_str:
            if name[:2] == "--":
                if opt.arg:
                    long_opt.append(name[2:]+'=')
                else:
                    long_opt.append(name[2:])
            elif name[:1] == "-":
                if opt.arg:
                    opt_str += name[1]+":"
                else:
                    opt_str += name[1]
            else:
                raise Exception("Option format invalid: %s" % name)
    
    try:
        opts, args = getopt.getopt(args, opt_str, long_opt)
    except getopt.GetoptError:
        raise
    
    for o,a in opts:
        for i in range(len(helper.options)):
            if o in helper.options[i].option_str:
                if helper.options[i].arg:
                    res[i] = a
                else:
                    res[i] = True
                break

    for i,_ in enumerate(res):
        if res[i] is None and not helper.options[i].optional:
            raise Exception("Required flag not passed, use any of the following: " + " ".join(helper.options[i].option_str))

    return res


def parse_script_options2(args, helper):
    """
    This function will parse given arguments based on a helper, it will generate nice error messages for every error it
    finds and is guarenteed to handel all internally thrown exceptions (so it cannot throw exceptions, it will stop
    execution instead). It will also parse the arguments into the right dataformat if it is provided.
    parse script options
    @param args: array of input arguments (typically from command line)
    @type args: L{list}
    @param helper: script helper
    @type helper: L{pytom.tools.script_helper.ScriptHelper}
    @return: result of parsing
    @rtype: L{list}
    """
    # To buffer all exceptions to present them in one batch to the user so all errors can be resolved in one go
    exception = ""

    res = []  # the result of parsing
    names = []
    arguments = []

    try:
        if '--help' in args:
            import sys
            print(helper)
            sys.exit()

        skip = False
        for n in range(len(args)):
            match = False

            if skip:
                skip = False
                continue

            for opt in helper.options:
                for name in opt.option_str:
                    if name == args[n]:
                        if opt.arg:
                            if n + 1 >= len(args):
                                exception += "Option {:s} requires arguments but none are given\n".format(name)
                                continue

                            a, e = parse_argument(args[n + 1], opt.arg_str, name)
                            if e != "": exception += e + "\n"
                            names.append(name)
                            arguments.append(a)
                            skip = True
                            match = True
                        else:
                            names.append(name)
                            arguments.append(True)
                            match = True

            if not match:
                exception += 'Could not match argument {:s} to any option\n'.format(args[n])

        if len(names) == 0:
            import sys
            print(helper)
            sys.exit()

        # Generate the resultslist in the right order
        for opt in helper.options:
            match = False
            for name in opt.option_str:
                try:
                    i = names.index(name)
                    if match: exception += 'Option {:s} was already defined (possibly as a different option string), multiple definitions for a single option is not supported\n'.format(
                        name)
                    match = True
                    if isinstance(arguments[i], list) and len(arguments[i]) == 1:
                        res.append(arguments[i][0])
                    else:
                        res.append(arguments[i])
                except:
                    pass

            if not match:
                if opt.required:
                    exception += "Option {:s} requires arguments but none are given\n".format(opt.option_str[0])
                    res.append(opt.default)
                else:
                    res.append(opt.default)

    except Exception as e:
        exception += "EXCEPTION:" + str(e) + "\n"

    if exception != "":
        import sys
        color1 = ""
        color2 = ""

        if sys.stdout.isatty():
            color1 = "\033[91m"
            color2 = "\033[0m"

        print(
                    "Some exception(s) while parsing the command line arguments:\n" + color1 + exception + color2 + "\nUse '--help' to get information on how this script can be called.")
        sys.exit()

    return res[:-1]  # Because the last option is by definition the help option


def parse_argument(inp, arg_str, option):
    """
    Parses the given input into the format as defined by the arg_str
    :param inp: The input string
    :param arg_str: The argument string
    :param option: The name of the option (used for creating error messages)
    :return: (result: list of values in the format specified, error message: str)
    :author: Douwe Schulte
    """
    import re
    if inp == "":
        return None, "The argument is zerolength for {:s}".format(option)
    if inp is None:
        return None, "The argument is None for {:s}".format(option)

    if arg_str == 'has arguments':
        return inp, ""
    else:
        res = []
        arg_str = re.escape(arg_str)
        reg = "^" + arg_str + "$"
        # construct a regex based on the arg_str
        for r in (("uint", r"(\d+)"), ("int", r"([-+]?\d+)"), ("float", r"([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]\d+)?)"),
                  ("string", r"(\S+?)"), ("file", r"(\S+?)"), ("directory", r"(\S+?)"),):
            reg = reg.replace(*r)
        match = re.match(reg, inp)

        if match:
            # retrieve the substrings and try to read them in
            substrings = re.findall("int|float|str|directory|file", arg_str)
            for n in range(len(substrings)):
                sub = match.group(n + 1)  # 0 is the full length match

                if substrings[n] == "uint":
                    try:
                        res.append(int(sub))
                    except:
                        return None, "The argument is not in the right format, there should be an uint at this position but it is not valid, for option {:s} and argument string \"{:s}\" with input {:s}".format(
                            option, arg_str, inp)
                elif substrings[n] == "int":
                    try:
                        res.append(int(sub))
                    except:
                        return None, "The argument is not in the right format, there should be an int at this position but it is not valid, for option {:s} and argument string \"{:s}\" with input {:s}".format(
                            option, arg_str, inp)
                elif substrings[n] == "float":
                    try:
                        res.append(float(sub))
                    except:
                        return None, "The argument is not in the right format, there should be a float at this position but it is not valid, for option {:s} and argument string \"{:s}\" with input {:s}".format(
                            option, arg_str, inp)

                elif substrings[n] == 'file':
                    if not checkFileExists(sub):
                        return None, f'File {sub} does not exist!'
                    else:
                        res.append(sub)

                elif substrings[n] == 'directory':
                    if not checkDirExists(sub):
                        return None, f'Directory {sub} does not exist!'
                    else:
                        res.append(sub)
                else:
                    res.append(sub)

            return res, ""

        return None, "The input does not match the given format, for option {:s} and argument string \"{:s}\" with input \"{:s}\"".format(
            option, arg_str, inp)

