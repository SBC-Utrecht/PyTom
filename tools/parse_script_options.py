'''
Created on Jun 29, 2011

@author: yuxiangchen
'''

import sys

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
