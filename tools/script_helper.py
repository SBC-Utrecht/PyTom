'''
Created on Jun 29, 2011

@author: yuxiangchen
'''

class ScriptArg:
    def __init__(self, arg_str, description=""):
        self.arg_str = arg_str
        self.description = description
    
    def __str__(self):
        return self.arg_str + "    " + self.description

class ScriptOption:
    """
    ScriptOption: Determine whether script option requires an argument or is optional 
    """
    def __init__(self, option_str, description="", arg=True, optional=False):
        """
	@param arg: requires argument?
	@type arg: bool
	@param optional: optional feature?
	@type optional: bool
        """
        if option_str.__class__ == list: # make it list
            self.option_str = option_str
        else:
            self.option_str = [option_str]
        self.description = description
        self.arg = arg
        self.optional = optional
    
    def __str__(self):
        name = ''
        for n in self.option_str:
            name += str(n)+', '
        name = name[:-2]

        return name + "    " + self.description + " (Is optional: " + ('Yes' if self.optional else 'No') + ";" + " Requires arguments: "+ ('Yes' if self.arg else 'No') + ")"

class ScriptHelper:
    """
    ScriptHelper: 
    """
    def __init__(self, name, **kwargs):
        self.name = name
        self.description = kwargs.get('description', None)
        args = kwargs.get('args', [])
        if args.__class__ == list: # make it list
            self.args = args
        else:
            self.args = [args]
        options = kwargs.get('options', [])
        if options.__class__ == list: # make it list
            self.options = options
        else:
            self.options = [options]
        self.examples = kwargs.get('examples', None)
        self.see_also = kwargs.get('see_also', None)
        self.version = kwargs.get('version', None)
        self.authors = kwargs.get('authors', None)
        self.license = kwargs.get('license', None)
    
    def __str__(self):
        help_str = ""
        help_str += "NAME\n    %s\n" % self.name
        if self.description is not None:
            help_str += "DESCRIPTION\n    %s\n" % self.description
        if len(self.args) > 0:
            help_str += "ARGUMENTS\n"
            for arg in self.args:
                help_str += "    " + str(arg) + "\n"
        if len(self.options) > 0:
            help_str += "OPTIONS\n"
            for opt in self.options:
                help_str += "    " + str(opt) + "\n"
        if self.examples is not None:
            help_str += "EXAMPLES\n    %s\n" % self.examples
        if self.see_also is not None:
            help_str += "SEE ALSO\n    %s\n" % self.see_also
        if self.version is not None:
            help_str += "VERSION\n    %s\n" % self.version
        if self.authors is not None:
            help_str += "AUTHORS\n    %s\n" % self.authors
        if self.license is not None:
            help_str += "LICENSE\n    %s\n" % self.license
                
        return help_str
