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


class ScriptOption2:
    """
    ScriptOption: Determine whether script option requires an argument or is optional
    """

    def __init__(self, option_str, description, arguments, required, default_value=None):
        """
        @param option_str: a list of names for the option ex. ['--particleList', '-pl']
        @type option_str: list(str)
        @param description: a description for the option, will be shown when the user types --help, some details will
            automatically be added (required or not and arguments or not and the format of the arguments)
        @param arguments: requires argument? choose between 'has arguments' and 'no arguments' or give the format of the
            arguments, consisting of possibly splitting characters and one or more of the following: 'uint', 'int',
            'float' and 'string'. Uint will only accept positive integers, int positive and negatice integers, float
            will accept any valid float (including exponent) and string will match anything that is not whitespace (so
            be careful). No whitespace is permitted in these format strings.
        @type arguments: str
        @param required: optional feature? choose between 'required' and 'optional'
        @type required: str
        @param default_value: the default value for the option, to be used if it is not provided by the user
        @type default_value: any
        """
        if option_str.__class__ == list:  # make it list
            self.option_str = option_str
        else:
            self.option_str = [option_str]

        self.description = description

        if arguments.lower() == 'no arguments':
            self.arg = False
        elif arguments.lower() == 'has arguments':
            self.arg = True
        else:
            # Handle invalid strings, whitespace is invalid because of the algorithm used in parse_script_options
            import re
            match = re.match(r"\s", arguments)
            if match:
                raise Exception(
                    "Argument strings should not contain spaces (or other whitespace) or it should be \"no arguments\" or \"has arguments\". For option {:s} with argument string \"{:s}\"".format(
                        self.option_str[0], arguments))
            self.arg = True

        self.arg_str = arguments.lower()

        if required.lower() == 'required':
            self.required = True
        elif required.lower() == 'optional':
            self.required = False
        else:
            raise Exception("The value for 'required' in ScriptOption should be \"required\" or \"optional\"")

        self.default = default_value

    def __str__(self):
        """Create a string for display in the help menu"""
        name = ''
        for n in self.option_str:
            name += str(n) + ', '
        name = name[:-2]

        arguments_details = ('requires arguments' + ((
                                                                 ' in the format: \"' + self.arg_str + '\"') if self.arg_str != 'has arguments' else '')) if self.arg else 'has no arguments'
        default_value_details = '' if self.default is None else ' with default value: \"' + str(self.default) + '\"'
        details = "(Is " + (
            'required' if self.required else 'optional') + " and " + arguments_details + default_value_details + ")"

        return name + "\n    " + self.description + "\n    " + details


class ScriptHelper2:
    """
    ScriptHelper2: A constructor to help in parsing command line arguments
    No help option has to be provided as this is generated automatically
    """

    def __init__(self, name, **kwargs):
        self.name = name
        self.description = kwargs.get('description', None)
        args = kwargs.get('args', [])
        if args.__class__ == list:  # make it list
            self.args = args
        else:
            self.args = [args]
        options = kwargs.get('options', [])
        if options.__class__ == list:  # make it list
            self.options = options
        else:
            self.options = [options]
        self.options.append(ScriptOption2(['--help'], 'Print this help.', 'no arguments', 'optional'))
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
                help_str += str(opt) + "\n"
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