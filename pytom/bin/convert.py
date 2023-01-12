#!/usr/bin/env pytom

"""
Created on October, 2019

@author: Douwe Schulte

To join the conversion scripts scattered all around the bin folder into one single script.
"""

import pytom.basic.files as f
import os

# The lookup tables for the conversion functions
extensions = all_extensions = ["em", "mrc", "ccp4", "pl", "meta", "pdb", "mmCIF", 'mrcs', 'st', 'mdoc', 'h5', 'star', 'txt', 'wimp', 'xml', 'xf']
special = ['pl']

lookuptable = []
for n in range(len(extensions)):
    lookuptable.append([])
    [lookuptable[n].append('') for m in range(len(extensions))]

for n, e1 in enumerate(extensions):
    for m, e2 in enumerate(extensions):
        if e1 == e2:
            lookuptable[n][m] = 'Same'
            continue
        func = f'{e1}2{e2}'
        if hasattr(f,func):
            lookuptable[n][m] = getattr(f,func)
        else:
            lookuptable[n][m] = None

# lookuptable = [
#     ["Same",        f.em2mrc,  f.em2ccp4],
#     [f.mrc2em,      "Same",    f.mrc2ccp4],
#     [f.ccp42em,     f.ccp42em, "Same"],
# ]
special_extensions = ["pdb", "mmCIF"]
special_lookuptable = [f.pdb2em, f.mmCIF2em, f.mrcs2mrc, f.mrcs2mrc, f.mdoc2meta]


def special(file, format, target, chaindata):
    """Converts a special fileformat to the given format"""
    in_ex = file.split('.')[-1]

    func = special_lookuptable[special_extensions.index(in_ex)]

    volume = func(file, chaindata[0], chaindata[1], invertDensity=chaindata[2], chain=chaindata[3])

    # Write it as an em file
    em_name  = f.name_to_format(file, target, 'em')
    f.write_em(em_name, volume)

    if format != "em":
        # Convert em file to target
        convertfile(em_name, format, target, chaindata)
        os.remove(em_name)


def convertfile(file, format, target, chaindata, subtomo_prefix=None, wedge_angles=None):
    """Convert a file to the given format, returns (statuscode, errormessage)"""
    in_ex = file.split('.')[-1]
    if in_ex in ("map", 'rec', 'st'): in_ex = "mrc"

    if in_ex == 'xml': in_ex = 'pl'

    # Check if func exists
    try:
        index_from = extensions.index(in_ex)
    except:
        return -1, "The given input file format is not supported: {:s}".format(in_ex)

    index_to = extensions.index(format)

    func = lookuptable[index_from][index_to]


    if in_ex in special_extensions and func is None:
        if chaindata is None:
            return 1, "For conversion from special format ({:s} and {:s}) chaindata is needed".format(", ".join([a for a in special_extensions[:-1]]), special_extensions[-1])
        special(file, format, target, chaindata)
        return 0, ""

    if in_ex == "txt" and format == "pl":
        try:
            f.txt2pl(f.name_to_format(file, target, format), file, subtomoPrefix=subtomo_prefix, wedgeAngle=wedge_angles)
            return 0,""

        except:
            return 1, "Exception while converting coordinatesfile {:s} to a particlelist".format(file)


    if (in_ex == "pl" and format == "star") or (in_ex == "star" and format == "xml") or (in_ex == "log" and format == "txt"):
        try:
            func(file, target, **chaindata)
            return 0,""

        except Exception as e:
            print(e)
            return 1, "Exception while converting {:s} file {:s} to a star file".format(in_ex, file)


    if func is None:
        return -1, "Converting from {:s} to {:s} is not implemented yet".format(in_ex, format)
    elif func == "Same":
        return -1, "Converting from {:s} to the same format seems not needed".format(in_ex)
    else:
        func(file, target, subtomo_prefix)
        return 0, ""


def print_errors(message):
    """To print to stderr and use fancy styling (if possible), for multiple errors"""
    import sys
    color1 = ""
    color2 = ""

    if sys.stdout.isatty():
        color1 = "\033[91m"
        color2 = "\033[0m"

    sys.stderr.write(color1 + message + color2 + "\nSee --help for help on how to use this function\n")
    sys.exit()


def print_single_error(message, warning=False):
    """
    To print to stderr and use fancy styling (if possible), on a single line

    :param message: The error message to show
    :type message: str
    :param warning: If the error is an error or a warning, changes the style
    :type warning: bool
    :return: void
    :returntype: void
    """
    import sys
    color = "93" if warning else "91"

    color1 = "WARNING: " if warning else "ERROR: "
    color2 = ""

    if sys.stdout.isatty():
        color1 = "\033[" + color + "m"
        color2 = "\033[0m"

    sys.stderr.write("{:s}{:s}{:s}\n".format(color1, message, color2))


def warn_if_file_exists(filename):
    if os.path.exists(filename):
        print_single_error("File: {:s} already exists will be overwritten".format(filename), warning=True)


def test_validity(filename, directory, target, format, chaindata):
    """Test the validity of the input arguments"""
    exceptions = ""

    # Determine the validity of the call
    if chaindata and chaindata[2] > 1:
        exceptions += "The only valid values for invertDensity are 0 and 1\n"

    if filename and directory:
        exceptions += "Only use one input at a time, only --file or --directory\n"

    if not format in all_extensions:
        exceptions += "The output format is not valid, choose between: {:s} and {:s}\n".format(
            ", ".join([a for a in all_extensions[:-1]]), all_extensions[-1])

    if not filename and not directory:
        exceptions += "There needs to be at least one input, so provide --file or --directory\n"

    if exceptions != "":
        print_errors(exceptions)


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Converts file(s) to the given output file, converts freely between {:s} and {:s} and converts '
                    'from {:s} and {:s}'.format(", ".join([a for a in all_extensions[:-1]]), all_extensions[-1], ", ".join([a for a in special_extensions[:-1]]), special_extensions[-1]),
        authors='Douwe Schulte',
        options=[ScriptOption2(
                     ['-f', '--file'], 'Filename will convert this single file to the target format, not to be used in '
                                       'combination with --directory', 'file', 'optional'),
                 ScriptOption2(
                     ['-d', '--directory'], 'A directory of files, will convert all possible files to the outputFormat '
                                            'specified, not te be used in combination with --file', 'directory', 'optional'),
                 ScriptOption2(
                     ['-t', '--targetPath'], 'Path to new file, as a relative or absolute path to a directory, the file'
                                             'name will be the same as the original file', 'directory', 'optional', '.'),
                 ScriptOption2(
                     ['-o', '--outputFormat'], 'The format of the output file, will keep the same name but only change '
                                               'the extension.', 'string', 'required'),

                 ScriptOption2(['--outname'], 'Filename output file. This name has presedence over automatic output filename determination.', 'string', 'optional', ''),

                 ScriptOption2(
                     ['-c', '--chainData'], 'Data needed for conversion from the special formats ({:s} and {:s}), in the format \
                     pixelsSize,cubeSize,invertDensity,chain (invertDensity is 1 or 0)'.format(", ".join([a for a in special_extensions[:-1]]), special_extensions[-1]),
                     'float,uint,uint,string', 'optional'),

                 ScriptOption2(
                     ['-s', '--prefix'], 'Data needed for the conversion from coordinates to a particlelist. '
                                         'Path and filename for subtomogram files (e.g., MyPath/particle_)',
                     'string', 'optional'),
                 ScriptOption2(
                     ['-w', '--wedgeAngles'], 'Data needed for the conversion from coordinates to a particlelist. '
                                              'Missing wedge angle(s) [counter-clock, clock] or single angle',
                     'has arguments', 'optional'),



                 ScriptOption2(['--prefixQuery'], 'Only files that match the prefix will be converted.'
                              , 'string', 'optional', ''),
                 ScriptOption2(['--suffixQuery'], 'Only files that end with suffix will be converted.'
                          , 'string', 'optional', ''),


                 ScriptOption2(['--pixelSize'], 'Pixelsize in Angstrom of original nanographs.', 'float', 'optional', 1),
                 ScriptOption2(['--binPyTom'], '(Linear) Binning factor of the pytom tomogram.', 'float', 'optional', 1),
                 ScriptOption2(['--binWarpM'], '(Linear) Binning factor of the warp/m volumes.', 'float', 'optional', 1),
                 ScriptOption2(['--tlt-file'], 'Either .tlt file or .rawtlt in case tilt angles have not been '
                                               'optimized during alignment. ', 'file', 'optional', ''),
                 ScriptOption2(['--sortedFolder'], 'Sorted Images are located in this folder (either sorted or sorted_ctf)', 'directory', 'optional', ''),
        ])

    #TODO write --filter to filter input files maybe on (glob) pattern or else on extension or similar

    filename, directory, target, format, outname, chaindata, subtomo_prefix, w, prefix, suffix, pixelsize, binningFactorPyTom,\
    binningFactorWarpM, tlt_file, sorted_folder = parse_script_options2(sys.argv[1:], helper)

    try:
        if w:
            split = w.split(',')
            if len(split) == 2:
                wedge_angles = []
                wedge_angles.append(float(split[0]))
                wedge_angles.append(float(split[1]))
            elif len(split) == 1:
                wedge_angles = [float(w),]*2
            else:
                raise Exception("The amount of angles given to wedge angle is not valid (has to be one or two).")
        else:
            wedge_angles = None
    except Exception as e:
        print_errors("The parsing of the wedge angle was not successful.\n{:s}".format(str(e)))




    # Parse the pattern

    # if pattern:
    #     # Escape the filter so that it can be compiled without issues
    #     pattern = re.escape(pattern)
    #     # Create the pattern matches
    #     pattern = pattern.replace('?', '??').replace('$', '.*?')
    #
    #     if len(pattern.split('/')) > 1:
    #         flags = pattern.split('/')[-1]
    #         pattern = '/'.join(pattern.split('/')[:-1])
    #     else:
    #         flags = ''
    #         pattern = pattern
    #
    #     if 'f' in flags:
    #         pattern = '^' + pattern + '$'
    #
    #     try:
    #         if 'i' in flags:
    #             input_filter = re.compile(pattern, re.IGNORECASE)
    #         else:
    #             input_filter = re.compile(pattern)
    #     except:
    #         print_errors("The parsing of the filter pattern was not successful.")
    # else:
    #     input_filter = None

    # Test for validity of the arguments passed, will stop execution if an error is found
    test_validity(filename, directory, target, format, chaindata)

    # TODO Rename this dictionary to something else, chaindata is an input argument for pdbs/CIF
    if chaindata is None:
        chaindata = {}
        chaindata['pixelsize'] = pixelsize
        chaindata['binningWarpM'] = binningFactorWarpM
        chaindata['binningPyTom'] = binningFactorPyTom
        outname = f'dummy.{format}' if (outname == '' and directory) else outname
        chaindata['outname'] = outname
        chaindata['tlt_file'] = tlt_file
        chaindata['sorted_folder'] = sorted_folder
        chaindata['wedgeAngles'] = wedge_angles

    if filename:
        warn_if_file_exists(f.name_to_format(filename, target, format) if outname == ''
                            else f.name_to_format(outname, target, format))

        num, ex = convertfile(filename, format, target, chaindata, subtomo_prefix, wedge_angles)

        # Print possible errors
        if num != 0:
            print_errors(ex)

    elif directory:
        import glob

        if prefix == '' and suffix == '':
            fileList = [os.path.join(directory, fname) for fname in os.listdir(directory) if not os.path.isdir(fname)]
        else:
            query = os.path.join(directory, prefix + '*' + suffix)
            fileList = sorted(glob.glob(query))


        if format == 'star':
            lookuptable[extensions.index('pl')][extensions.index('star')](directory, target, pixelsize=pixelsize,
                                                                          binningPyTom=binningFactorPyTom,
                                                                          binningWarpM=binningFactorWarpM,
                                                                          outname=outname)
        elif format == 'xml':
            lookuptable[extensions.index('star')][extensions.index('pl')](directory, target, pixelsize=pixelsize,
                                                                          binningPyTom=binningFactorPyTom,
                                                                          binningWarpM=binningFactorWarpM,
                                                                          outname=outname)
        else:
            for filename in fileList:
                warn_if_file_exists(f.name_to_format(filename, target, format))

                try:
                    num, ex = convertfile(filename, format, target, chaindata, subtomo_prefix, wedge_angles)
                except Exception as e:
                    num = 1
                    ex = "CONVERSION EXCEPTION " + str(e)

                # Print the result, ignores status code -1
                if num == 0:
                    print("Converted {:s} to {:s}".format(filename, f.name_to_format(filename, target, format)))
                elif num == 1:
                    print_single_error("File: {:s} gave error: {:s}".format(filename, ex))