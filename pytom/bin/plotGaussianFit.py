#!/usr/bin/env pytom

"""
Created on Jun 27, 2011
Updated on Jun 21, 2021

@author: yuxiangchen, Marten Chaillet
"""
import sys
import os


if __name__ == '__main__':
    # parse command line arguments with ScriptHelper2

    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.plotting.plistQuality import parse_linker_file_w_gt, plist_quality_gaussian_fit, \
        plist_quality_ground_truth, read_ground_truth, select_ground_truth
    from pytom.basic.structures import ParticleList

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Fit a Gaussian to particle list peak.',
        authors='Marten Chaillet (based on Yuxiang Chen)',
        options=[
            ScriptOption2(['-f', '--file'], 'Particle list xml file', 'file', 'required'),
            ScriptOption2(['-o', '--output'], 'Filename of output or directory for output.', 'string',
                          'optional', '.'),
            ScriptOption2(['-g', '--ground_truth'], 'File with simulated ground truth data.', 'file', 'optional'),
            ScriptOption2(['--gt_template_id'], 'Template id to select from ground truth file in case providing a '
                                                'single particle list with no linker information.',
                          'string', 'optional'),
            ScriptOption2(['--gt_zoffset'], 'Z offset for ground truth locations.', 'int',
                          'optional', 0),
            ScriptOption2(['--gt_tolerance'], 'Distance tolerance for correct match between estimated particle '
                                              'position and ground truth position', 'int', 'optional', 5),
            ScriptOption2(['-n', '--numberBins'], 'Numbers of bins to plot the score distribution', 'int',
                          'optional', 30),
            ScriptOption2(['-p', '--gaussianPeak'], 'Index of the peak (p < nbins) to be taken as a starting point '
                                                    'for searching the mean of the particle population', 'int',
                          'optional'),
            ScriptOption2(['--forcePeak'], 'Force the given peak index to be the mean', 'no arguments', 'optional'),
            ScriptOption2(['--cropPlot'], 'Crops distribution on y-axis to not show full height, '
                                          'can make particle peak more clearly visible', 'no arguments', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    input_file, output, gt_file, gt_template_id, gt_zoffset, gt_distance_tolerance, number_of_bins, peak_index, \
        force_peak, crop_plot = options

    # parse peak index
    if ((peak_index is None) or (peak_index > number_of_bins)) and gt_file is None:
        print(' - fall back to default peak index estimate, nbins / 2')
        peak_index = number_of_bins // 2
        if force_peak:
            print(' - WARNING: peak is forced without specifying a peak index, the default peak estimate at bin number '
                  'nbins/2 will be enforced.')

    # select function to call based on input extension
    input_name = os.path.splitext(os.path.split(input_file)[1])[0]
    input_ext = os.path.splitext(input_file)[1]
    if input_ext == '.txt':
        """
        Linker file should be of the following format, every line contains a template id (corresponding to ground 
        truth file) and a particle list file that can be loaded from the current directory :
        
    particle_id particle_list
    1BXN 1BXN_20px_10A_3.5um_cpu_plist.xml
    1QVR 1QVR_30px_10A_3.5um_cpu_plist.xml
    1S3X 1S3X_12px_10A_3.5um_cpu_plist.xml
    ...
    1U6G 1U6G_22px_10A_3.5um_cpu_plist.xml
        """
        if os.path.exists(output) and os.path.isdir(output):
            parse_linker_file_w_gt(input_file, input_name, gt_file, gt_zoffset, output, number_of_bins)
        else:
            print('Invalid output directory.')
            sys.exit(0)
    elif input_ext == '.xml':
        # load particle list
        particle_list = ParticleList()
        particle_list.fromXMLFile(input_file)
        output_figure_name = output if output != '.' else None
        # evaluate with ground truth or not, depending on if its provided
        if gt_file is not None:
            ground_truth_names, ground_truth_locations, _ = read_ground_truth(gt_file, gt_zoffset)
            ground_truth = select_ground_truth(ground_truth_names, ground_truth_locations, gt_template_id)
            if len(ground_truth) != 0:
                plist_quality_ground_truth(particle_list, ground_truth, output_figure_name=output_figure_name,
                                           nbins=number_of_bins, position_tolerance=gt_distance_tolerance)
            else:
                print('No locations of this particle id in the ground truth file. Exiting.')
                sys.exit(0)
        else:
            plist_quality_gaussian_fit(particle_list, peak_index, force_peak=force_peak, crop_hist=crop_plot,
                                       num_bins=number_of_bins, output_figure_name=output_figure_name)
    else:
        print('Invalid input. Provide either a linker .txt file or an particle list in .xml format.')
        sys.exit(0)
