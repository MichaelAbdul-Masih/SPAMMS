
__version__ = '1.2.0'

import numpy as np
import glob
import os
import shutil
import sys
import itertools
import time
import functools
from scipy.interpolate import splrep, splev
from scipy import stats
from scipy.special import erf
from scipy.signal import fftconvolve
import scipy.optimize as so
import phoebe
from phoebe import u,c
import math
from tqdm import tqdm, trange
import settings
from astropy.constants import R_sun, M_sun, G
import getopt
# import tarfile
# from io import BytesIO




def read_input_file(input_file):
    print('Reading input file...')
    lines = tuple(open(input_file, 'r'))
    object_type_ind = [i for i in range(len(lines)) if lines[i].startswith('object_type')][0]
    object_type = lines[object_type_ind].split('=')[1].strip()
    if object_type == 'single':
        return read_s_input_file(input_file)
    elif object_type == 'contact_binary':
        return read_cb_input_file(input_file)
    elif object_type == 'binary':
        return read_b_input_file(input_file)
    else:
        raise ValueError('object_type must be one of the following [single, binary, contact_binary]')


def read_cb_input_file(input_file):
    lines = tuple(open(input_file, 'r'))

    object_type_ind = [i for i in range(len(lines)) if lines[i].startswith('object_type')][0]
    object_type = lines[object_type_ind].split('=')[1].strip()

    ntriangles = 5000
    ntriangles_ind = [i for i in range(len(lines)) if lines[i].startswith('ntriangles')]
    if len(ntriangles_ind) >= 1:
        ntriangles = lines[ntriangles_ind[0]].split('=')[1].strip()

    path_to_obs_spectra_ind = [i for i in range(len(lines)) if lines[i].startswith('path_to_obs_spectra')][0]
    path_to_obs_spectra = lines[path_to_obs_spectra_ind].split('=')[1].strip()
    if not path_to_obs_spectra.endswith('/'): path_to_obs_spectra += '/'
    if path_to_obs_spectra == 'None/': path_to_obs_spectra = 'None'

    output_directory_ind = [i for i in range(len(lines)) if lines[i].startswith('output_directory')][0]
    output_directory = lines[output_directory_ind].split('=')[1].strip()
    if not output_directory.endswith('/'): output_directory += '/'

    path_to_grid_ind = [i for i in range(len(lines)) if lines[i].startswith('path_to_grid')][0]
    path_to_grid = lines[path_to_grid_ind].split('=')[1].strip()
    if not path_to_grid.endswith('/'): path_to_grid += '/'

    try:
        grid_type_ind = [i for i in range(len(lines)) if lines[i].startswith('grid_type')][0]
        grid_type = lines[grid_type_ind].split('=')[1].strip()
        if grid_type not in ['FW', 'FWNN', 'T', 'K']:
            grid_type = 'FW'
    except:
        grid_type = 'FW'


    fit_params = ['fillout_factor', 'teff_primary', 'teff_secondary', 'period', 'sma', 'inclination', 'q', 't0', 'async_primary', 'async_secondary', 'gamma', 'v_macro', 'v_micro', 'metallicity', 'alpha_enhancement']
    abundance_params = ['he_abundances', 'cno_abundances']

    fit_param_values = {'async_primary':1.0, 'async_secondary':1.0, 'v_macro':0.0, 'v_micro':10.0, 'metallicity':1.000, 'alpha_enhancement':0.0}
    if grid_type == 'K':
        fit_param_values['metallicity'] = 0.000
    abund_param_values = {}
    io_dict = {'object_type':object_type, 'ntriangles':ntriangles, 'path_to_obs_spectra':path_to_obs_spectra, 'output_directory':output_directory, 'grid_type':grid_type, 'path_to_grid':path_to_grid, 'input_file':input_file, 'rad_bound':False}
    try:
        times_ind = [i for i in range(len(lines)) if lines[i].startswith('times')][0]
        times = lines[times_ind].split('=')[1].strip()
        io_dict['times'] = arg_parse(times)
    except:
        pass

    for param in fit_params:
        if param in fit_param_values.keys():
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except:
                fit_param_values[param] = arg_parse(str(fit_param_values[param]))
        else:
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except IndexError:
                # Handle the IndexError here
                print(f"Error: {param} not found in input file")

    if grid_type == 'K':
        alpha = np.where(np.array(fit_param_values['alpha_enhancement']) >= 1, 1, 0)
        fit_param_values['alpha_enhancement'] = list(set(alpha))

    abund_param_values['he_abundances'] = [0.06, 0.1, 0.15, 0.2]
    abund_param_values['cno_abundances'] = [6.5, 7.0, 7.5, 8.0, 8.5]
    abund_param_values['lp_bins'] = 161
    for param in abundance_params:
        try:
            arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
            abund_param_values[param] = arg_parse(arg)
        except:
            abund_param_values[param] = arg_parse(str(abund_param_values[param]))

    if type(abund_param_values['lp_bins']) is list:
        abund_param_values['lp_bins'] = int(abund_param_values['lp_bins'][0])

    abund_param_values['interpolate_abundances'] = False

    # interp_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('interpolate_abundances')][0]].split('=')[1].strip()
    # if interp_arg == 'False' or interp_arg == '0':
    #     abund_param_values['interpolate_abundances'] = False
    # else:
    #     abund_param_values['interpolate_abundances'] = True

    if grid_type in ['FWNN', 'FW']:
        line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_line_list =')][0]].split('=')[1].strip()
        line_list = parse_line_list(line_list_arg)
    elif grid_type in ['T', 'K']:
        line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_wavelength_range =')][0]].split('=')[1].strip()
        line_list = parse_wavelength_range(line_list_arg)


    return fit_param_values, abund_param_values, line_list, io_dict


def read_b_input_file(input_file):
    lines = tuple(open(input_file, 'r'))

    object_type_ind = [i for i in range(len(lines)) if lines[i].startswith('object_type')][0]
    object_type = lines[object_type_ind].split('=')[1].strip()

    ntriangles = 5000
    ntriangles_ind = [i for i in range(len(lines)) if lines[i].startswith('ntriangles')]
    if len(ntriangles_ind) >= 1:
        ntriangles = lines[ntriangles_ind[0]].split('=')[1].strip()

    path_to_obs_spectra_ind = [i for i in range(len(lines)) if lines[i].startswith('path_to_obs_spectra')][0]
    path_to_obs_spectra = lines[path_to_obs_spectra_ind].split('=')[1].strip()
    if not path_to_obs_spectra.endswith('/'): path_to_obs_spectra += '/'
    if path_to_obs_spectra == 'None/': path_to_obs_spectra = 'None'

    output_directory_ind = [i for i in range(len(lines)) if lines[i].startswith('output_directory')][0]
    output_directory = lines[output_directory_ind].split('=')[1].strip()
    if not output_directory.endswith('/'): output_directory += '/'

    path_to_grid_ind = [i for i in range(len(lines)) if lines[i].startswith('path_to_grid')][0]
    path_to_grid = lines[path_to_grid_ind].split('=')[1].strip()
    if not path_to_grid.endswith('/'): path_to_grid += '/'

    try:
        grid_type_ind = [i for i in range(len(lines)) if lines[i].startswith('grid_type')][0]
        grid_type = lines[grid_type_ind].split('=')[1].strip()
        if grid_type not in ['FW', 'FWNN', 'T', 'K']:
            grid_type = 'FW'
    except:
        grid_type = 'FW'

    try:
        dist_ind = [i for i in range(len(lines)) if lines[i].startswith('distortion')][0]
        dist = lines[dist_ind].split('=')[1].strip()
        if dist not in ['rotstar', 'roche', 'sphere']:
            dist = 'roche'
    except:
        dist = 'roche'

    fit_params = ['r_equiv_primary', 'r_equiv_secondary', 'teff_primary', 'teff_secondary', 'period', 'sma', 'inclination', 'q', 't0', 'async_primary', 'async_secondary', 'pitch_primary', 'pitch_secondary', 'yaw_primary', 'yaw_secondary', 'gamma', 'v_macro', 'v_micro', 'metallicity', 'alpha_enhancement']
    abundance_params = ['he_abundances', 'cno_abundances']

    fit_param_values = {'async_primary':1.0, 'async_secondary':1.0, 'pitch_primary':0.0, 'pitch_secondary':0.0, 'yaw_primary':0.0, 'yaw_secondary':0.0, 'v_macro':0.0, 'v_micro':10.0, 'metallicity':1.0, 'alpha_enhancement':0.0}
    if grid_type == 'K':
        fit_param_values['metallicity'] = 0.00
    abund_param_values = {}
    io_dict = {'object_type':object_type, 'ntriangles':ntriangles, 'path_to_obs_spectra':path_to_obs_spectra, 'output_directory':output_directory, 'path_to_grid':path_to_grid, 'grid_type':grid_type, 'input_file':input_file, 'distortion':dist, 'rad_bound':False}
    try:
        times_ind = [i for i in range(len(lines)) if lines[i].startswith('times')][0]
        times = lines[times_ind].split('=')[1].strip()
        io_dict['times'] = arg_parse(times)
    except:
        pass

    for param in fit_params:
        if param in fit_param_values.keys():
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except:
                fit_param_values[param] = arg_parse(str(fit_param_values[param]))
        else:
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except IndexError:
                # Handle the IndexError here
                print(f"Error: {param} not found in input file")

    if grid_type == 'K':
        alpha = np.where(np.array(fit_param_values['alpha_enhancement']) >= 1, 1, 0)
        fit_param_values['alpha_enhancement'] = list(set(alpha))

    abund_param_values['he_abundances'] = [0.06, 0.1, 0.15, 0.2]
    abund_param_values['cno_abundances'] = [6.5, 7.0, 7.5, 8.0, 8.5]
    abund_param_values['lp_bins'] = 161
    for param in abundance_params:
        try:
            arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
            abund_param_values[param] = arg_parse(arg)
        except:
            abund_param_values[param] = arg_parse(str(abund_param_values[param]))

    if type(abund_param_values['lp_bins']) is list:
        abund_param_values['lp_bins'] = int(abund_param_values['lp_bins'][0])

    abund_param_values['interpolate_abundances'] = False
    # interp_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('interpolate_abundances')][0]].split('=')[1].strip()
    # if interp_arg == 'False' or interp_arg == '0':
    #     abund_param_values['interpolate_abundances'] = False
    # else:
    #     abund_param_values['interpolate_abundances'] = True

    if grid_type in ['FWNN', 'FW']:
        line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_line_list =')][0]].split('=')[1].strip()
        line_list = parse_line_list(line_list_arg)
    elif grid_type in ['T', 'K']:
        line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_wavelength_range =')][0]].split('=')[1].strip()
        line_list = parse_wavelength_range(line_list_arg)

    return fit_param_values, abund_param_values, line_list, io_dict


def read_s_input_file(input_file):
    lines = tuple(open(input_file, 'r'))

    object_type_ind = [i for i in range(len(lines)) if lines[i].startswith('object_type')][0]
    object_type = lines[object_type_ind].split('=')[1].strip()

    ntriangles = 5000
    ntriangles_ind = [i for i in range(len(lines)) if lines[i].startswith('ntriangles')]
    if len(ntriangles_ind) >= 1:
        ntriangles = lines[ntriangles_ind[0]].split('=')[1].strip()

    path_to_obs_spectra_ind = [i for i in range(len(lines)) if lines[i].startswith('path_to_obs_spectra')][0]
    path_to_obs_spectra = lines[path_to_obs_spectra_ind].split('=')[1].strip()
    if not path_to_obs_spectra.endswith('/'): path_to_obs_spectra += '/'
    if path_to_obs_spectra == 'None/': path_to_obs_spectra = 'None'

    output_directory_ind = [i for i in range(len(lines)) if lines[i].startswith('output_directory')][0]
    output_directory = lines[output_directory_ind].split('=')[1].strip()
    if not output_directory.endswith('/'): output_directory += '/'

    path_to_grid_ind = [i for i in range(len(lines)) if lines[i].startswith('path_to_grid')][0]
    path_to_grid = lines[path_to_grid_ind].split('=')[1].strip()
    if not path_to_grid.endswith('/'): path_to_grid += '/'

    try:
        grid_type_ind = [i for i in range(len(lines)) if lines[i].startswith('grid_type')][0]
        grid_type = lines[grid_type_ind].split('=')[1].strip()
        if grid_type not in ['FW', 'FWNN', 'T', 'K']:
            grid_type = 'FW'
    except:
        grid_type = 'FW'

    try:
        dist_ind = [i for i in range(len(lines)) if lines[i].startswith('distortion')][0]
        dist = lines[dist_ind].split('=')[1].strip()
        if dist not in ['rotstar', 'roche', 'sphere']:
            dist = 'rotstar'
    except:
        dist = 'rotstar'

    try:
        gd_ind = [i for i in range(len(lines)) if lines[i].startswith('gravity_darkening')][0]
        gd = lines[gd_ind].split('=')[1].strip()
        if gd not in ['VZ', 'EL']:
            gd = 'VZ'
    except:
        gd = 'VZ'


    fit_params = ['teff', 'rotation_rate', 'requiv', 'inclination', 'mass', 't0', 'gamma']
    fit_params_alt = ['teff', 'vsini', 'rotation_rate', 'v_crit_frac', 'requiv', 'r_pole', 'inclination', 'mass', 't0', 'gamma', 'v_macro', 'sigma_R', 'sigma_T',  'v_micro', 'metallicity', 'alpha_enhancement']
    abundance_params = ['he_abundances', 'cno_abundances']

    fit_param_values = {'v_macro':0.0, 'sigma_R':0.0, 'sigma_T':0.0, 'v_micro':10.0, 'metallicity':1.0, 'alpha_enhancement':0.0}
    if grid_type == 'K':
        fit_param_values['metallicity'] = 0.00
    abund_param_values = {}
    io_dict = {'object_type':object_type, 'ntriangles':ntriangles, 'path_to_obs_spectra':path_to_obs_spectra, 'output_directory':output_directory, 'path_to_grid':path_to_grid, 'grid_type':grid_type, 'input_file':input_file, 'distortion':dist, 'gravity_darkening':gd, 'rad_bound':False}
    try:
        times_ind = [i for i in range(len(lines)) if lines[i].startswith('times')][0]
        times = lines[times_ind].split('=')[1].strip()
        io_dict['times'] = arg_parse(times)
    except:
        pass

    for param in fit_params_alt:
        if param in fit_param_values.keys():
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except:
                fit_param_values[param] = arg_parse(str(fit_param_values[param]))
        elif param in ['vsini', 'rotation_rate', 'v_crit_frac', 'requiv', 'r_pole', 'inclination']:
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except:
                fit_param_values[param] = [-1.0]
        else:
            try:
                arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
                fit_param_values[param] = arg_parse(arg)
            except IndexError:
                # Handle the IndexError here
                print(f"Error: {param} not found in input file")

    if (max(fit_param_values['sigma_T']) > 0 or max(fit_param_values['sigma_R']) > 0) and np.any(np.array(fit_param_values['v_macro']) > 0):
        print('vmacro and sigma_R/sigma_T cannot both be greater than 0. Defaulting to sigma_R/sigma_T and turning vmacro off')
        fit_param_values['v_macro'] = [-1.0]
    elif (max(fit_param_values['sigma_T']) > 0 or max(fit_param_values['sigma_R']) > 0):
        fit_param_values['v_macro'] = [-1.0]

    if grid_type == 'K':
        alpha = np.where(np.array(fit_param_values['alpha_enhancement']) >= 1, 1, 0)
        fit_param_values['alpha_enhancement'] = list(set(alpha))
    else:
        fit_param_values.pop('alpha_enhancement', None)
    
    if grid_type not in ['K', 'T']:
        fit_param_values.pop('metallicity', None)
        fit_param_values.pop('v_micro', None)

    abund_param_values['he_abundances'] = [0.06, 0.1, 0.15, 0.2]
    abund_param_values['cno_abundances'] = [6.5, 7.0, 7.5, 8.0, 8.5]
    abund_param_values['lp_bins'] = 161
    for param in abundance_params:
        try:
            arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
            abund_param_values[param] = arg_parse(arg)
        except:
            abund_param_values[param] = arg_parse(str(abund_param_values[param]))

    if type(abund_param_values['lp_bins']) is list:
        abund_param_values['lp_bins'] = int(abund_param_values['lp_bins'][0])

    abund_param_values['interpolate_abundances'] = False
    # interp_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('interpolate_abundances')][0]].split('=')[1].strip()
    # if interp_arg == 'False' or interp_arg == '0':
    #     abund_param_values['interpolate_abundances'] = False
    # else:
    #     abund_param_values['interpolate_abundances'] = True

    if grid_type in ['FWNN', 'FW']:
        line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_line_list =')][0]].split('=')[1].strip()
        line_list = parse_line_list(line_list_arg)
    elif grid_type in ['T', 'K']:
        line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_wavelength_range =')][0]].split('=')[1].strip()
        line_list = parse_wavelength_range(line_list_arg)

    return fit_param_values, abund_param_values, line_list, io_dict


def arg_parse(arg):
    if arg.startswith('('):
        dstr = [float(i) for i in arg.strip('()').split(',')]
        value = np.linspace(dstr[0], dstr[1], int(dstr[2]))
    elif arg.startswith('['):
        value = [float(i) for i in arg.strip('[]').split(',')]
    else:
        value = [float(arg)]
    return list(value)


def parse_line_list(arg):
    arg = arg.strip('[]').split(',')
    line_list = [i.strip().strip("'") for i in arg]
    return line_list


def parse_wavelength_range(arg):
    arg = arg.strip('[]').split(',')
    wave_arg_list = [i.strip().strip("'").split('-') for i in arg]
    return np.array(wave_arg_list, dtype='float')


def setup_output_directory(io_dict):
    '''
    Sets up output directory
    '''
    print('Setting up output directory...')
    output_directory = io_dict['output_directory']
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
    os.mkdir(output_directory)
    shutil.copy(io_dict['input_file'], output_directory + 'input.txt')
    try:
        shutil.copytree(io_dict['path_to_obs_spectra'], output_directory + 'input_spectra')
    except:
        pass

    print('Output Directory:  %s' %output_directory)


def check_input_spectra(io_dict):
    spec_files = glob.glob(io_dict['path_to_obs_spectra'] + '*_spec.txt')
    hjd_files = glob.glob(io_dict['path_to_obs_spectra'] + '*_hjd.txt')

    for spec_file in spec_files:
        spec = np.loadtxt(spec_file).T
        n_spec = len(spec[1:])

        expected_hjd_file = spec_file[:-9] + '_hjd.txt'

        try:
            hjd = np.loadtxt(expected_hjd_file, ndmin=1)
            if len(hjd) != n_spec:
                raise ValueError('Number of times in HJD file does not match number of Spectra in file: %s' %spec_file)

        except:
            raise IOError('Expected HJD file not found:  %s' %expected_hjd_file)

    spec_files_corenames = [i[:-9] for i in spec_files]
    hjd_files_corenames = [i[:-8] for i in hjd_files]
    dif_set = list(set(spec_files_corenames).symmetric_difference(set(hjd_files_corenames)))
    if len(dif_set) !=0:
        raise ValueError('Mismatch between Spectra and HJD core filenames %s' %dif_set)

    print('Input spectra checks complete')


def check_grid(times, abund_param_values, io_dict, grid_entries, run_dictionary):

    if io_dict['object_type'] == 'contact_binary':
        cb = run_cb_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
    elif io_dict['object_type'] == 'binary':
        cb = run_b_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
    else:
        if io_dict['distortion'] in ['rotstar', 'sphere']:
            cb = run_s_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
        else:
            cb = run_sb_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
    combs, mode_combs = determine_tgr_combinations(cb, io_dict, run_dictionary)

    missing_combs = [i for i in combs if i not in grid_entries]
    return missing_combs


def get_obs_spec_and_times(io_dict):
    if io_dict['path_to_obs_spectra'] == 'None':
        times = io_dict['times']
        obs_specs = None
        return np.array(times), obs_specs
    else:
        spec_files = glob.glob(io_dict['path_to_obs_spectra'] + '*_spec.txt')
        hjd_files = glob.glob(io_dict['path_to_obs_spectra'] + '*_hjd.txt')
        spec_files.sort()
        hjd_files.sort()
        times = []
        for hjd_file in hjd_files:
            times.extend(np.loadtxt(hjd_file, ndmin=1))
        obs_specs = {}
        for spec_file in spec_files:
            x = np.loadtxt(spec_file).T
            w = x[0]
            f = x[1:]
            times_temp = np.loadtxt(spec_file[:-8] + 'hjd.txt', ndmin=1)
            for i in range(len(f)):
                dic = {'wavelength':w, 'flux':f[i]}
                obs_specs[str(round(times_temp[i], 13)).ljust(13, '0')] = dic
        return np.array(times), obs_specs



def create_runs_and_ids(fit_param_values):
    keys = []
    values = []
    for key, value in fit_param_values.items():
        keys.append(key)
        values.append(value)

    runs = list(itertools.product(*values))
    run_dictionaries = [dict(zip(keys,i)) for i in runs]
    for i, j in enumerate(run_dictionaries):
        j['run_id'] = i

    run_ids = range(len(run_dictionaries))

    dictionary_of_run_dicts = dict(zip(run_ids, run_dictionaries))

    return run_dictionaries


def rpole_to_requiv(r_pole, vrot, n=5000, return_r_equator=False):
    '''
    r_pole - is the polar radius in units of solar radius
    vrot   - is the rotational velocity as a percentage of the critical velocity (value between 0 and 1)
    n      - is the number of theta points that we will use to plot the surface
    '''
    r_pole*=1.0
    vrot*=1.0
    n+=1

    # create an array of angles ranging from theta = 0 (upper pole) to theta = pi (lower pole)
    theta = np.linspace(0,np.pi, n)

    # to solve Eq. 6, we need to find the roots of the 3rd order polynomial. a, b, c and d are the coefficients of the polynomial
    a = (vrot / r_pole - vrot**3/(3.*r_pole))**2 * np.sin(theta)**2
    b = 0
    c = -3.0
    d = 3.*r_pole
    # we'll collect our radii at each theta in the rs array
    rs = []

    # for each theta we find the roots and check to make sure it is positive and not complex and dump them into rs
    for i in range(len(a)):
        x = np.roots([a[i], b, c, d])
        inds = np.iscomplex(x) == False
        y = x.real[inds]
        ind = np.argmin((r_pole*1.25 - y)**2)
        if y[ind] > 0:
            rs.append(y[ind])
        else:
            rs.append(y[ind]/2.*-1.0)

    # convert rs and thetas to xs and ys
    rs = np.array(rs)
    xs = rs * np.sin(theta)
    ys = rs * np.cos(theta)

    # now we can calculate the volume of our system using the disk integration method
    # we'll create a new equally spaced y array and a corresponding x array
    new_ys = np.linspace(-1*r_pole, r_pole, n)
    new_xs = np.interp(new_ys, ys[::-1], xs[::-1])

    # calculate the volume
    vol = 0
    dy = (2*r_pole) / (n-1)
    for i in range(1, n):
        r = (new_xs[i] + new_xs[i-1])/2
        vol += np.pi * r**2 *dy

    # calculate r_equiv from the calculated total volume
    r_equiv = (vol * 3./(4.*np.pi))**(1./3)

    if return_r_equator:
        return r_equiv, max(rs)
    else:
        return r_equiv


def func_requiv_to_rpole(rpole, vrot, requiv):
    return np.sqrt((rpole_to_requiv(rpole, vrot) - requiv)**2)


def requiv_to_rpole(requiv, vrot):
    '''
    requiv - is the volumetric equivalent spherical radius in units of solar radius
    vrot   - is the rotational velocity as a percentage of the critical velocity (value between 0 and 1)
    '''
    res = so.minimize(func_requiv_to_rpole, np.array([requiv]), args = (vrot,requiv), bounds=[(requiv*0.5,requiv*1.1)])
    return res.x[0]


def func_requiv_to_rpole_abs_units(rpole, vrot, requiv, mass):
    v_crit = calc_critical_velocity(mass, rpole)
    v_percent_crit = vrot/v_crit
    return np.sqrt((rpole_to_requiv(rpole, v_percent_crit) - requiv)**2)


def requiv_to_rpole_abs_units(requiv, vrot, mass):
    '''
    requiv - is the volumetric equivalent spherical radius in units of solar radius
    vrot   - is the rotational velocity in km/s
    mass   - is the mass in solar masses
    '''

    res = so.minimize(func_requiv_to_rpole_abs_units, np.array([requiv]), args = (vrot,requiv,mass), bounds=[(requiv*0.5,requiv*1.1)])
    return res.x[0]


def calc_critical_velocity(M, r_pole):
    '''
    M      - mass in units of solar mass
    r_pole - polar radius in units of solar radius
    Calculates critical velocity given M and R
    '''
    # We convert all values to km, kg and s so that the final rotational velocity is in km/s
    v_crit = np.sqrt(2./3. * G.to('km3/(kg s2)') * M*M_sun.to('kg') / (r_pole*R_sun.to('km')))
    return v_crit.value


def rotation_rate_to_period(v, r):
    if v == 0:
        P = 9999999999999
    else:
        P = ((2 * np.pi * r * R_sun.to('km').value) / v)/(24*60*60)
    return P


def run_cb_phoebe_model(times, abund_param_values, io_dict, run_dictionary):
    start_time_prog_1 = time.time()
    logger = phoebe.logger(clevel='ERROR')

    cb = phoebe.default_binary(contact_binary = True)

    cb.flip_constraint('pot@contact_envelope', 'requiv')
    cb.flip_constraint('fillout_factor', 'pot@contact_envelope')

    # cb['pot@contact_envelope@component'].set_value()
    cb['fillout_factor@component'].set_value(value = run_dictionary['fillout_factor'])
    cb['gravb_bol'].set_value_all(value=1.0)
    cb['irrad_frac_refl_bol'].set_value_all(value=1.0)
    cb['teff@primary'].set_value(value = run_dictionary['teff_primary'])
    cb['teff@secondary'].set_value(value = run_dictionary['teff_secondary'])
    cb['period@binary'].set_value(value = run_dictionary['period'])
    cb['sma@binary'].set_value(value = run_dictionary['sma'])
    cb['q@binary'].set_value(value = run_dictionary['q'])
    if phoebe_ver < 2.2:
        cb['ntriangles'].set_value_all(value = io_dict['ntriangles'])
    else:
        cb['ntriangles'].set_value(value = io_dict['ntriangles'])
    cb['incl'].set_value(value = run_dictionary['inclination'])

    t = list(times)

    cb.add_dataset('lc', times=t, dataset='lc01')
    cb.add_dataset('rv', times=t, dataset='rv01')
    cb.add_dataset('mesh', times=t, dataset='mesh01')
    cb.add_dataset('orb', times=t, dataset='orb01')
    if phoebe_ver < 2.2:
        cb['ld_func'].set_value_all(value = 'logarithmic')
    else:
        cb['ld_mode'].set_value_all(value = 'manual')
        cb['ld_func'].set_value_all(value = 'logarithmic')
        cb['ld_mode_bol'].set_value_all(value = 'manual')
        cb['ld_func_bol'].set_value_all(value = 'logarithmic')

    cb['atm'].set_value_all(value='blackbody')
    cb['include_times'] = 't0_ref@binary'
    cb.flip_constraint('t0_ref@binary', 't0_supconj')
    cb['t0_ref@binary@component'].set_value(value = run_dictionary['t0'])

    cb['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'xs', 'ys', 'zs', 'vxs', 'vys', 'vzs', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
    cb.run_compute()

    execution_time = time.time() - start_time_prog_1
    # print execution_time
    return cb


def run_b_phoebe_model(times, abund_param_values, io_dict, run_dictionary):
    start_time_prog_1 = time.time()
    logger = phoebe.logger(clevel='ERROR')

    b = phoebe.default_binary()

    b['gravb_bol'].set_value_all(value=1.0)
    b['irrad_frac_refl_bol'].set_value_all(value=1.0)
    b['requiv@primary'].set_value(value = run_dictionary['r_equiv_primary'])
    b['requiv@secondary'].set_value(value = run_dictionary['r_equiv_secondary'])
    b['teff@primary'].set_value(value = run_dictionary['teff_primary'])
    b['teff@secondary'].set_value(value = run_dictionary['teff_secondary'])
    b['period@binary'].set_value(value = run_dictionary['period'])
    b['sma@binary'].set_value(value = run_dictionary['sma'])
    b['q@binary'].set_value(value = run_dictionary['q'])

    b['distortion_method'].set_value_all(value = io_dict['distortion'])

    b['ntriangles'].set_value_all(value = io_dict['ntriangles'])

    b['incl@binary'].set_value(value = run_dictionary['inclination'])
    b['syncpar@primary'].set_value(value = run_dictionary['async_primary'])
    b['syncpar@secondary'].set_value(value = run_dictionary['async_secondary'])
    b['pitch@primary'].set_value(value = run_dictionary['pitch_primary'])
    b['pitch@secondary'].set_value(value = run_dictionary['pitch_secondary'])
    b['yaw@primary'].set_value(value = run_dictionary['yaw_primary'])
    b['yaw@secondary'].set_value(value = run_dictionary['yaw_secondary'])

    t = list(times)

    b.add_dataset('lc', times=t, dataset='lc01')
    b.add_dataset('rv', times=t, dataset='rv01')
    b.add_dataset('mesh', times=t, dataset='mesh01')
    b.add_dataset('orb', times=t, dataset='orb01')
    if phoebe_ver < 2.2:
        b['ld_func'].set_value_all(value = 'logarithmic')
    else:
        b['ld_mode'].set_value_all(value = 'manual')
        b['ld_func'].set_value_all(value = 'logarithmic')
        b['ld_mode_bol'].set_value_all(value = 'manual')
        b['ld_func_bol'].set_value_all(value = 'logarithmic')

    b['atm'].set_value_all(value='blackbody')
    b['include_times'] = 't0_ref@binary'
    b.flip_constraint('t0_ref@binary', 't0_supconj')
    b['t0_ref@binary@component'].set_value(value = run_dictionary['t0'])

    b['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'xs', 'ys', 'zs', 'vxs', 'vys', 'vzs', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
    b.run_compute()

    execution_time = time.time() - start_time_prog_1
    # print execution_time
    return b


def run_s_phoebe_model(times, abund_param_values, io_dict, run_dictionary):
    start_time_prog_1 = time.time()
    logger = phoebe.logger(clevel='ERROR')

    s = phoebe.default_star()

    s['teff@component'].set_value(value = run_dictionary['teff'])
    s['gravb_bol'].set_value(value = 1.0)
    s['irrad_frac_refl_bol'].set_value(value = 1.0)

    s['distortion_method'].set_value(value = io_dict['distortion'])

    s['mass@component'].set_value(value = run_dictionary['mass'])

    if run_dictionary['r_pole'] != -1:
        v_crit = calc_critical_velocity(run_dictionary['mass'], run_dictionary['r_pole'])
        if run_dictionary['v_crit_frac'] != -1:
            vrot = v_crit * run_dictionary['v_crit_frac']
            v_percent_crit = run_dictionary['v_crit_frac']
        elif run_dictionary['rotation_rate'] != -1:
            vrot = run_dictionary['rotation_rate']
            v_percent_crit = vrot / v_crit
        else:
            vrot = run_dictionary['vsini'] / (np.sin(run_dictionary['inclination'] * np.pi/180.))
            v_percent_crit = vrot / v_crit

        if vrot == 0:
            s['distortion_method'].set_value(value='sphere')
            s['requiv@component'].set_value(value = run_dictionary['r_pole'])
            s['period@component'].set_value(value = 9999999999999)
        else:
            # calculate r_equiv given r_pole and v_percent_crit:
            r_equiv, r_equator = rpole_to_requiv(run_dictionary['r_pole'], v_percent_crit, n=5000, return_r_equator=True)
            s['requiv@component'].set_value(value = r_equiv)

            # calculate period from v_rot and r_equator
            period = rotation_rate_to_period(vrot, r_equator)
            s['period@component'].set_value(value = period)

    else:
        s['requiv@component'].set_value(value = run_dictionary['requiv'])
        if run_dictionary['rotation_rate'] == 0:
            s['distortion_method'].set_value('sphere')
        elif run_dictionary['rotation_rate'] == -1 and run_dictionary['v_crit_frac'] != -1:
            # calculate r_pole given r_equiv and v_percent_crit; calc r_equator:
            r_pole = requiv_to_rpole(run_dictionary['requiv'], run_dictionary['v_crit_frac'])
            junk, r_equator = rpole_to_requiv(r_pole, run_dictionary['v_crit_frac'], n=5000, return_r_equator=True)
            # calc v_crit from mass and r_pole and then v_rot from v_crit
            v_crit = calc_critical_velocity(run_dictionary['mass'], r_pole)
            vrot = v_crit * run_dictionary['v_crit_frac']
            period = rotation_rate_to_period(vrot, r_equator)
        else:
            if run_dictionary['rotation_rate'] == -1:
                vrot = run_dictionary['vsini'] / (np.sin(run_dictionary['inclination'] * np.pi/180.))
            else:
                vrot = run_dictionary['rotation_rate']

            r_pole = requiv_to_rpole_abs_units(run_dictionary['requiv'], vrot, run_dictionary['mass'])
            v_crit = calc_critical_velocity(run_dictionary['mass'], r_pole)
            v_percent_crit = vrot / v_crit
            junk, r_equator = rpole_to_requiv(r_pole, v_percent_crit, n=5000, return_r_equator=True)
            period = rotation_rate_to_period(vrot, r_equator)
        s['period@component'].set_value(value = period)

    if run_dictionary['inclination'] == -1 and run_dictionary['rotation_rate'] == -1:
        s['incl@binary'].set_value(value = np.arcsin(run_dictionary['vsini'] / vrot) * 180./np.pi)
    elif run_dictionary['inclination'] == -1 and run_dictionary['rotation_rate'] != -1:
        s['incl@component'].set_value(value = np.arcsin(run_dictionary['vsini'] / run_dictionary['rotation_rate']) * 180./np.pi)
    else:
        s['incl@component'].set_value(value = run_dictionary['inclination'])

    s['ntriangles'].set_value(value = io_dict['ntriangles'])

    t = list(times)

    s.add_dataset('lc', times=t, dataset='lc01')
    s.add_dataset('rv', times=t, dataset='rv01')
    s.add_dataset('mesh', times=t, dataset='mesh01')
    s.add_dataset('orb', times=t, dataset='orb01')
    if phoebe_ver < 2.2:
        s['ld_func'].set_value_all(value = 'logarithmic')
    else:
        s['ld_mode'].set_value_all(value = 'manual')
        s['ld_func'].set_value_all(value = 'logarithmic')
        s['ld_mode_bol'].set_value(value = 'manual')
        s['ld_func_bol'].set_value(value = 'logarithmic')

    s['atm'].set_value(value='blackbody')
    s['include_times'] = 't0@system'
    s['t0@system'].set_value(value = run_dictionary['t0'])

    s['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'xs', 'ys', 'zs', 'vxs', 'vys', 'vzs', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
    s.run_compute()

    execution_time = time.time() - start_time_prog_1
    # print execution_time
    return s


def run_sb_phoebe_model(times, abund_param_values, io_dict, run_dictionary):
    start_time_prog_1 = time.time()
    logger = phoebe.logger(clevel='ERROR')

    b = phoebe.default_binary()

    b['teff@primary'].set_value(value = run_dictionary['teff'])
    b['gravb_bol'].set_value_all(value=1.0)
    b['irrad_frac_refl_bol'].set_value_all(value=1.0)

    b['distortion_method'].set_value_all(value = io_dict['distortion'])

    b.flip_constraint('mass@primary', 'sma@binary')
    b.flip_constraint('mass@secondary', 'q@binary')
    b['mass@component@primary'].set_value(value = run_dictionary['mass'])
    b['mass@component@secondary'].set_value(999)
    b['period@binary'].set_value(value = 99999999)

    if run_dictionary['r_pole'] != -1:
        # calculate v_{%c} from M, r_pole and v_rot:
        v_crit = calc_critical_velocity(run_dictionary['mass'], run_dictionary['r_pole'])
        if run_dictionary['v_crit_frac'] != -1:
            vrot = v_crit * run_dictionary['v_crit_frac']
            v_percent_crit = run_dictionary['v_crit_frac']
        elif run_dictionary['rotation_rate'] != -1:
            vrot = run_dictionary['rotation_rate']
            v_percent_crit = vrot / v_crit
        else:
            vrot = run_dictionary['vsini'] / (np.sin(run_dictionary['inclination'] * np.pi/180.))
            v_percent_crit = vrot / v_crit

        if vrot == 0:
            b['distortion_method'].set_value_all('sphere')
            b['requiv@primary'].set_value(value = run_dictionary['r_pole'])
        else:
            # calculate r_equiv given r_pole and v_percent_crit:
            r_equiv, r_equator = rpole_to_requiv(run_dictionary['r_pole'], v_percent_crit, n=5000, return_r_equator=True)
            b['requiv@primary'].set_value(value = r_equiv)

            # calculate period from v_rot and r_equator
            period = rotation_rate_to_period(vrot, r_equator)
            b.flip_constraint('period@primary', 'syncpar@primary')
            b['period@primary'].set_value(value = period)

    else:
        b['requiv@primary'].set_value(value = run_dictionary['requiv'])

        if run_dictionary['rotation_rate'] == 0:
            b['distortion_method'].set_value_all('sphere')
        elif run_dictionary['rotation_rate'] == -1 and run_dictionary['v_crit_frac'] != -1:
            # calculate r_pole given r_equiv and v_percent_crit; calc r_equator:
            r_pole = requiv_to_rpole(run_dictionary['requiv'], run_dictionary['v_crit_frac'])
            junk, r_equator = rpole_to_requiv(r_pole, run_dictionary['v_crit_frac'], n=5000, return_r_equator=True)
            # calc v_crit from mass and r_pole and then v_rot from v_crit
            v_crit = calc_critical_velocity(run_dictionary['mass'], r_pole)
            vrot = v_crit * run_dictionary['v_crit_frac']
            period = rotation_rate_to_period(vrot, r_equator)
        else:
            if run_dictionary['rotation_rate'] == -1:
                vrot = run_dictionary['vsini'] / (np.sin(run_dictionary['inclination'] * np.pi/180.))
            else:
                vrot = run_dictionary['rotation_rate']

            r_pole = requiv_to_rpole_abs_units(run_dictionary['requiv'], vrot, run_dictionary['mass'])
            v_crit = calc_critical_velocity(run_dictionary['mass'], r_pole)
            v_percent_crit = vrot / v_crit
            junk, r_equator = rpole_to_requiv(r_pole, v_percent_crit, n=5000, return_r_equator=True)
            period = rotation_rate_to_period(vrot, r_equator)

        b.flip_constraint('period@primary', 'syncpar@primary')
        b['period@primary'].set_value(value = period)

    if run_dictionary['inclination'] == -1 and run_dictionary['rotation_rate'] == -1:
        b['incl@binary'].set_value(value = np.arcsin(run_dictionary['vsini'] / vrot) * 180./np.pi)
    elif run_dictionary['inclination'] == -1 and run_dictionary['rotation_rate'] != -1:
        b['incl@binary'].set_value(value = np.arcsin(run_dictionary['vsini'] / run_dictionary['rotation_rate']) * 180./np.pi)
    else:
        b['incl@binary'].set_value(value = run_dictionary['inclination'])

    b['ntriangles'].set_value_all(value = io_dict['ntriangles'])

    t = list(times)

    b.add_dataset('lc', times=t, dataset='lc01')
    b.add_dataset('rv', times=t, dataset='rv01')
    b.add_dataset('mesh', times=t, dataset='mesh01')
    b.add_dataset('orb', times=t, dataset='orb01')
    if phoebe_ver < 2.2:
        b['ld_func'].set_value_all(value = 'logarithmic')
    else:
        b['ld_mode'].set_value_all(value = 'manual')
        b['ld_func'].set_value_all(value = 'logarithmic')
        b['ld_mode_bol'].set_value_all(value = 'manual')
        b['ld_func_bol'].set_value_all(value = 'logarithmic')

    b['atm'].set_value_all(value='blackbody')
    b['include_times'] = 't0_ref@binary'
    b.flip_constraint('t0_ref@binary', 't0_supconj')
    b['t0_ref@binary@component'].set_value(value = -24999999)

    b['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'xs', 'ys', 'zs', 'vxs', 'vys', 'vzs', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
    b.run_compute()

    execution_time = time.time() - start_time_prog_1
    # print execution_time
    return b


def update_output_directories(times, abund_param_values, io_dict, run_dictionary):
    model_path = io_dict['output_directory'] + 'Model_' + str(run_dictionary['run_id']).zfill(4)
    os.mkdir(model_path)
    with open(model_path + '/model_info.txt', 'w') as file:
        for key, value in io_dict.items():
            file.write('%s:%s\n' % (key, value))
        for key, value in run_dictionary.items():
            file.write('%s:%s\n' % (key, value))
    if abund_param_values['interpolate_abundances']:
        print('abundance interpolation is not supported yet.')
    he_abundances = [i for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
    cno_abundances = [j for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
    # he_abundances = [0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2]
    # cno_abundances = [6.5, 6.5, 6.5, 6.5, 7.0, 7.0, 7.0, 7.0, 7.5, 7.5, 7.5, 7.5, 8.0, 8.0, 8.0, 8.0, 8.5, 8.5, 8.5, 8.5]
    if io_dict['grid_type'] in ['FW']:
        for i in range(len(he_abundances)):
            os.mkdir(model_path + '/He' + str(he_abundances[i]) + '_CNO' + str(cno_abundances[i]))
    return model_path


def calc_spec_by_phase(mesh_vals, hjd, model_path, lines, abund_param_values, lines_dic, io_dict, run_dictionary):
    nw = []
    nf = []

    # print 'assigning spectra'
    for line in lines:
        assign_and_calc_abundance(mesh_vals, hjd, model_path, abund_param_values, lines_dic, io_dict, run_dictionary, line)


def assign_and_calc_abundance(mesh_vals, hjd, model_path, abund_param_values, lines_dic, io_dict, run_dictionary, line):
    start_time = time.time()

    if io_dict['grid_type'] in ['FW', 'FWNN']:
        he_abundances = [i for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
        cno_abundances = [j for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]

    if io_dict['grid_type'] == 'FW':
        ws, star_profs, wind_profs = assign_spectra_interp_FW(mesh_vals, line, lines_dic, io_dict, abund_param_values, run_dictionary)

        waves = []
        phots = []
        lp_bins = abund_param_values['lp_bins']
        for i in range(int(len(ws[0])/lp_bins)):
            wavg_single, phot_avg_single = calc_flux_optimize(ws[:,lp_bins*i:lp_bins*(i+1)], star_profs[:,lp_bins*i:lp_bins*(i+1)], wind_profs[:,lp_bins*i:lp_bins*(i+1)], mesh_vals)
            waves.append(wavg_single)
            phot_normalized = phot_avg_single/phot_avg_single[-1]
            if run_dictionary['v_macro'] > 0:
                phot_normalized = macroBroad(wavg_single, phot_normalized, run_dictionary['v_macro'])
            phots.append(phot_normalized)
            np.savetxt(model_path + '/He' + str(he_abundances[i]) + '_CNO' + str(cno_abundances[i]) + '/hjd' + str(round(hjd, 13)).ljust(13, '0') + '_' + line + '.txt', np.array([wavg_single, phot_normalized]).T)

    elif io_dict['grid_type'] == 'FWNN':
        for i in range(len(he_abundances)):
            ws, star_profs, wind_profs = assign_spectra_FWNN(mesh_vals, he_abundances[i], cno_abundances[i], line, lines_dic)

            wavelength, flux = calc_flux_optimize(ws, star_profs, wind_profs, mesh_vals)
            # wavg_single, phot_avg_single = calc_flux_optimize(ws[:,lp_bins*i:lp_bins*(i+1)], star_profs[:,lp_bins*i:lp_bins*(i+1)], wind_profs[:,lp_bins*i:lp_bins*(i+1)], mesh_vals)
            flux_normalized = flux/flux[-1]
            if run_dictionary['v_macro'] > 0:
                phot_normalized = macroBroad(wavelength, flux_normalized, run_dictionary['v_macro'])
            np.savetxt(model_path + '/He' + str(he_abundances[i]) + '_CNO' + str(cno_abundances[i]) + '/hjd' + str(round(hjd, 13)).ljust(13, '0') + '_' + line + '.txt', np.array([wavelength, flux_normalized]).T)
        
    elif io_dict['grid_type'] in ['K', 'T']:
        ws, star_profs, cont_profs = assign_spectra_interp_TK(mesh_vals, line, lines_dic, io_dict, run_dictionary)

        wavelength, flux, flux_cont = calc_flux_TK(ws, star_profs, cont_profs, mesh_vals)
        flux_normalized = flux/flux_cont

        if run_dictionary['v_macro'] > 0:
            phot_normalized = macroBroad(wavelength, flux_normalized, run_dictionary['v_macro'])
        
        wave_range_string = '%0.2f-%0.2f' % (line[0], line[1])
        np.savetxt(model_path + '/hjd' + str(round(hjd, 13)).ljust(13, '0') + '_' + wave_range_string + '.txt', np.array([wavelength, flux_normalized]).T)


    print(time.time() - start_time)


def macroBroad(xdata, ydata, vmacro):
    """
    Edited broadening routine from http://dx.doi.org/10.5281/zenodo.10013

      This broadens the data by a given macroturbulent velocity.
    It works for small wavelength ranges. I need to make a better
    version that is accurate for large wavelength ranges! Sorry
    for the terrible variable names, it was copied from
    convol.pro in AnalyseBstar (Karolien Lefever)
    """
    # Make the kernel
    sq_pi = np.sqrt(np.pi)
    lambda0 = np.median(xdata)
    xspacing = xdata[1] - xdata[0]
    mr = vmacro * lambda0 / c
    ccr = 2 / (sq_pi * mr)

    px = np.arange(-len(xdata) / 2, len(xdata) / 2 + 1) * xspacing
    pxmr = abs(px) / mr
    profile = ccr * (np.exp(-pxmr ** 2) + sq_pi * pxmr * (erf(pxmr) - 1.0))

    # Extend the xy axes to avoid edge-effects
    before = ydata[-profile.size / 2 + 1:]
    after = ydata[:profile.size / 2]
    extended = np.r_[before, ydata, after]

    first = xdata[0] - float(int(profile.size / 2.0 + 0.5)) * xspacing
    last = xdata[-1] + float(int(profile.size / 2.0 + 0.5)) * xspacing
    x2 = np.linspace(first, last, extended.size)

    conv_mode = "valid"

    # Do the convolution
    newydata = fftconvolve(extended, profile / profile.sum(), mode=conv_mode)

    return newydata


def assign_spectra(mesh_vals, line, lines_dic, io_dict):
    ts = np.around(mesh_vals['teffs'] / 1000.0) * 1000.0
    lgs = np.around(mesh_vals['loggs']*10.) / 10.
    rads = np.around(mesh_vals['rs'] * 4.0) / 4.0

    ws = []
    star_profs = []
    wind_profs = []
    start_time = time.time()
    for i in tqdm(range(len(ts))):
        w, st, wi = lookup_line_profs_from_dic_FW(ts[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic)
        ws.append(w)
        star_profs.append(st)
        wind_profs.append(np.array(wi))
    elapsed_time = time.time() - start_time
    # print 'Average iterations per second: ' + str(len(ts) / elapsed_time)
    ws = dopler_shift(np.array(ws), np.array([mesh_vals['rvs']]*len(ws[0])).T)
    ws=np.array(ws, dtype='float')
    return np.array(ws), np.array(star_profs), np.array(wind_profs)


def assign_spectra_interp_FW(mesh_vals, line, lines_dic, io_dict, abund_param_values, run_dictionary):
    ts = mesh_vals['ts']
    tls = mesh_vals['tls']
    tus = mesh_vals['tus']
    w1s = mesh_vals['w1s']
    w2s = mesh_vals['w2s']
    lgs = mesh_vals['lgs']
    rads = mesh_vals['rads']

    ws = []
    star_low_profs = []
    wind_low_profs = []
    star_high_profs = []
    wind_high_profs = []
    start_time = time.time()
    for i in tqdm(range(len(ts))):
        w, stl, wil = lookup_line_profs_from_dic_FW(tls[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic)
        wu, stu, wiu = lookup_line_profs_from_dic_FW(tus[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic)
        ws.append(w)
        star_low_profs.append(stl)
        wind_low_profs.append(np.array(wil))
        star_high_profs.append(stu)
        wind_high_profs.append(wiu)

    n_tot_bins = len(abund_param_values['he_abundances']) * len(abund_param_values['cno_abundances']) * abund_param_values['lp_bins']

    w1s = np.array([w1s]* n_tot_bins).T
    w2s = np.array([w2s]* n_tot_bins).T
    #When you fall directly on a grid temperature, w1==w2==0.  check for this with w3
    w3s = w1s + w2s == 0

    star_profs = np.array(star_low_profs) * w1s + np.array(star_high_profs) * w2s + np.array(star_high_profs) * w3s
    wind_profs = np.array(wind_low_profs) * w1s + np.array(wind_high_profs) * w2s + np.array(wind_high_profs) * w3s

    # Macro
    if run_dictionary['v_macro'] == -1:
        if run_dictionary['sigma_R'] > 0:
            sigma_R = np.random.normal(0, run_dictionary['sigma_R'], size=star_profs.shape[0])
        else:
            sigma_R = np.ones(star_profs.shape[0]) * run_dictionary['sigma_R']
        v_R = mesh_vals['mus'] * sigma_R

        if run_dictionary['sigma_T'] > 0:
            sigma_T = np.random.normal(0, run_dictionary['sigma_T'], size=star_profs.shape[0])
        else:
            sigma_T = np.ones(star_profs.shape[0]) * run_dictionary['sigma_T']
        theta_T = np.random.uniform(0, 2*np.pi, size=star_profs.shape[0])
        theta_mu = np.arccos(mesh_vals['mus'])
        v_T = sigma_T * np.sin(theta_mu) * np.cos(theta_T)
    else:
        v_R, v_T = 0.0, 0.0

    elapsed_time = time.time() - start_time
    # print 'Average iterations per second: ' + str(len(ts) / elapsed_time)
    # print mesh_vals['rvs']
    ws = dopler_shift(np.array(ws), np.array([mesh_vals['rvs'] + v_R + v_T]*len(ws[0])).T)
    ws=np.array(ws, dtype='float')
    return np.array(ws), np.array(star_profs), np.array(wind_profs)


def assign_spectra_FWNN(mesh_vals, he, cno, line, lines_dic):
    ts = mesh_vals['teffs']
    lgs = mesh_vals['loggs']
    rads = mesh_vals['rs']
    cnos = np.ones_like(ts) * cno
    hes = np.ones_like(ts) * he
    mus = mesh_vals['mus']

    w = lines_dic[line]['wavelength']
    ws = np.tile(w, (len(ts), 1))
    start_time = time.time()

    # ORDER OF INPUTS FOR NN: teff, logg, r, cno, he, mu
    input_array = np.array([ts, lgs, rads, cnos, hes, mus]).T
    input_array -= lines_dic['mean']
    input_array /= lines_dic['std']
    
    star_profs = 10**np.array(lines_dic[line]['phot'].predict(input_array, batch_size=10000))
    wind_profs = 10**np.array(lines_dic[line]['wind'].predict(input_array, batch_size=10000))
    wind_profs -= star_profs
    # the NN calculates the wind at its radius, but spamms calculates based on stellar radius.  we'll put it into stellar radius to avoid issues
    wind_profs /= 112**2
    # wind_profs = 10**np.array(lines_dic[line]['wind'].predict(input_array))
    # wind_profs = np.nan_to_num(wind_profs, nan=0.0)
    wind_profs[wind_profs < 0] = 0.0
    # wind_profs = np.zeros_like(wind_profs, dtype='float')
    # wind_profs = 10**np.array(lines_dic[line]['wind'].predict(input_array))
    # star_profs = np.zeros_like(wind_profs, dtype='float')


    elapsed_time = time.time() - start_time
    # print 'Average iterations per second: ' + str(len(ts) / elapsed_time)
    # print mesh_vals['rvs']
    ws = dopler_shift(np.array(ws), np.array([mesh_vals['rvs']]*len(ws[0])).T)
    ws=np.array(ws, dtype='float')
    return np.array(ws), star_profs, wind_profs


def assign_spectra_interp_TK(mesh_vals, wvrange, ranges_dic, io_dict, run_dictionary):
    ts = mesh_vals['ts']
    tls = mesh_vals['tls']
    tus = mesh_vals['tus']
    lgs = mesh_vals['lgs']
    lgls = mesh_vals['lgls']
    lgus = mesh_vals['lgus']
    w1ts = mesh_vals['w1ts']
    w2ts = mesh_vals['w2ts']
    w1gs = mesh_vals['w1gs']
    w2gs = mesh_vals['w2gs'] # in spec by phase

    ws = []
    star_lowlow_profs = []
    star_lowhigh_profs = []
    star_highlow_profs = []
    star_highhigh_profs = []
    cont_lowlow_profs = []
    cont_lowhigh_profs = []
    cont_highlow_profs = []
    cont_highhigh_profs = []

    start_time = time.time()
    for i in tqdm(range(len(ts))):
        w, stlgl, ctlgl = lookup_line_profs_from_dic_TK(tls[i], lgls[i], mesh_vals['mus'][i], mesh_vals['viss'][i], tuple(wvrange), ranges_dic, io_dict, run_dictionary)
        w, stlgu, ctlgu = lookup_line_profs_from_dic_TK(tls[i], lgus[i], mesh_vals['mus'][i], mesh_vals['viss'][i], tuple(wvrange), ranges_dic, io_dict, run_dictionary)
        w, stugl, ctugl = lookup_line_profs_from_dic_TK(tus[i], lgls[i], mesh_vals['mus'][i], mesh_vals['viss'][i], tuple(wvrange), ranges_dic, io_dict, run_dictionary)
        w, stugu, ctugu = lookup_line_profs_from_dic_TK(tus[i], lgus[i], mesh_vals['mus'][i], mesh_vals['viss'][i], tuple(wvrange), ranges_dic, io_dict, run_dictionary)
        
        ws.append(w)
        star_lowlow_profs.append(stlgl)
        star_lowhigh_profs.append(stlgu)
        star_highlow_profs.append(stugl)
        star_highhigh_profs.append(stugu)
        cont_lowlow_profs.append(ctlgl)
        cont_lowhigh_profs.append(ctlgu)
        cont_highlow_profs.append(ctugl)
        cont_highhigh_profs.append(ctugu)

    n_tot_bins = len(w)
    w1ts = np.array([w1ts]* n_tot_bins).T
    w2ts = np.array([w2ts]* n_tot_bins).T
    w1gs = np.array([w1gs]* n_tot_bins).T
    w2gs = np.array([w2gs]* n_tot_bins).T

    #When you fall directly on a grid temperature, w1==w2==0.  check for this with w3
    w3ts = w1ts + w2ts == 0
    w3gs = w1gs + w2gs == 0

    star_profs = np.array(star_lowlow_profs) * (w1ts*w1gs) + np.array(star_lowhigh_profs) * (w1ts*w2gs) + np.array(star_highlow_profs) * (w2ts*w1gs) + np.array(star_highhigh_profs) * (w2ts*w2gs) + \
        np.array(star_lowlow_profs) * (w3ts*w1gs) + np.array(star_lowhigh_profs) * (w3ts*w2gs) + \
            np.array(star_lowlow_profs) * (w1ts*w3gs) + np.array(star_highlow_profs) * (w2ts*w3gs) + \
                    np.array(star_lowlow_profs) * (w3ts*w3gs)
    
    cont_profs = np.array(cont_lowlow_profs) * (w1ts*w1gs) + np.array(cont_lowhigh_profs) * (w1ts*w2gs) + np.array(cont_highlow_profs) * (w2ts*w1gs) + np.array(cont_highhigh_profs) * (w2ts*w2gs) + \
        np.array(cont_lowlow_profs) * (w3ts*w1gs) + np.array(cont_lowhigh_profs) * (w3ts*w2gs) + \
            np.array(cont_lowlow_profs) * (w1ts*w3gs) + np.array(cont_highlow_profs) * (w2ts*w3gs) + \
                    np.array(cont_lowlow_profs) * (w3ts*w3gs)



    elapsed_time = time.time() - start_time
    # print 'Average iterations per second: ' + str(len(ts) / elapsed_time)
    # print mesh_vals['rvs']
    ws = dopler_shift(np.array(ws), np.array([mesh_vals['rvs']]*len(ws[0])).T)
    ws=np.array(ws, dtype='float')

    return np.array(ws), np.array(star_profs), np.array(cont_profs)



def dopler_shift(w, rv):
    c = 299792.458
    return w*c/(c-rv)


def lookup_line_profs_from_dic_FW(t, g, r, m, v, line, lines_dic):
    combination = 'T' + str(int(t)) + '_G' + str(g) + '_R' + format(r, '.2f')
    w = lines_dic[line]['wavelength'][combination]
    if v == 0:
        return w, np.zeros_like(w, dtype='float'), np.zeros_like(w, dtype='float')
    wlfr = 121.585278

    pray_phot = np.sqrt(1 - m**2)
    pray_wind = np.sqrt(wlfr**2 - (np.sqrt(wlfr**2 - 1)/wlfr * m * wlfr)**2)
    pray_phot_norm = pray_phot/wlfr
    pray_wind_norm = pray_wind/wlfr

    ind = int(pray_phot*100)
    indw = int(pray_wind_norm*100)

    # print filename, ind

    upper = lines_dic[line]['phot'][combination][ind+1]
    lower = lines_dic[line]['phot'][combination][ind]
    upperw = lines_dic[line]['wind'][combination][indw+1]
    lowerw = lines_dic[line]['wind'][combination][indw]

    rise = upper - lower
    risew = upperw - lowerw

    run = (pray_phot*100)%1
    runw = (pray_wind_norm*100)%1

    star_prof = lower + rise*run
    wind_prof = lowerw + risew*runw
    return w, star_prof, wind_prof


def lookup_line_profs_from_dic_TK(t, g, m, v, wvrange, ranges_dic, io_dict, run_dictionary):
    if io_dict['grid_type'] == 'K':
        metallicity_sign = 'p' if run_dictionary['metallicity'] >=0 else 'm'
        alpha_str = 'a' if run_dictionary['alpha_enhancement'] == 1 else ''
        metallicity_str = alpha_str + metallicity_sign + format(abs(run_dictionary['metallicity']), '.3f')

        combination = 'T' + str(int(t)) + '_G' + str(format(g, '.2f')) + '_M' + metallicity_str + '_V' + str(int(run_dictionary['v_micro']))
    elif io_dict['grid_type'] == 'T':
        combination = 'T' + str(int(t)) + '_G' + str(format(g, '.2f')) + '_Z' + format(abs(run_dictionary['metallicity']), '.3f') + '_V' + str(int(run_dictionary['v_micro']))

    w = ranges_dic[tuple(wvrange)]['wavelength'][combination]
    if v == 0:
        return w, np.zeros_like(w, dtype='float'), np.zeros_like(w, dtype='float')

    ind = int(m*100)

    upper_phot = ranges_dic[tuple(wvrange)]['phot'][combination][ind+1]
    lower_phot = ranges_dic[tuple(wvrange)]['phot'][combination][ind]
    upper_photcon = ranges_dic[tuple(wvrange)]['phot_cont'][combination][ind+1]
    lower_photcon = ranges_dic[tuple(wvrange)]['phot_cont'][combination][ind]

    w1 = m*100 - ind # 
    w2 = (ind+1) - m*100

    star_prof = w1*upper_phot + w2*lower_phot
    cont_prof = w1*upper_photcon + w2*lower_photcon

    return w, star_prof, cont_prof


def calc_flux_TK(ws, star_profs, cont_profs, mesh_vals):
    viss = mesh_vals['viss']
    areas = mesh_vals['areas']
    mus = mesh_vals['mus']
    Ro = 695700000

    factor_phot = viss * mus * areas / Ro**2

    star_profs *= np.array([factor_phot]*len(star_profs[0])).T
    cont_profs *= np.array([factor_phot]*len(cont_profs[0])).T


    w_min = min(ws[:,0])
    w_min = math.floor(w_min*10)/10
    w_max = max(ws[:,-1])
    w_max = math.ceil(w_max*10)/10

    wave = np.arange(w_min, w_max, 0.01)
    I_star = np.array([np.interp(wave, ws[i], star_profs[i]) for i in range(len(ws))])
    I_cont = np.array([np.interp(wave, ws[i], cont_profs[i]) for i in range(len(ws))])

    # I_star/=I_cont

    indi = np.argsort(I_star[:,0])
    indi = indi[::-1]
    flux = np.sum(I_star, axis=0)
    flux_cont = np.sum(I_cont, axis=0)
    i = 0
    while max(I_star[:,0][indi][i:]/flux[0]) > 0.05:
        i += 1
        flux = np.sum(I_star[indi][i:], axis=0)
        flux_cont = np.sum(I_cont[indi][i:], axis=0)
    return wave, flux, flux_cont


def calc_flux_optimize(ws, star_profs, wind_profs, mesh_vals):
    viss = mesh_vals['viss']
    areas = mesh_vals['areas']
    mus = mesh_vals['mus']
    rs_sol = mesh_vals['rs_sol']
    Ro = 695700000

    factor_phot = viss * mus * areas / Ro**2
    factor_wind = (mus > 0) * mus * areas / (Ro)**2 * 112**2

    star_profs *= np.array([factor_phot]*len(star_profs[0])).T
    wind_profs *= np.array([factor_wind]*len(wind_profs[0])).T

    star_profs += wind_profs


    w_min = min(ws[:,0])
    w_min = math.floor(w_min*10)/10
    w_max = max(ws[:,-1])
    w_max = math.ceil(w_max*10)/10

    wave = np.arange(w_min, w_max, 0.01)
    I_star = np.array([np.interp(wave, ws[i], star_profs[i]) for i in range(len(ws))])
    # for i in range(len(ws)):
    #     I_star.append(np.interp(wave, ws[i], star_profs[i]))

    # I_star = np.array(I_star)
    indi = np.argsort(I_star[:,0])
    indi = indi[::-1]
    flux = np.sum(I_star, axis=0)
    i = 0
    while max(I_star[:,0][indi][i:]/flux[0]) > 0.05:
        i += 1
        flux = np.sum(I_star[indi][i:], axis=0)
    return wave, flux


def calc_flux(ws, ws_all, star_profs, wind_profs, mesh_vals):
    viss = mesh_vals['viss']
    areas = mesh_vals['areas']
    mus = mesh_vals['mus']
    rs_sol = mesh_vals['rs_sol']

    w_min = min(ws_all[:,0])
    w_min = math.floor(w_min*10)/10
    w_max = max(ws_all[:,-1])
    w_max = math.ceil(w_max*10)/10

    w = ws.T
    w = np.insert(w, 0, [w_min] * len(w[0]), axis=0)
    w = np.insert(w, len(w), [w_max] * len(w[0]), axis=0)
    ws = w.T

    f = np.array(star_profs).T
    f = np.insert(f, 0, f[0], axis=0)
    f = np.insert(f, len(f), f[-1], axis=0)
    star_profs = f.T

    f = np.array(wind_profs).T
    f = np.insert(f, 0, f[0], axis=0)
    f = np.insert(f, len(f), f[-1], axis=0)
    wind_profs = f.T

    wave = np.arange(w_min, w_max, 0.01)
    I_star = []
    I_wind = []
    for i in range(len(ws)):
        I_star.append(np.interp(wave, ws[i], star_profs[i]))
        I_wind.append(np.interp(wave, ws[i], wind_profs[i]))
    factor_phot = viss * mus * areas / rs_sol**2
    factor_wind = (mus > 0) * mus * areas / (rs_sol)**2 * 120**2
    flux_phot = np.sum(I_star * np.array([factor_phot]*len(wave)).T, axis=0)
    flux_wind = np.sum(I_wind * np.array([factor_wind]*len(wave)).T, axis=0)
    flux = flux_wind + flux_phot
    return wave, flux


def apply_rad_bound(io_dict, rads, teffs, loggs):
    # grab the models in the grid
    f = glob.glob(io_dict['path_to_grid'] + '*')
    grid_points = [i.split('/')[-1] for i in f]

    # grab just the teff logg combinations:
    teff_g_combos = list(set([i.split('_R')[0] for i in grid_points]))
    teff_g_combos.sort()

    # grab the radii available for each teff logg combination
    rads_per_combo = [[float(i.split('_R')[-1]) for i in grid_points if i.startswith(j)] for j in teff_g_combos]

    # define the max and min radii dictionaries
    max_rads_dict = {}
    min_rads_dict = {}
    for i, combo in enumerate(teff_g_combos):
        max_rads_dict[combo] = max(rads_per_combo[i])
        min_rads_dict[combo] = min(rads_per_combo[i])

    # test to see if the grid is sparse or regular
    if len(list(set([value for key, value in min_rads_dict.items()]))) == 1:
        min_rads = min_rads_dict[teff_g_combos[0]]
    else:
        min_rads = np.array([min_rads_dict['T%s_G%s'%(int(teffs[i]), loggs[i])] for i in range(len(rads))])

    if len(list(set([value for key, value in max_rads_dict.items()]))) == 1:
        max_rads = max_rads_dict[teff_g_combos[0]]
    else:
        max_rads = np.array([max_rads_dict['T%s_G%s'%(int(teffs[i]), loggs[i])] for i in range(len(rads))])

    # define the new radii
    rads = rads * (rads <= max_rads) + max_rads*(rads > max_rads)
    rads = rads * (rads >= min_rads) + min_rads*(rads < min_rads)
    return rads



def spec_by_phase_cb(cb, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = cb['times@dataset@lc'].value

    if io_dict['grid_type'] == 'FW':
        combs, mode_combs = determine_tgr_combinations(cb, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)

    elif io_dict['grid_type'] == 'FWNN':
        lines_dic = line_dictionary_structure_FWNN(line_list, io_dict)

    elif io_dict['grid_type'] in ['K', 'T']:
        combs, mode_combs = determine_tgr_combinations(cb, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = wavelength_range_dictionary_structure_TK(combs, line_list, io_dict)

    
    rv_primary = cb['rvs@model@primary@rv'].value
    rv_secondary = cb['rvs@model@secondary@rv'].value

    rvs_primary_dic = {}
    rvs_secondary_dic = {}
    for i in range(len(times)):
        rvs_primary_dic[times[i]] = rv_primary[i]
        rvs_secondary_dic[times[i]] = rv_secondary[i]

    for hjd in times:
        phcb = cb['%09.6f'%hjd]
        teffs = np.concatenate([phcb['teffs@primary'].get_value(), phcb['teffs@secondary'].get_value()])
        loggs = np.concatenate([phcb['loggs@primary'].get_value(), phcb['loggs@secondary'].get_value()])
        xs = np.concatenate([phcb['us@primary'].get_value(), phcb['us@secondary'].get_value()])
        ys = np.concatenate([phcb['vs@primary'].get_value(), phcb['vs@secondary'].get_value()])
        zs = np.concatenate([phcb['ws@primary'].get_value(), phcb['ws@secondary'].get_value()])
        rvs = np.concatenate([phcb['rvs@primary@mesh'].get_value(), phcb['rvs@secondary@mesh'].get_value()])

        rvs_prim = phcb['rvs@primary@mesh'].get_value(unit=u.km/u.s)
        rvs_sec = phcb['rvs@secondary@mesh'].get_value(unit=u.km/u.s)

        # rvs_prim = phcb['vws@primary'].get_value(unit=u.km/u.s) * -1.0
        # rvs_sec = phcb['vws@secondary'].get_value(unit=u.km/u.s) * -1.0

        rv_prim_async = (rvs_prim - rvs_primary_dic[hjd]) * run_dictionary['async_primary'] + rvs_primary_dic[hjd]
        rv_sec_async = (rvs_sec - rvs_secondary_dic[hjd]) * run_dictionary['async_secondary'] + rvs_secondary_dic[hjd]

        rvs = np.concatenate([rv_prim_async, rv_sec_async])
        rvs += run_dictionary['gamma']

        # vzs = np.concatenate([phcb['vws@primary'].get_value(unit=u.km/u.s), phcb['vws@secondary'].get_value(unit=u.km/u.s)])
        # rvs = vzs * -1.0
        mus = np.concatenate([phcb['mus@primary'].get_value(), phcb['mus@secondary'].get_value()])
        viss = np.concatenate([phcb['visibilities@primary'].get_value(), phcb['visibilities@secondary'].get_value()])
        areas = np.concatenate([phcb['areas@primary'].get_value(unit=u.m**2), phcb['areas@secondary'].get_value(unit=u.m**2)])

        abs_intens = np.concatenate([phcb['abs_intensities@primary@lc01'].get_value(), phcb['abs_intensities@secondary@lc01'].get_value()])

        ldints = np.concatenate([phcb['ldint@primary@lc01'].get_value(), phcb['ldint@secondary@lc01'].get_value()])

        rs = np.concatenate([phcb['rs@primary'].get_value(), phcb['rs@secondary'].get_value()])

        rs_sol = rs * 695700000         # meters

        start_time = time.time()

        if io_dict['grid_type'] == 'FW':
            ts = np.around(np.array(teffs) / 1000.0) * 1000.0
            tls = np.floor(teffs / 1000.0) * 1000.0
            tus = np.ceil(teffs / 1000.0) * 1000.0
            w1s = (tus - teffs)/1000.0
            w2s = (teffs - tls)/1000.0
            lgs = np.around(loggs*10.) / 10.
            rads = np.around(rs * 4.0) / 4.0

            if io_dict['rad_bound']:
                rads = apply_rad_bound(io_dict, rads, ts, lgs)

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1s':w1s, 'w2s':w2s, 'lgs':lgs, 'rads':rads}

        elif io_dict['grid_type'] == 'FWNN':
            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        elif io_dict['grid_type'] == 'K':
            ts = np.where(teffs <= 13000, np.around(teffs / 250) * 250, np.around(teffs / 1000) * 1000)
            tls = np.where(teffs <= 13000, np.floor(teffs / 250) * 250, np.floor(teffs / 1000) * 1000)
            tus = np.where(teffs <= 13000, np.ceil(teffs / 250) * 250, np.ceil(teffs / 1000) * 1000)

            lgs = np.around(loggs / 0.5) * 0.5
            lgls = np.floor(loggs / 0.5) * 0.5
            lgus = np.ceil(loggs / 0.5) * 0.5

            w1ts = np.where(ts <= 13000, (tus - ts)/250.0, (tus - ts)/1000.0)
            w2ts = np.where(ts <= 13000, (ts - tls)/250.0, (ts - tls)/1000.0)
            w1gs = (lgus - lgs)/0.5
            w2gs = (lgs - lgls)/0.5

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}
        
        elif io_dict['grid_type'] == 'T':
            ts = np.where(teffs < 27500, np.around(teffs / 1000) * 1000, np.around(teffs / 2500) * 2500)
            tls = np.where(teffs < 27500, np.floor(teffs / 1000) * 1000, np.floor(teffs / 2500) * 2500)
            tus = np.where(teffs < 27500, np.ceil(teffs / 1000) * 1000, np.ceil(teffs / 2500) * 2500)

            lgs = np.around(loggs / 0.25) * 0.25
            lgls = np.floor(loggs / 0.25) * 0.25
            lgus = np.ceil(loggs / 0.25) * 0.25

            w1ts = np.where(ts < 27500, (tus - ts)/1000.0, (tus - ts)/2500.0)
            w2ts = np.where(ts < 27500, (ts - tls)/1000.0, (ts - tls)/2500.0)
            w1gs = (lgus - lgs)/0.25
            w2gs = (lgs - lgls)/0.25

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}


        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict, run_dictionary)
        # print time.time() - start_time


def spec_by_phase_b(b, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = b['times@dataset@lc'].value

    if io_dict['grid_type'] == 'FW':
        combs, mode_combs = determine_tgr_combinations(b, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)

    elif io_dict['grid_type'] == 'FWNN':
        lines_dic = line_dictionary_structure_FWNN(line_list, io_dict)

    elif io_dict['grid_type'] in ['K', 'T']:
        combs, mode_combs = determine_tgr_combinations(b, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = wavelength_range_dictionary_structure_TK(combs, line_list, io_dict)

    
    rv_primary = b['rvs@model@primary@rv'].value
    rv_secondary = b['rvs@model@secondary@rv'].value

    rvs_primary_dic = {}
    rvs_secondary_dic = {}
    for i in range(len(times)):
        rvs_primary_dic[times[i]] = rv_primary[i]
        rvs_secondary_dic[times[i]] = rv_secondary[i]

    for hjd in times:
        phb = b['%09.6f'%hjd]
        teffs = np.concatenate([phb['teffs@primary'].get_value(), phb['teffs@secondary'].get_value()])
        loggs = np.concatenate([phb['loggs@primary'].get_value(), phb['loggs@secondary'].get_value()])
        xs = np.concatenate([phb['us@primary'].get_value(), phb['us@secondary'].get_value()])
        ys = np.concatenate([phb['vs@primary'].get_value(), phb['vs@secondary'].get_value()])
        zs = np.concatenate([phb['ws@primary'].get_value(), phb['ws@secondary'].get_value()])
        rvs = np.concatenate([phb['rvs@primary@mesh'].get_value(), phb['rvs@secondary@mesh'].get_value()])

        # rvs_prim = phb['vws@primary'].get_value(unit=u.km/u.s) * -1.0
        # rvs_sec = phb['vws@secondary'].get_value(unit=u.km/u.s) * -1.0
        # rvs = np.concatenate([rvs_prim, rvs_sec])
        rvs += run_dictionary['gamma']

        mus = np.concatenate([phb['mus@primary'].get_value(), phb['mus@secondary'].get_value()])
        viss = np.concatenate([phb['visibilities@primary'].get_value(), phb['visibilities@secondary'].get_value()])
        areas = np.concatenate([phb['areas@primary'].get_value(unit=u.m**2), phb['areas@secondary'].get_value(unit=u.m**2)])

        abs_intens = np.concatenate([phb['abs_intensities@primary@lc01'].get_value(), phb['abs_intensities@secondary@lc01'].get_value()])

        ldints = np.concatenate([phb['ldint@primary@lc01'].get_value(), phb['ldint@secondary@lc01'].get_value()])

        rs = np.concatenate([phb['rs@primary'].get_value(), phb['rs@secondary'].get_value()])

        rs_sol = rs * 695700000         # meters

        start_time = time.time()

        if io_dict['grid_type'] == 'FW':
            ts = np.around(np.array(teffs) / 1000.0) * 1000.0
            tls = np.floor(teffs / 1000.0) * 1000.0
            tus = np.ceil(teffs / 1000.0) * 1000.0
            w1s = (tus - teffs)/1000.0
            w2s = (teffs - tls)/1000.0
            lgs = np.around(loggs*10.) / 10.
            rads = np.around(rs * 4.0) / 4.0

            if io_dict['rad_bound']:
                rads = apply_rad_bound(io_dict, rads, ts, lgs)

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1s':w1s, 'w2s':w2s, 'lgs':lgs, 'rads':rads}

        elif io_dict['grid_type'] == 'FWNN':
            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        elif io_dict['grid_type'] == 'K':
            ts = np.where(teffs <= 13000, np.around(teffs / 250) * 250, np.around(teffs / 1000) * 1000)
            tls = np.where(teffs <= 13000, np.floor(teffs / 250) * 250, np.floor(teffs / 1000) * 1000)
            tus = np.where(teffs <= 13000, np.ceil(teffs / 250) * 250, np.ceil(teffs / 1000) * 1000)

            lgs = np.around(loggs / 0.5) * 0.5
            lgls = np.floor(loggs / 0.5) * 0.5
            lgus = np.ceil(loggs / 0.5) * 0.5

            w1ts = np.where(ts <= 13000, (tus - ts)/250.0, (tus - ts)/1000.0)
            w2ts = np.where(ts <= 13000, (ts - tls)/250.0, (ts - tls)/1000.0)
            w1gs = (lgus - lgs)/0.5
            w2gs = (lgs - lgls)/0.5

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}
        
        elif io_dict['grid_type'] == 'T':
            ts = np.where(teffs < 27500, np.around(teffs / 1000) * 1000, np.around(teffs / 2500) * 2500)
            tls = np.where(teffs < 27500, np.floor(teffs / 1000) * 1000, np.floor(teffs / 2500) * 2500)
            tus = np.where(teffs < 27500, np.ceil(teffs / 1000) * 1000, np.ceil(teffs / 2500) * 2500)

            lgs = np.around(loggs / 0.25) * 0.25
            lgls = np.floor(loggs / 0.25) * 0.25
            lgus = np.ceil(loggs / 0.25) * 0.25

            w1ts = np.where(ts < 27500, (tus - ts)/1000.0, (tus - ts)/2500.0)
            w2ts = np.where(ts < 27500, (ts - tls)/1000.0, (ts - tls)/2500.0)
            w1gs = (lgus - lgs)/0.25
            w2gs = (lgs - lgls)/0.25

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}


        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict, run_dictionary)
        # print time.time() - start_time


def spec_by_phase_s(s, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = s['times@dataset@lc'].value

    if io_dict['grid_type'] == 'FW':
        combs, mode_combs = determine_tgr_combinations(s, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)
    
    elif io_dict['grid_type'] == 'FWNN':
        lines_dic = line_dictionary_structure_FWNN(line_list, io_dict)
    
    elif io_dict['grid_type'] in ['K', 'T']:
        combs, mode_combs = determine_tgr_combinations(s, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = wavelength_range_dictionary_structure_TK(combs, line_list, io_dict)


    for hjd in times:
        s_t = s['%09.6f'%hjd]
        teffs = s_t['teffs'].get_value()
        loggs = s_t['loggs'].get_value()
        xs = s_t['us'].get_value()
        ys = s_t['vs'].get_value()
        zs = s_t['ws'].get_value()
        rvs = s_t['rvs@mesh'].get_value()

        # vzs = s_t['vws'].get_value(unit=u.km/u.s)
        # rvs = vzs * -1.0
        if run_dictionary['rotation_rate'] == 0:
            rvs = np.zeros_like(rvs)
            # print('zeroed')
        rvs += run_dictionary['gamma']
        mus = s_t['mus'].get_value()
        viss = s_t['visibilities'].get_value()
        areas = s_t['areas'].get_value(unit=u.m**2)

        abs_intens = s_t['abs_intensities@lc01'].get_value()

        ldints = s_t['ldint@lc01'].get_value()

        rs = s_t['rs'].get_value()

        rs_sol = rs * 695700000         # meters

        start_time = time.time()

        if io_dict['grid_type'] == 'FW':
            ts = np.around(np.array(teffs) / 1000.0) * 1000.0
            tls = np.floor(teffs / 1000.0) * 1000.0
            tus = np.ceil(teffs / 1000.0) * 1000.0
            w1s = (tus - teffs)/1000.0
            w2s = (teffs - tls)/1000.0
            lgs = np.around(loggs*10.) / 10.
            rads = np.around(rs * 4.0) / 4.0

            if io_dict['rad_bound']:
                rads = apply_rad_bound(io_dict, rads, ts, lgs)

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1s':w1s, 'w2s':w2s, 'lgs':lgs, 'rads':rads}

        elif io_dict['grid_type'] == 'FWNN':
            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        elif io_dict['grid_type'] == 'K':
            ts = np.where(teffs <= 13000, np.around(teffs / 250) * 250, np.around(teffs / 1000) * 1000)
            tls = np.where(teffs <= 13000, np.floor(teffs / 250) * 250, np.floor(teffs / 1000) * 1000)
            tus = np.where(teffs <= 13000, np.ceil(teffs / 250) * 250, np.ceil(teffs / 1000) * 1000)

            lgs = np.around(loggs / 0.5) * 0.5
            lgls = np.floor(loggs / 0.5) * 0.5
            lgus = np.ceil(loggs / 0.5) * 0.5

            w1ts = np.where(ts <= 13000, (tus - ts)/250.0, (tus - ts)/1000.0)
            w2ts = np.where(ts <= 13000, (ts - tls)/250.0, (ts - tls)/1000.0)
            w1gs = (lgus - lgs)/0.5
            w2gs = (lgs - lgls)/0.5

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}
        
        elif io_dict['grid_type'] == 'T':
            ts = np.where(teffs < 27500, np.around(teffs / 1000) * 1000, np.around(teffs / 2500) * 2500)
            tls = np.where(teffs < 27500, np.floor(teffs / 1000) * 1000, np.floor(teffs / 2500) * 2500)
            tus = np.where(teffs < 27500, np.ceil(teffs / 1000) * 1000, np.ceil(teffs / 2500) * 2500)

            lgs = np.around(loggs / 0.25) * 0.25
            lgls = np.floor(loggs / 0.25) * 0.25
            lgus = np.ceil(loggs / 0.25) * 0.25

            w1ts = np.where(ts < 27500, (tus - ts)/1000.0, (tus - ts)/2500.0)
            w2ts = np.where(ts < 27500, (ts - tls)/1000.0, (ts - tls)/2500.0)
            w1gs = (lgus - lgs)/0.25
            w2gs = (lgs - lgls)/0.25

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}


        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict, run_dictionary)
        # print time.time() - start_time


def spec_by_phase_sb(s, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = s['times@dataset@lc'].value

    if io_dict['grid_type'] == 'FW':
        combs, mode_combs = determine_tgr_combinations(s, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)
        
        lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)
    
    elif io_dict['grid_type'] == 'FWNN':
        lines_dic = line_dictionary_structure_FWNN(line_list, io_dict)

    elif io_dict['grid_type'] in ['K', 'T']:
        combs, mode_combs = determine_tgr_combinations(s, io_dict, run_dictionary)

        # check to see if the grid is complete
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]
        missing_combs = [i for i in combs if i not in grid_entries]
        if len(missing_combs) > 0:
            print('WARNING: input grid entries missing.')
            print(missing_combs)

        lines_dic = wavelength_range_dictionary_structure_TK(combs, line_list, io_dict)

    for hjd in times:
        s_t = s['%09.6f'%hjd]
        if io_dict['gravity_darkening'] == 'EL':
            teffs = Espinosa_Lara_2011_gd_grid(s, s_t, run_dictionary)
        else:
            teffs = s_t['teffs@primary'].get_value()
        loggs = s_t['loggs@primary'].get_value()
        xs = s_t['us@primary'].get_value()
        ys = s_t['vs@primary'].get_value()
        zs = s_t['ws@primary'].get_value()
        rvs = s_t['rvs@primary@mesh'].get_value()

        # vzs = s_t['vws'].get_value(unit=u.km/u.s)
        # rvs = vzs * -1.0
        if run_dictionary['rotation_rate'] == 0:
            rvs = np.zeros_like(rvs)
            # print('zeroed')
        rvs += run_dictionary['gamma']
        mus = s_t['mus@primary'].get_value()
        viss = s_t['visibilities@primary'].get_value()
        areas = s_t['areas@primary'].get_value(unit=u.m**2)

        abs_intens = s_t['abs_intensities@primary@lc01'].get_value()

        ldints = s_t['ldint@primary@lc01'].get_value()

        rs = s_t['rs@primary'].get_value()

        rs_sol = rs * 695700000         # meters

        start_time = time.time()

        if io_dict['grid_type'] == 'FW':
            ts = np.around(np.array(teffs) / 1000.0) * 1000.0
            tls = np.floor(teffs / 1000.0) * 1000.0
            tus = np.ceil(teffs / 1000.0) * 1000.0
            w1s = (tus - teffs)/1000.0
            w2s = (teffs - tls)/1000.0
            lgs = np.around(loggs*10.) / 10.
            rads = np.around(rs * 4.0) / 4.0

            if io_dict['rad_bound']:
                rads = apply_rad_bound(io_dict, rads, ts, lgs)

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1s':w1s, 'w2s':w2s, 'lgs':lgs, 'rads':rads}

        elif io_dict['grid_type'] == 'FWNN':
            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        elif io_dict['grid_type'] == 'K':
            ts = np.where(teffs <= 13000, np.around(teffs / 250) * 250, np.around(teffs / 1000) * 1000)
            tls = np.where(teffs <= 13000, np.floor(teffs / 250) * 250, np.floor(teffs / 1000) * 1000)
            tus = np.where(teffs <= 13000, np.ceil(teffs / 250) * 250, np.ceil(teffs / 1000) * 1000)

            lgs = np.around(loggs / 0.5) * 0.5
            lgls = np.floor(loggs / 0.5) * 0.5
            lgus = np.ceil(loggs / 0.5) * 0.5

            w1ts = np.where(ts <= 13000, (tus - ts)/250.0, (tus - ts)/1000.0)
            w2ts = np.where(ts <= 13000, (ts - tls)/250.0, (ts - tls)/1000.0)
            w1gs = (lgus - lgs)/0.5
            w2gs = (lgs - lgls)/0.5

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}
        
        elif io_dict['grid_type'] == 'T':
            ts = np.where(teffs < 27500, np.around(teffs / 1000) * 1000, np.around(teffs / 2500) * 2500)
            tls = np.where(teffs < 27500, np.floor(teffs / 1000) * 1000, np.floor(teffs / 2500) * 2500)
            tus = np.where(teffs < 27500, np.ceil(teffs / 1000) * 1000, np.ceil(teffs / 2500) * 2500)

            lgs = np.around(loggs / 0.25) * 0.25
            lgls = np.floor(loggs / 0.25) * 0.25
            lgus = np.ceil(loggs / 0.25) * 0.25

            w1ts = np.where(ts < 27500, (tus - ts)/1000.0, (tus - ts)/2500.0)
            w2ts = np.where(ts < 27500, (ts - tls)/1000.0, (ts - tls)/2500.0)
            w1gs = (lgus - lgs)/0.25
            w2gs = (lgs - lgls)/0.25

            mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol, 'ts':ts, 'tls':tls, 'tus':tus, 'w1ts':w1ts, 'w2ts':w2ts, 'lgs':lgs, 'lgls':lgls, 'lgus':lgus, 'w1gs':w1gs, 'w2gs':w2gs}


        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict, run_dictionary)
        # print time.time() - start_time


def cartesian_to_spherical(x, y, z, x_c = 0, y_c = 0, z_c = 0):
    x -= x_c
    y -= y_c
    z -= z_c

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.abs(z) / r)
    theta -= (theta > np.pi) * np.pi
    phi = np.arctan(y/x)

    return r, theta, phi


def interpolate_psi_grid(v_linear_crit_frac):
    psi_grid = np.load('psi_grid.npy')
    v_linear_percent_crit = v_linear_crit_frac * 100

    if (v_linear_percent_crit).is_integer():
        psis = psi_grid[int(v_linear_percent_crit)]
    else:
        v_upper = int(np.ceil(v_linear_percent_crit))
        v_lower = int(np.floor(v_linear_percent_crit))

        weight_lower = v_upper - v_linear_percent_crit
        weight_upper = v_linear_percent_crit - v_lower

        psis = psi_grid[v_upper] * weight_upper + psi_grid[v_lower] * weight_lower
    return psis


def Espinosa_Lara_2011_gd_grid(s_full, s, run_dictionary):

    loggs = s['loggs@primary'].get_value()
    # grab total surface area
    areas = s['areas@primary'].get_value()
    area = np.sum(areas)

    # grab r, theta, phi
    x = s['xs@primary'].get_value()
    y = s['ys@primary'].get_value()
    z = s['zs@primary'].get_value()

    # x_c = s_full['xs@orb@primary'].get_value()[0]
    # y_c = s_full['zs@orb@primary'].get_value()[0]
    # z_c = s_full['ys@orb@primary'].get_value()[0]

    x_c = 0
    y_c = 0
    z_c = 0

    rs, thetas, phis = cartesian_to_spherical(x, y, z, x_c, y_c, z_c)

    # calculate psi
    psis_grid = interpolate_psi_grid(run_dictionary['v_crit_frac'])
    thetas_grid = np.load('theta_grid.npy')
    psis = np.interp(thetas, thetas_grid, psis_grid)

#     Ts = (T_eff**4 * area / (4 * np.pi * r_equator**2))**(1./4.) * (1/r_prime**4 + w**4 * r_prime**2 * np.sin(thetas)**2 - 2 * w**2 * np.sin(thetas)**2 / r_prime)**(1./8.) * np.sqrt(np.tan(psis)/np.tan(thetas))
    Ts = (run_dictionary['teff']**4 * area / (4 * np.pi * c.G.to('solRad3/(solMass s2)').value * run_dictionary['mass']))**(1./4.) * np.sqrt(np.tan(psis)/np.tan(thetas)) * (10**loggs / c.R_sun.to('cm').value)**(1./4.)

    return Ts


def determine_tgr_combinations(cb, io_dict, run_dictionary):
    times = cb['times@dataset@lc'].value
    teffs = []
    loggs = []
    rs = []
    if io_dict['object_type'] == 'contact_binary' or io_dict['object_type'] == 'binary':
        for i in times:
            teff = np.concatenate([cb['teffs@primary@%09.6f'%i].get_value(), cb['teffs@secondary@%09.6f'%i].get_value()])
            logg = np.concatenate([cb['loggs@primary@%09.6f'%i].get_value(), cb['loggs@secondary@%09.6f'%i].get_value()])
            r = np.concatenate([cb['rs@primary@%09.6f'%i].get_value(), cb['rs@secondary@%09.6f'%i].get_value()])
            teffs.extend(teff)
            loggs.extend(logg)
            rs.extend(r)
    elif io_dict['object_type'] == 'single':
        if io_dict['distortion'] in ['rotstar', 'sphere']:
            # if len(times) > 1:
            for i in times:
                phcb = cb['%09.6f'%i]
                teff = phcb['teffs'].get_value()
                logg = phcb['loggs'].get_value()
                r = phcb['rs'].get_value()
                teffs.extend(teff)
                loggs.extend(logg)
                rs.extend(r)
            # else:
            #     teffs.extend(cb['teffs'].get_value())
            #     loggs.extend(cb['loggs'].get_value())
            #     rs.extend(cb['rs'].get_value())
        else:
            for i in times:
                phcb = cb['%09.6f'%i]
                logg = phcb['loggs@primary'].get_value()
                r = phcb['rs@primary'].get_value()
                if io_dict['gravity_darkening'] == 'EL':
                    teff = Espinosa_Lara_2011_gd_grid(cb, phcb, run_dictionary)
                else:
                    teff = phcb['teffs@primary'].get_value()
                # teff = phcb['teffs@primary'].get_value()
                teffs.extend(teff)
                loggs.extend(logg)
                rs.extend(r)

    teffs = np.array(teffs)
    loggs = np.array(loggs)
    rs = np.array(rs)

    if io_dict['grid_type'] == 'FW':
        ts = np.around(teffs / 1000.0) * 1000.0
        tls = np.floor(teffs / 1000.0) * 1000.0
        tus = np.ceil(teffs / 1000.0) * 1000.0
        lgs = np.around(loggs*10.) / 10.
        lgls = np.floor(loggs / 0.1) * 0.1
        lgus = np.ceil(loggs / 0.1) * 0.1
        rads = np.around(rs * 4.0) / 4.0
        if io_dict['rad_bound']:
            rads = apply_rad_bound(io_dict, rads, ts, lgs)

        combinations = ['T' + str(int(tls[i])) + '_G' + format(lgs[i], '.1f') + '_R' + format(rads[i], '.2f') for i in range(len(ts))]
        combinations.extend(['T' + str(int(tls[i])) + '_G' + format(lgls[i], '.1f') + '_R' + format(rads[i], '.2f') for i in range(len(ts))])
        combinations.extend(['T' + str(int(tls[i])) + '_G' + format(lgus[i], '.1f') + '_R' + format(rads[i], '.2f') for i in range(len(ts))])
        combinations.extend(['T' + str(int(tus[i])) + '_G' + format(lgls[i], '.1f') + '_R' + format(rads[i], '.2f') for i in range(len(ts))])
        combinations.extend(['T' + str(int(tus[i])) + '_G' + format(lgus[i], '.1f') + '_R' + format(rads[i], '.2f') for i in range(len(ts))])
        
    
    elif io_dict['grid_type'] == 'K':
        ts = np.where(teffs <= 13000, np.around(teffs / 250) * 250, np.around(teffs / 1000) * 1000)
        tls = np.where(teffs <= 13000, np.floor(teffs / 250) * 250, np.floor(teffs / 1000) * 1000)
        tus = np.where(teffs <= 13000, np.ceil(teffs / 250) * 250, np.ceil(teffs / 1000) * 1000)

        lgs = np.around(loggs / 0.5) * 0.5
        lgls = np.floor(loggs / 0.5) * 0.5
        lgus = np.ceil(loggs / 0.5) * 0.5

        metallicity_sign = 'p' if run_dictionary['metallicity'] >=0 else 'm'
        alpha_str = 'a' if run_dictionary['alpha_enhancement'] == 1 else ''
        metallicity_str = alpha_str + metallicity_sign + format(abs(run_dictionary['metallicity']), '.3f')

        combinations = ['T' + str(int(ts[i])) + '_G' + format(lgs[i], '.2f') + '_M' + metallicity_str + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))]
        combinations.extend(['T' + str(int(tls[i])) + '_G' + format(lgls[i], '.2f') + '_M' + metallicity_str + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])
        combinations.extend(['T' + str(int(tls[i])) + '_G' + format(lgus[i], '.2f') + '_M' + metallicity_str + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])
        combinations.extend(['T' + str(int(tus[i])) + '_G' + format(lgls[i], '.2f') + '_M' + metallicity_str + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])
        combinations.extend(['T' + str(int(tus[i])) + '_G' + format(lgus[i], '.2f') + '_M' + metallicity_str + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))]) #f"{lgus[i]:.1f}"


    # TLUSTY
    elif io_dict['grid_type'] == 'T':
        ts = np.where(teffs < 27500, np.around(teffs / 1000) * 1000, np.around(teffs / 2500) * 2500)
        tls = np.where(teffs < 27500, np.floor(teffs / 1000) * 1000, np.floor(teffs / 2500) * 2500)
        tus = np.where(teffs < 27500, np.ceil(teffs / 1000) * 1000, np.ceil(teffs / 2500) * 2500)

        lgs = np.around(loggs / 0.25) * 0.25
        lgls = np.floor(loggs / 0.25) * 0.25
        lgus = np.ceil(loggs / 0.25) * 0.25

        combinations = ['T' + str(int(ts[i])) + '_G' + format(lgs[i], '.2f') + '_Z' + format(abs(run_dictionary['metallicity']), '.3f') + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))]
        combinations.extend(['T' + str(int(tls[i])) + '_G' + format(lgls[i], '.2f') + '_Z' + format(abs(run_dictionary['metallicity']), '.3f') + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])
        combinations.extend(['T' + str(int(tls[i])) + '_G' + format(lgus[i], '.2f') + '_Z' + format(abs(run_dictionary['metallicity']), '.3f') + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])
        combinations.extend(['T' + str(int(tus[i])) + '_G' + format(lgls[i], '.2f') + '_Z' + format(abs(run_dictionary['metallicity']), '.3f') + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])
        combinations.extend(['T' + str(int(tus[i])) + '_G' + format(lgus[i], '.2f') + '_Z' + format(abs(run_dictionary['metallicity']), '.3f') + '_V' + str(int(run_dictionary['v_micro'])) for i in range(len(ts))])

    return list(set(combinations)), max(set(combinations), key=combinations.count)



def line_dictionary_structure(combinations, lines, io_dict):
    lines_dic = {}
    for line in lines:
        wl = {}
        wi = {}
        ph = {}
        for i in trange(len(combinations), desc='line_dic', leave=False):
            # for comb in combinations:
            filename = io_dict['path_to_grid'] + combinations[i] + '/' + line
            w = np.load(filename + '_wl.npy')
            wind = np.load(filename + 'wind_101.npy')
            phot = np.load(filename + 'phot_101.npy')
            wl[combinations[i]] = w
            wi[combinations[i]] = wind
            ph[combinations[i]] = phot
        line_dic = {'wavelength': wl, 'wind':wi, 'phot':ph}
        lines_dic[line] = line_dic
    return lines_dic


def line_dictionary_structure_FWNN(lines, io_dict):
    import keras

    lines_dic = {}
    for line in lines:
        filepath = io_dict['path_to_grid'] + line + '/'
        wl = np.load(filepath + 'wnew_' + line + '.npy')
        wi = keras.saving.load_model(filepath + 'winds_%s_model.keras'%line)
        ph = keras.saving.load_model(filepath + 'phots_%s_model.keras'%line)
        line_dic = {'wavelength': wl, 'wind':wi, 'phot':ph}
        lines_dic[line] = line_dic
    mean, std = np.loadtxt(io_dict['path_to_grid'] + 'norm_array.txt')
    lines_dic['mean'] = mean
    lines_dic['std'] = std
    return lines_dic


def wavelength_range_dictionary_structure_TK(combinations, wavelength_dict, io_dict):
    ranges_dic = {}

    filepath = io_dict['path_to_grid']

    for wvrange in wavelength_dict:
        wl = {}
        ph = {}
        ph_cont = {}

        for combination in combinations:
            # w = load_tar_npy(filepath + combination + '.tar.gz', '_wave.npy')
            # phot = load_tar_npy(filepath + combination + '.tar.gz', '_intens.npy')
            # phot_cont = load_tar_npy(filepath + combination + '.tar.gz', '_contintens.npy')
            w = np.load(filepath + combination + '/' + combination + '_wave.npy')
            phot = np.load(filepath + combination + '/' + combination + '_intens.npy')
            phot_cont = np.load(filepath + combination + '/' + combination + '_contintens.npy')

            mask_range = (w >= wvrange[0]) & (w <= wvrange[1])
            w_range = w[mask_range]
            phot_range = phot[:, mask_range]
            phot_cont_range = phot_cont[:, mask_range]

            wl[combination] = w_range
            ph[combination] = phot_range
            ph_cont[combination] = phot_cont_range

        range_dic = {'wavelength': wl, 'phot': ph, 'phot_cont': ph_cont}
        ranges_dic[tuple(wvrange)] = range_dic

    return ranges_dic

'''
def load_tar_npy(tarball, file_extension):
    with tarfile.open(tarball, 'r') as tar:
        array_file = BytesIO()
        array_file.write(tar.extractfile('./%s%s' %(tarball.split('/')[-1][:-7], file_extension)).read())
        array_file.seek(0)
        w = np.load(array_file, allow_pickle=True)
    return w
'''

def interp_line_dictionary_structure_new(combinations, lines, io_dict, mode_combs, abund_param_values):
    lines_dic = {}
    lp_bins = abund_param_values['lp_bins']
    n_abundance_combinations = len(abund_param_values['he_abundances']) * len(abund_param_values['cno_abundances'])
    for line in lines:
        wl = {}
        wi = {}
        ph = {}
        filename = io_dict['path_to_grid'] + mode_combs + '/' + line
        w_ref = np.load(filename + '_wl.npy')
        for i in trange(len(combinations), desc='interp_dict', leave=False):
            filename = io_dict['path_to_grid'] + combinations[i] + '/' + line
            w = np.load(filename + '_wl.npy')
            wind = np.load(filename + 'wind_101.npy')
            phot = np.load(filename + 'phot_101.npy')
            winds = []
            phots = []
            for j in range(n_abundance_combinations):
                w_low = w_ref[lp_bins*j:lp_bins*(j+1)]

                w_high = w[lp_bins*j:lp_bins*(j+1)]
                w_high = np.insert(np.insert(w_high, 0, 0), len(w_high)+1, 99999)

                wind_high = wind[:,lp_bins*j:lp_bins*(j+1)]
                wind_high = np.insert(wind_high, 0, wind_high.T[0], axis=1)
                wind_high = np.insert(wind_high, -1, wind_high.T[-1], axis=1)

                phot_high = phot[:,lp_bins*j:lp_bins*(j+1)]
                phot_high = np.insert(phot_high, 0, phot_high.T[0], axis=1)
                phot_high = np.insert(phot_high, -1, phot_high.T[-1], axis=1)

                wind_new = [np.interp(w_low, w_high, wind_high[k]) for k in range(len(wind_high))]
                phot_new = [np.interp(w_low, w_high, phot_high[k]) for k in range(len(phot_high))]

                winds.extend(np.array(wind_new).T)
                phots.extend(np.array(phot_new).T)

            wl[combinations[i]] = w_ref
            wi[combinations[i]] = np.array(winds).T
            ph[combinations[i]] = np.array(phots).T
        line_dic = {'wavelength': wl, 'wind':wi, 'phot':ph}
        lines_dic[line] = line_dic
    return lines_dic


def PFGS_checks(io_dict, times, line_list):
    mods = glob.glob(io_dict['output_directory'] + '*')
    for mod in mods:
        abunds = glob.glob(mod + '/*')
        for abund in abunds:
            out_lines = glob.glob(abund + '/*')
            if len(out_lines) < (len(times) * len(line_list)):
                print('some models failed to run in %s' %abund)


def calc_chi2(obs, exp, w = None, lb = None):
    obs = np.array(obs)
    exp = np.array(exp)
    if w is not None:
        inds = [i for i in range(len(w)) if lb[0] <= w[i] and w[i] <= lb[1]]
        chi2 = np.sum((obs[inds] - exp[inds])**2 / exp[inds])
    else:
        chi2 = np.sum((obs - exp)**2 / exp)
    return chi2


def correct_obs_exp(obs_wave, obs_flux, exp_wave, exp_flux):
    min_w = min(exp_wave)
    max_w = max(exp_wave)
    wavelength_corrected = np.array([i for i in obs_wave if min_w <= i <= max_w])
    obs_flux_corrected = np.interp(wavelength_corrected, obs_wave, obs_flux)
    exp_flux_corrected = np.interp(wavelength_corrected, exp_wave, exp_flux)
    return wavelength_corrected, obs_flux_corrected, exp_flux_corrected


def fw_stitch(w1, f1, w2, f2):
    # To be used with Fastwind Spectral lines.  This function stitches together spectral lines by adding the
    # normalized fluxes (minus 1) and then renormalizing (adding 1 back).  If there is no overlap region then
    # the lines are just appended
    finwave = []
    finflux = []
    if w1[-1] < w2[0]:
        #if no overlap then
        finwave.extend(w1)
        finwave.extend(w2)
        finflux.extend(f1)
        finflux.extend(f2)
    else:
        # determine overlap and then use the resample_linear function to get the two sections to have the same
        # wavelength array  Then multiply the fluxes together and return final wavelength and flux arrays
        woverlap = [i for i in w1 if i >= w2[0]]
        woverlap.extend([i for i in w2 if i <= w1[-1]])
        woverlap = list(set(woverlap))
        woverlap.sort()
        finwave = [i for i in w1 if i <w2[0]]
        finwave.extend(woverlap)
        finwave.extend([i for i in w2 if i > w1[-1]])
        finflux = [f1[i] for i in range(len(w1)) if w1[i] < w2[0]]
        nf1 = np.interp(woverlap, w1, f1)
        nf2 = np.interp(woverlap, w2, f2)
        #nf1 = resample_linear(w1, f1, woverlap)
        #nf2 = resample_linear(w2, f2, woverlap)
        foverlap = np.array(nf1) * np.array(nf2)
        finflux.extend(foverlap)
        finflux.extend([f2[i] for i in range(len(w2)) if w2[i] > w1[-1]])
    return finwave, finflux


def calc_chi2_per_model_cb_new(line_list, abund_param_values, obs_specs, run_dictionary, model_path):
    chi_method = 'spec'
    line_bounds = settings.line_bounds()
    abund_dic = settings.abundance_dictionary()

    abunds = glob.glob(model_path + '/*')
    he_abunds = list(set([i.split('/')[-1].split('_')[0].strip('He') for i in abunds]))
    he_abunds.sort()
    c_abunds = n_abunds = o_abunds = list(set([i.split('/')[-1].split('_')[1].strip('CNO') for i in abunds]))
    c_abunds.sort()
    n_abunds.sort()
    o_abunds.sort()

    abund_combo = list(set([i.split('/')[-1] for i in abunds]))
    abund_combo.sort()
    calculated_lines = glob.glob(abunds[0] + '/*')
    hjds = [calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0') for calc_line in calculated_lines]
    hjds = list(set(hjds))
    hjds.sort()

    # bounds = [line_bounds[line] + np.array([-5, 5]) for line in line_list]
    bounds = [line_bounds[line] + np.array([0, 0]) for line in line_list]
    obs_wave = obs_specs[hjds[0]]['wavelength']
    inds = [i for b in bounds for i, val in  enumerate(obs_wave) if val >= b[0] and val <= b[1]]
    inds = list(set(inds))
    inds.sort()
    final_wave = obs_wave[inds]

    update_c = [i for i in line_list if i in abund_dic['C']]
    update_n = [i for i in line_list if i in abund_dic['N']]
    update_o = [i for i in line_list if i in abund_dic['O']]
    line_calc_dic = {}

    mod_dic = {}
    for abund in abund_combo:
        a_d = {}
        for hjd in hjds:
            hjd_d = {}
            for line in line_list:
                w,f = np.loadtxt(model_path + '/' + abund + '/hjd' + hjd + '_' + line + '.txt').T
                l_d = {'w':w, 'f':f}
                hjd_d[line] = l_d
            a_d[hjd] = hjd_d
        mod_dic[abund] = a_d

    chi_array = []
    for he in he_abunds:
        for line in line_list:
            line_calc_dic[line] = 'He' + str(he) + '_CNO7.5'
        for c in c_abunds:
            for line in update_c:
                line_calc_dic[line] = 'He' + str(he) + '_CNO' + str(c)
            for n in n_abunds:
                for line in update_n:
                    line_calc_dic[line] = 'He' + str(he) + '_CNO' + str(n)
                for o in o_abunds:
                    for line in update_o:
                        line_calc_dic[line] = 'He' + str(he) + '_CNO' + str(o)

                    chi2 = 0
                    for hjd in hjds:
                        exp_wave_all, exp_flux_all = [], []
                        for line in line_list:
                            exp_line_wave = mod_dic[line_calc_dic[line]][hjd][line]['w']
                            exp_line_flux = mod_dic[line_calc_dic[line]][hjd][line]['f']
                            exp_wave_all.append(exp_line_wave)
                            exp_flux_all.append(exp_line_flux)
                        ind_order = np.argsort([wave[0] for wave in exp_wave_all])
                        exp_wave = exp_wave_all[ind_order[0]]
                        exp_flux = exp_flux_all[ind_order[0]]
                        for j in range(1, len(ind_order)):
                            exp_wave, exp_flux = fw_stitch(exp_wave, exp_flux, exp_wave_all[ind_order[j]], exp_flux_all[ind_order[j]])

                        obs_wave = obs_specs[hjd]['wavelength']
                        obs_flux = obs_specs[hjd]['flux']
                        obs_flux_final = np.interp(final_wave, obs_wave, obs_flux)

                        exp_flux_final = np.interp(final_wave, exp_wave, exp_flux)
                        chi2 += calc_chi2(obs_flux_final, exp_flux_final)

                    fillout_factor = run_dictionary['fillout_factor']
                    teff_primary = run_dictionary['teff_primary']
                    teff_secondary = run_dictionary['teff_secondary']
                    period = run_dictionary['period']
                    sma = run_dictionary['sma']
                    inclination = run_dictionary['inclination']
                    q = run_dictionary['q']
                    t0 = run_dictionary['t0']
                    async_primary = run_dictionary['async_primary']
                    async_secondary = run_dictionary['async_secondary']
                    gamma = run_dictionary['gamma']
                    run_id = run_dictionary['run_id']
                    chi2_info = [chi2, fillout_factor, teff_primary, teff_secondary, period, sma, q, inclination, gamma, t0, async_primary, async_secondary, float(he), float(c), float(n), float(o), run_id]
                    chi_array.append(chi2_info)
    return chi_array


def calc_chi2_per_model_cb(abund_param_values, obs_specs, run_dictionary, model_path):
    chi_method = 'line'
    line_bounds = settings.line_bounds()
    abunds = glob.glob(model_path + '/*')
    chi_array = []
    for abund in abunds:
        abund_combo = abund.split('/')[-1]
        calculated_lines = glob.glob(abund + '/*')
        chi2 = 0
        if chi_method == 'spec':
            hjds = [calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0') for calc_line in calculated_lines]
            hjds = list(set(hjds))
            for hjd in hjds:
                exp_wave_all, exp_flux_all = [], []
                calculated_lines = glob.glob(abund + '/hjd' + hjd + '*')
                for calc_line in calculated_lines:
                    exp_line_wave, exp_line_flux = np.loadtxt(calc_line).T
                    exp_wave_all.append(exp_line_wave)
                    exp_flux_all.append(exp_line_flux)
                ind_order = np.argsort([wave[0] for wave in exp_wave_all])
                exp_wave = exp_wave_all[ind_order[0]]
                exp_flux = exp_flux_all[ind_order[0]]
                for j in range(1, len(ind_order)):
                    exp_wave, exp_flux = fw_stitch(exp_wave, exp_flux, exp_wave_all[ind_order[j]], exp_flux_all[ind_order[j]])

                obs_wave = obs_specs[hjd]['wavelength']
                obs_flux = obs_specs[hjd]['flux']

                exp_flux_final = np.interp(obs_wave, exp_wave, exp_flux)
                chi2 += calc_chi2(obs_flux, exp_flux_final)

        elif chi_method == 'line':
            for calc_line in calculated_lines:
                exp_wave, exp_flux = np.loadtxt(calc_line).T
                hjd = calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0')
                line_name = calc_line.split('_')[-1].split('.')[0]
                obs_wave = obs_specs[hjd]['wavelength']
                obs_flux = obs_specs[hjd]['flux']
                wavelength_corrected, obs_flux_corrected, exp_flux_corrected = correct_obs_exp(obs_wave, obs_flux, exp_wave, exp_flux)
                chi2 += calc_chi2(obs_flux_corrected, exp_flux_corrected, wavelength_corrected, line_bounds[line_name])
        fillout_factor = run_dictionary['fillout_factor']
        teff_primary = run_dictionary['teff_primary']
        teff_secondary = run_dictionary['teff_secondary']
        period = run_dictionary['period']
        sma = run_dictionary['sma']
        inclination = run_dictionary['inclination']
        q = run_dictionary['q']
        t0 = run_dictionary['t0']
        async_primary = run_dictionary['async_primary']
        async_secondary = run_dictionary['async_secondary']
        gamma = run_dictionary['gamma']
        run_id = run_dictionary['run_id']
        he = float(abund_combo.split('_')[0].strip('He'))
        cno = float(abund_combo.split('_')[1].strip('CNO'))
        chi2_info = [chi2, fillout_factor, teff_primary, teff_secondary, period, sma, q, inclination, gamma, t0, async_primary, async_secondary, he, cno, run_id]
        chi_array.append(chi2_info)
    return chi_array


def calc_chi2_per_model_b(abund_param_values, obs_specs, run_dictionary, model_path):
    chi_method = 'line'
    line_bounds = settings.line_bounds()
    abunds = glob.glob(model_path + '/*')
    chi_array = []
    for abund in abunds:
        abund_combo = abund.split('/')[-1]
        calculated_lines = glob.glob(abund + '/*')
        chi2 = 0
        if chi_method == 'spec':
            hjds = [calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0') for calc_line in calculated_lines]
            hjds = list(set(hjds))
            for hjd in hjds:
                exp_wave_all, exp_flux_all = [], []
                calculated_lines = glob.glob(abund + '/hjd' + hjd + '*')
                for calc_line in calculated_lines:
                    exp_line_wave, exp_line_flux = np.loadtxt(calc_line).T
                    exp_wave_all.append(exp_line_wave)
                    exp_flux_all.append(exp_line_flux)
                ind_order = np.argsort([wave[0] for wave in exp_wave_all])
                exp_wave = exp_wave_all[ind_order[0]]
                exp_flux = exp_flux_all[ind_order[0]]
                for j in range(1, len(ind_order)):
                    exp_wave, exp_flux = fw_stitch(exp_wave, exp_flux, exp_wave_all[ind_order[j]], exp_flux_all[ind_order[j]])

                obs_wave = obs_specs[hjd]['wavelength']
                obs_flux = obs_specs[hjd]['flux']

                exp_flux_final = np.interp(obs_wave, exp_wave, exp_flux)
                chi2 += calc_chi2(obs_flux, exp_flux_final)

        elif chi_method == 'line':
            for calc_line in calculated_lines:
                exp_wave, exp_flux = np.loadtxt(calc_line).T
                hjd = calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0')
                line_name = calc_line.split('_')[-1].split('.')[0]
                obs_wave = obs_specs[hjd]['wavelength']
                obs_flux = obs_specs[hjd]['flux']
                wavelength_corrected, obs_flux_corrected, exp_flux_corrected = correct_obs_exp(obs_wave, obs_flux, exp_wave, exp_flux)
                chi2 += calc_chi2(obs_flux_corrected, exp_flux_corrected, wavelength_corrected, line_bounds[line_name])
        r_equiv_primary = run_dictionary['r_equiv_primary']
        r_equiv_secondary = run_dictionary['r_equiv_secondary']
        teff_primary = run_dictionary['teff_primary']
        teff_secondary = run_dictionary['teff_secondary']
        period = run_dictionary['period']
        sma = run_dictionary['sma']
        inclination = run_dictionary['inclination']
        q = run_dictionary['q']
        t0 = run_dictionary['t0']
        async_primary = run_dictionary['async_primary']
        async_secondary = run_dictionary['async_secondary']
        pitch_primary = run_dictionary['pitch_primary']
        pitch_secondary = run_dictionary['pitch_secondary']
        yaw_primary = run_dictionary['yaw_primary']
        yaw_secondary = run_dictionary['yaw_secondary']
        gamma = run_dictionary['gamma']
        run_id = run_dictionary['run_id']
        he = float(abund_combo.split('_')[0].strip('He'))
        cno = float(abund_combo.split('_')[1].strip('CNO'))
        chi2_info = [chi2, r_equiv_primary, r_equiv_secondary, teff_primary, teff_secondary, period, sma, q, inclination, gamma, t0, async_primary, async_secondary, pitch_primary, pitch_secondary, yaw_primary, yaw_secondary, he, cno, run_id]
        chi_array.append(chi2_info)
    return chi_array


def calc_chi2_per_model_s(abund_param_values, obs_specs, run_dictionary, model_path):
    chi_method = 'line'
    line_bounds = settings.line_bounds()
    abunds = glob.glob(model_path + '/*')
    chi_array = []
    for abund in abunds:
        abund_combo = abund.split('/')[-1]
        calculated_lines = glob.glob(abund + '/*')
        chi2 = 0
        if chi_method == 'spec':
            hjds = [calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0') for calc_line in calculated_lines]
            hjds = list(set(hjds))
            for hjd in hjds:
                exp_wave_all, exp_flux_all = [], []
                calculated_lines = glob.glob(abund + '/hjd' + hjd + '*')
                for calc_line in calculated_lines:
                    exp_line_wave, exp_line_flux = np.loadtxt(calc_line).T
                    exp_wave_all.append(exp_line_wave)
                    exp_flux_all.append(exp_line_flux)
                ind_order = np.argsort([wave[0] for wave in exp_wave_all])
                exp_wave = exp_wave_all[ind_order[0]]
                exp_flux = exp_flux_all[ind_order[0]]
                for j in range(1, len(ind_order)):
                    exp_wave, exp_flux = fw_stitch(exp_wave, exp_flux, exp_wave_all[ind_order[j]], exp_flux_all[ind_order[j]])

                obs_wave = obs_specs[hjd]['wavelength']
                obs_flux = obs_specs[hjd]['flux']

                exp_flux_final = np.interp(obs_wave, exp_wave, exp_flux)
                chi2 += calc_chi2(obs_flux, (exp_flux_final-1)*1.0+1)

        elif chi_method == 'line':
            for calc_line in calculated_lines:
                exp_wave, exp_flux = np.loadtxt(calc_line).T
                hjd = calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0')
                line_name = calc_line.split('_')[-1].split('.')[0]
                obs_wave = obs_specs[hjd]['wavelength']
                obs_flux = obs_specs[hjd]['flux']
                wavelength_corrected, obs_flux_corrected, exp_flux_corrected = correct_obs_exp(obs_wave, obs_flux, exp_wave, exp_flux)
                chi2 += calc_chi2(obs_flux_corrected, exp_flux_corrected, wavelength_corrected, line_bounds[line_name])

        teff = run_dictionary['teff']
        rotation_rate = run_dictionary['rotation_rate']
        mass = run_dictionary['mass']
        r = run_dictionary['requiv']
        inclination = run_dictionary['inclination']
        t0 = run_dictionary['t0']
        gamma = run_dictionary['gamma']
        run_id = run_dictionary['run_id']
        he = float(abund_combo.split('_')[0].strip('He'))
        cno = float(abund_combo.split('_')[1].strip('CNO'))
        chi2_info = [chi2, teff, rotation_rate, mass, r, inclination, gamma, t0, he, cno, run_id]
        chi_array.append(chi2_info)
    return chi_array


def calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path):
    chi_method = 'spec'
    line_bounds = settings.line_bounds()
    abund_dic = settings.abundance_dictionary()

    abunds = glob.glob(model_path + '/He*')
    he_abunds = list(set([i.split('/')[-1].split('_')[0].strip('He') for i in abunds]))
    he_abunds.sort()
    c_abunds = n_abunds = o_abunds = list(set([i.split('/')[-1].split('_')[1].strip('CNO') for i in abunds]))
    c_abunds.sort()
    n_abunds.sort()
    o_abunds.sort()

    abund_combo = list(set([i.split('/')[-1] for i in abunds]))
    abund_combo.sort()
    calculated_lines = glob.glob(abunds[0] + '/*')
    hjds = [calc_line.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0') for calc_line in calculated_lines]
    hjds = list(set(hjds))
    hjds.sort()

    # bounds = [line_bounds[line] + np.array([-5, 5]) for line in line_list]
    bounds = [line_bounds[line] + np.array([0, 0]) for line in line_list]
    obs_wave = obs_specs[hjds[0]]['wavelength']
    inds = [i for b in bounds for i, val in  enumerate(obs_wave) if val >= b[0] and val <= b[1]]
    inds = list(set(inds))
    inds.sort()
    final_wave = obs_wave[inds]

    update_c = [i for i in line_list if i in abund_dic['C']]
    update_n = [i for i in line_list if i in abund_dic['N']]
    update_o = [i for i in line_list if i in abund_dic['O']]
    line_calc_dic = {}

    mod_dic = {}
    for abund in abund_combo:
        a_d = {}
        for hjd in hjds:
            hjd_d = {}
            for line in line_list:
                w,f = np.loadtxt(model_path + '/' + abund + '/hjd' + hjd + '_' + line + '.txt').T
                l_d = {'w':w, 'f':f}
                hjd_d[line] = l_d
            a_d[hjd] = hjd_d
        mod_dic[abund] = a_d

    chi_array = []
    for he in he_abunds:
        for line in line_list:
            line_calc_dic[line] = 'He' + str(he) + '_CNO7.5'
        for c in c_abunds:
            for line in update_c:
                line_calc_dic[line] = 'He' + str(he) + '_CNO' + str(c)
            for n in n_abunds:
                for line in update_n:
                    line_calc_dic[line] = 'He' + str(he) + '_CNO' + str(n)
                for o in o_abunds:
                    for line in update_o:
                        line_calc_dic[line] = 'He' + str(he) + '_CNO' + str(o)

                    chi2 = 0
                    for hjd in hjds:
                        exp_wave_all, exp_flux_all = [], []
                        for line in line_list:
                            exp_line_wave = mod_dic[line_calc_dic[line]][hjd][line]['w']
                            exp_line_flux = mod_dic[line_calc_dic[line]][hjd][line]['f']
                            exp_wave_all.append(exp_line_wave)
                            exp_flux_all.append(exp_line_flux)
                        ind_order = np.argsort([wave[0] for wave in exp_wave_all])
                        exp_wave = exp_wave_all[ind_order[0]]
                        exp_flux = exp_flux_all[ind_order[0]]
                        for j in range(1, len(ind_order)):
                            exp_wave, exp_flux = fw_stitch(exp_wave, exp_flux, exp_wave_all[ind_order[j]], exp_flux_all[ind_order[j]])

                        obs_wave = obs_specs[hjd]['wavelength']
                        obs_flux = obs_specs[hjd]['flux']
                        obs_flux_final = np.interp(final_wave, obs_wave, obs_flux)

                        exp_flux_final = np.interp(final_wave, exp_wave, exp_flux)
                        chi2 += calc_chi2(obs_flux_final, exp_flux_final)
                    if io_dict['object_type'] == 'contact_binary':
                        fillout_factor = run_dictionary['fillout_factor']
                        teff_primary = run_dictionary['teff_primary']
                        teff_secondary = run_dictionary['teff_secondary']
                        period = run_dictionary['period']
                        sma = run_dictionary['sma']
                        inclination = run_dictionary['inclination']
                        q = run_dictionary['q']
                        t0 = run_dictionary['t0']
                        async_primary = run_dictionary['async_primary']
                        async_secondary = run_dictionary['async_secondary']
                        gamma = run_dictionary['gamma']
                        v_macro = run_dictionary['v_macro']
                        run_id = run_dictionary['run_id']
                        chi2_info = [chi2, fillout_factor, teff_primary, teff_secondary, period, sma, q, inclination, gamma, v_macro, t0, async_primary, async_secondary, float(he), float(c), float(n), float(o), run_id, 1]

                    elif io_dict['object_type'] == 'binary':
                        r_equiv_primary = run_dictionary['r_equiv_primary']
                        r_equiv_secondary = run_dictionary['r_equiv_secondary']
                        teff_primary = run_dictionary['teff_primary']
                        teff_secondary = run_dictionary['teff_secondary']
                        period = run_dictionary['period']
                        sma = run_dictionary['sma']
                        inclination = run_dictionary['inclination']
                        q = run_dictionary['q']
                        t0 = run_dictionary['t0']
                        async_primary = run_dictionary['async_primary']
                        async_secondary = run_dictionary['async_secondary']
                        pitch_primary = run_dictionary['pitch_primary']
                        pitch_secondary = run_dictionary['pitch_secondary']
                        yaw_primary = run_dictionary['yaw_primary']
                        yaw_secondary = run_dictionary['yaw_secondary']
                        gamma = run_dictionary['gamma']
                        v_macro = run_dictionary['v_macro']
                        run_id = run_dictionary['run_id']
                        chi2_info = [chi2, r_equiv_primary, r_equiv_secondary, teff_primary, teff_secondary, period, sma, q, inclination, gamma, v_macro, t0, async_primary, async_secondary, pitch_primary, pitch_secondary, yaw_primary, yaw_secondary, float(he), float(c), float(n), float(o), run_id, 1]

                    elif io_dict['object_type'] == 'single':
                        teff = run_dictionary['teff']
                        mass = run_dictionary['mass']
                        r = run_dictionary['requiv']
                        r_pole = run_dictionary['r_pole']
                        if run_dictionary['inclination'] == -1:
                            inclination = np.arcsin(run_dictionary['vsini'] / run_dictionary['rotation_rate']) * 180/np.pi
                        else:
                            inclination = run_dictionary['inclination']
                        t0 = run_dictionary['t0']
                        gamma = run_dictionary['gamma']
                        v_macro = run_dictionary['v_macro']
                        run_id = run_dictionary['run_id']
                        v_crit_frac = run_dictionary['v_crit_frac']
                        if v_crit_frac == -1:
                            if run_dictionary['vsini'] == -1:
                                vsini = run_dictionary['rotation_rate'] * np.sin(run_dictionary['inclination'] * np.pi/180.)
                                rotation_rate = run_dictionary['rotation_rate']
                            if run_dictionary['rotation_rate'] == -1:
                                rotation_rate = run_dictionary['vsini'] / np.sin(run_dictionary['inclination'] * np.pi/180.)
                                vsini = run_dictionary['vsini']
                            else:
                                vsini = run_dictionary['vsini']
                                rotation_rate = run_dictionary['rotation_rate']
                        else:
                            rotation_rate = run_dictionary['rotation_rate']
                            vsini = run_dictionary['vsini']
                        chi2_info = [chi2, teff, vsini, rotation_rate, v_crit_frac, mass, r, r_pole, inclination, gamma, v_macro, t0, float(he), float(c), float(n), float(o), run_id, 1]

                    chi_array.append(chi2_info)
    return chi_array


def calc_chi2_per_model_TK(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path):

    calculated_regions = glob.glob(model_path + '/hjd*')
    hjds = [calc_region.split('/')[-1].split('_')[0].strip('hjd').ljust(13, '0') for calc_region in calculated_regions]
    hjds = list(set(hjds))
    hjds.sort()

    bounds = [line for line in line_list]
    obs_wave = obs_specs[hjds[0]]['wavelength']
    inds = [i for b in bounds for i, val in  enumerate(obs_wave) if val >= b[0] and val <= b[1]]
    inds = list(set(inds))
    inds.sort()
    final_wave = obs_wave[inds]

    line_calc_dic = {}

    mod_dic = {}
    for hjd in hjds:
        hjd_d = {}
        for line in line_list:
            wave_range_string = '%0.2f-%0.2f' % (line[0], line[1])
            w,f = np.loadtxt(model_path + '/hjd' + hjd + '_' + wave_range_string + '.txt').T
            l_d = {'w':w, 'f':f}
            hjd_d[wave_range_string] = l_d
        mod_dic[hjd] = hjd_d

    chi_array = []
    chi2 = 0
    for hjd in hjds:
        exp_wave_all, exp_flux_all = [], []
        for line in line_list:
            wave_range_string = '%0.2f-%0.2f' % (line[0], line[1])
            exp_line_wave = mod_dic[hjd][wave_range_string]['w']
            exp_line_flux = mod_dic[hjd][wave_range_string]['f']
            exp_wave_all.append(exp_line_wave)
            exp_flux_all.append(exp_line_flux)
        ind_order = np.argsort([wave[0] for wave in exp_wave_all])
        exp_wave = exp_wave_all[ind_order[0]]
        exp_flux = exp_flux_all[ind_order[0]]
        for j in range(1, len(ind_order)):
            exp_wave, exp_flux = fw_stitch(exp_wave, exp_flux, exp_wave_all[ind_order[j]], exp_flux_all[ind_order[j]])

        obs_wave = obs_specs[hjd]['wavelength']
        obs_flux = obs_specs[hjd]['flux']
        obs_flux_final = np.interp(final_wave, obs_wave, obs_flux)

        exp_flux_final = np.interp(final_wave, exp_wave, exp_flux)
        chi2 += calc_chi2(obs_flux_final, exp_flux_final)

    if io_dict['object_type'] == 'contact_binary':
        fillout_factor = run_dictionary['fillout_factor']
        teff_primary = run_dictionary['teff_primary']
        teff_secondary = run_dictionary['teff_secondary']
        period = run_dictionary['period']
        sma = run_dictionary['sma']
        inclination = run_dictionary['inclination']
        q = run_dictionary['q']
        t0 = run_dictionary['t0']
        async_primary = run_dictionary['async_primary']
        async_secondary = run_dictionary['async_secondary']
        gamma = run_dictionary['gamma']
        v_macro = run_dictionary['v_macro']
        run_id = run_dictionary['run_id']
        metallicity = run_dictionary['metallicity']
        if io_dict['grid_type'] == 'K':
            alpha_enhancement = run_dictionary['alpha_enhancement']
            chi2_info = [chi2, fillout_factor, teff_primary, teff_secondary, period, sma, q, inclination, gamma, v_macro, t0, async_primary, async_secondary, metallicity, alpha_enhancement, run_id, 1]
        else:
            chi2_info = [chi2, fillout_factor, teff_primary, teff_secondary, period, sma, q, inclination, gamma, v_macro, t0, async_primary, async_secondary, metallicity, run_id, 1]

    elif io_dict['object_type'] == 'binary':
        r_equiv_primary = run_dictionary['r_equiv_primary']
        r_equiv_secondary = run_dictionary['r_equiv_secondary']
        teff_primary = run_dictionary['teff_primary']
        teff_secondary = run_dictionary['teff_secondary']
        period = run_dictionary['period']
        sma = run_dictionary['sma']
        inclination = run_dictionary['inclination']
        q = run_dictionary['q']
        t0 = run_dictionary['t0']
        async_primary = run_dictionary['async_primary']
        async_secondary = run_dictionary['async_secondary']
        pitch_primary = run_dictionary['pitch_primary']
        pitch_secondary = run_dictionary['pitch_secondary']
        yaw_primary = run_dictionary['yaw_primary']
        yaw_secondary = run_dictionary['yaw_secondary']
        gamma = run_dictionary['gamma']
        v_macro = run_dictionary['v_macro']
        run_id = run_dictionary['run_id']
        metallicity = run_dictionary['metallicity']
        if io_dict['grid_type'] == 'K':
            alpha_enhancement = run_dictionary['alpha_enhancement']
            chi2_info = [chi2, r_equiv_primary, r_equiv_secondary, teff_primary, teff_secondary, period, sma, q, inclination, gamma, v_macro, t0, async_primary, async_secondary, pitch_primary, pitch_secondary, yaw_primary, yaw_secondary, metallicity, alpha_enhancement, run_id, 1]
        else:
            chi2_info = [chi2, r_equiv_primary, r_equiv_secondary, teff_primary, teff_secondary, period, sma, q, inclination, gamma, v_macro, t0, async_primary, async_secondary, pitch_primary, pitch_secondary, yaw_primary, yaw_secondary, metallicity, run_id, 1]

    elif io_dict['object_type'] == 'single':
        teff = run_dictionary['teff']
        mass = run_dictionary['mass']
        r = run_dictionary['requiv']
        r_pole = run_dictionary['r_pole']
        if run_dictionary['inclination'] == -1:
            inclination = np.arcsin(run_dictionary['vsini'] / run_dictionary['rotation_rate']) * 180/np.pi
        else:
            inclination = run_dictionary['inclination']
        t0 = run_dictionary['t0']
        gamma = run_dictionary['gamma']
        v_macro = run_dictionary['v_macro']
        run_id = run_dictionary['run_id']
        v_crit_frac = run_dictionary['v_crit_frac']
        if v_crit_frac == -1:
            if run_dictionary['vsini'] == -1:
                vsini = run_dictionary['rotation_rate'] * np.sin(run_dictionary['inclination'] * np.pi/180.)
                rotation_rate = run_dictionary['rotation_rate']
            if run_dictionary['rotation_rate'] == -1:
                rotation_rate = run_dictionary['vsini'] / np.sin(run_dictionary['inclination'] * np.pi/180.)
                vsini = run_dictionary['vsini']
            else:
                vsini = run_dictionary['vsini']
                rotation_rate = run_dictionary['rotation_rate']
        else:
            rotation_rate = run_dictionary['rotation_rate']
            vsini = run_dictionary['vsini']
        metallicity = run_dictionary['metallicity']
        if io_dict['grid_type'] == 'K':
            alpha_enhancement = run_dictionary['alpha_enhancement']
            chi2_info = [chi2, teff, vsini, rotation_rate, v_crit_frac, mass, r, r_pole, inclination, gamma, v_macro, t0, metallicity, alpha_enhancement, run_id, 1]
        else:
            chi2_info = [chi2, teff, vsini, rotation_rate, v_crit_frac, mass, r, r_pole, inclination, gamma, v_macro, t0, metallicity, run_id, 1]

    chi_array.append(chi2_info)
    return chi_array


def PFGS(times, abund_param_values, line_list, io_dict, obs_specs, run_dictionary):
    print('starting...' + str(run_dictionary['run_id']))
    model_path = update_output_directories(times, abund_param_values, io_dict, run_dictionary)
    if io_dict['object_type'] == 'contact_binary':
        if io_dict['grid_type'] in ['FW', 'FWNN']:
            chi_array_fail = chi_array = [[9999, run_dictionary['fillout_factor'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], -1, -1, -1, -1, run_dictionary['run_id'], 1]]
        elif io_dict['grid_type'] == 'K':
            chi_array_fail = chi_array = [[9999, run_dictionary['fillout_factor'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], run_dictionary['metallicity'], run_dictionary['alpha_enhancement'], run_dictionary['run_id'], 1]]
        elif io_dict['grid_type'] == 'T':
            chi_array_fail = chi_array = [[9999, run_dictionary['fillout_factor'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], run_dictionary['metallicity'], run_dictionary['run_id'], 1]]
        try:
            cb = run_cb_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
            spec_by_phase_cb(cb, line_list, abund_param_values, io_dict, run_dictionary, model_path)
            if obs_specs == None:
                chi_array = chi_array_fail
            else:
                if io_dict['grid_type'] in ['FW', 'FWNN']:
                    chi_array = calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
                elif io_dict['grid_type'] in ['K', 'T']:
                    chi_array = calc_chi2_per_model_TK(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
        except FileNotFoundError:
            print('\nFileNotFoundError: At least one patch falls outside of the specified grid.  This model will be skipped.  To prevent this in the future, run grid checks first by passing "-c" when running SPAMMS to make sure that all of the models fall within the grid.')
            chi_array = chi_array_fail

    elif io_dict['object_type'] == 'binary':
        if io_dict['grid_type'] in ['FW', 'FWNN']:
            chi_array_fail = [[9999, run_dictionary['r_equiv_primary'], run_dictionary['r_equiv_secondary'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], run_dictionary['pitch_primary'], run_dictionary['pitch_secondary'], run_dictionary['yaw_primary'], run_dictionary['yaw_secondary'], -1, -1, -1, -1, run_dictionary['run_id'], 1]]
        elif io_dict['grid_type'] == 'K':
            chi_array_fail = [[9999, run_dictionary['r_equiv_primary'], run_dictionary['r_equiv_secondary'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], run_dictionary['pitch_primary'], run_dictionary['pitch_secondary'], run_dictionary['yaw_primary'], run_dictionary['yaw_secondary'], run_dictionary['metallicity'], run_dictionary['alpha_enhancement'], run_dictionary['run_id'], 1]]
        elif io_dict['grid_type'] == 'T':
            chi_array_fail = [[9999, run_dictionary['r_equiv_primary'], run_dictionary['r_equiv_secondary'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], run_dictionary['pitch_primary'], run_dictionary['pitch_secondary'], run_dictionary['yaw_primary'], run_dictionary['yaw_secondary'], run_dictionary['metallicity'], run_dictionary['run_id'], 1]]
        try:
            b = run_b_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
            spec_by_phase_b(b, line_list, abund_param_values, io_dict, run_dictionary, model_path)
            if obs_specs == None:
                chi_array = chi_array_fail
            else:
                if io_dict['grid_type'] in ['FW', 'FWNN']:
                    chi_array = calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
                elif io_dict['grid_type'] in ['K', 'T']:
                    chi_array = calc_chi2_per_model_TK(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
        except FileNotFoundError:
            print('\nFileNotFoundError: At least one patch falls outside of the specified grid.  This model will be skipped.  To prevent this in the future, run grid checks first by passing "-c" when running SPAMMS to make sure that all of the models fall within the grid.')
            chi_array = chi_array_fail
    elif io_dict['object_type'] == 'single':
        if io_dict['grid_type'] in ['FW', 'FWNN']:
            chi_array_fail = [[9999, run_dictionary['teff'], run_dictionary['vsini'], run_dictionary['rotation_rate'], run_dictionary['v_crit_frac'], run_dictionary['mass'], run_dictionary['requiv'], run_dictionary['r_pole'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], -1, -1, -1, -1, run_dictionary['run_id'], 1]]
        elif io_dict['grid_type'] == 'K':
            chi_array_fail = [[9999, run_dictionary['teff'], run_dictionary['vsini'], run_dictionary['rotation_rate'], run_dictionary['v_crit_frac'], run_dictionary['mass'], run_dictionary['requiv'], run_dictionary['r_pole'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['metallicity'], run_dictionary['alpha_enhancement'], run_dictionary['run_id'], 1]]
        elif io_dict['grid_type'] == 'T':
            chi_array_fail = [[9999, run_dictionary['teff'], run_dictionary['vsini'], run_dictionary['rotation_rate'], run_dictionary['v_crit_frac'], run_dictionary['mass'], run_dictionary['requiv'], run_dictionary['r_pole'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['v_macro'], run_dictionary['t0'], run_dictionary['metallicity'], run_dictionary['run_id'], 1]]
        try:
            if io_dict['distortion'] in ['rotstar', 'sphere']:
                s = run_s_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
                spec_by_phase_s(s, line_list, abund_param_values, io_dict, run_dictionary, model_path)
            elif io_dict['distortion'] in ['roche']:
                s = run_sb_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
                spec_by_phase_sb(s, line_list, abund_param_values, io_dict, run_dictionary, model_path)
            if obs_specs == None:
                chi_array = chi_array_fail
            else:
                if io_dict['grid_type'] in ['FW', 'FWNN']:
                    chi_array = calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
                elif io_dict['grid_type'] in ['K', 'T']:
                    chi_array = calc_chi2_per_model_TK(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
        except FileNotFoundError:
            print('\nFileNotFoundError: At least one patch falls outside of the specified grid.  This model will be skipped.  To prevent this in the future, run grid checks first by passing "-c" when running SPAMMS to make sure that all of the models fall within the grid.')
            chi_array = chi_array_fail

    return chi_array


def main():
    phoebe.mpi.off()

    run_checks = False
    rad_bound = False
    input_file = 'input.txt'
    opts, args = getopt.getopt(sys.argv[1:], 'i:n:bc', ['input=', 'n_cores=', 'bound', 'checks'])
    for opt, arg in opts:
        if opt in ('-i', '--input'):
            input_file = str(arg)
        if opt in ('-b', '--bound'):
            rad_bound = True
        if opt in ('-n', '--n_cores'):
            n_cores = int(str(arg))
        if opt in ('-c', '--checks'):
            run_checks = True

    try:
        from schwimmbad import MultiPool
        pool = MultiPool(processes = n_cores)
    except:
        pass

    if 'pool' in locals():
        MPI = True
    else:
        MPI = False




    fit_param_values, abund_param_values, line_list, io_dict = read_input_file(input_file)
    io_dict['rad_bound'] = rad_bound
    if run_checks == False:
        setup_output_directory(io_dict)
        check_input_spectra(io_dict)
    times, obs_specs = get_obs_spec_and_times(io_dict)

    run_dictionaries = create_runs_and_ids(fit_param_values)
    # run_dictionary = run_dictionaries[0]
    # chi2 = run_phoebe_model(times, abund_param_values, io_dict, run_dictionary)

    if run_checks:
        grid = glob.glob(io_dict['path_to_grid'] + 'T*')
        grid_entries = [i.split('/')[-1] for i in grid]

        if MPI:
            combs = list(pool.map(functools.partial(check_grid, times, abund_param_values, io_dict, grid_entries), run_dictionaries))
            # pool.close()
        else:
            combs = list(map(functools.partial(check_grid, times, abund_param_values, io_dict, grid_entries), run_dictionaries))

        flat_list = [i for j in combs for i in j]
        final_list = list(set(flat_list))
        final_list.sort()

        if len(final_list) > 0:
            error_string = 'The chosen parameters result in patches that fall outside of the specified grid. The missing grid points are: \n{}'.format(final_list)
            raise ValueError('Failed to pass grid checks: \n{}'.format(error_string))
        else:
            print('Grid checks complete, no issues detected')
            sys.exit()



    print('hello')

    if MPI:
        chi2 = list(pool.map(functools.partial(PFGS, times, abund_param_values, line_list, io_dict, obs_specs), run_dictionaries))
        pool.close()
    else:
        chi2 = list(map(functools.partial(PFGS, times, abund_param_values, line_list, io_dict, obs_specs), run_dictionaries))
        # flat_list = [i for j in chi2 for i in j]
        # print(chi2)
    chi_full_array = []
    for i in chi2:
        chi_full_array.extend(i)

    # if obs_specs != None:
    #     # print len(chi_full_array)
    #
    if io_dict['object_type'] == 'contact_binary':
        if io_dict['grid_type'] in ['FW', 'FWNN']:
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.3f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.3f %0.2f %0.2f %0.2f %i %i', header = 'chi2 fillout_factor teff_primary teff_secondary period sma q inclination gamma v_macro t0 async_primary async_secondary he c n o run_id run_success')
        elif io_dict['grid_type'] == 'K':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.3f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.3f %0.2f %i %i', header = 'chi2 fillout_factor teff_primary teff_secondary period sma q inclination gamma v_macro t0 async_primary async_secondary metallicity alpha_enhancement run_id run_success')
        elif io_dict['grid_type'] == 'T':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.3f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.3f %i %i', header = 'chi2 fillout_factor teff_primary teff_secondary period sma q inclination gamma v_macro t0 async_primary async_secondary metallicity run_id run_success')
    elif io_dict['object_type'] == 'binary':
        if io_dict['grid_type'] in ['FW', 'FWNN']:
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.2f %0.2f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.1f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.2f %i %i', header = 'chi2 r_equiv_primary r_equiv_secondary teff_primary teff_secondary period sma q inclination gamma v_macro t0 async_primary async_secondary pitch_primary pitch_secondary yaw_primary yaw_secondary he c n o run_id run_success')
        elif io_dict['grid_type'] == 'K':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.2f %0.2f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.1f %0.1f %0.1f %0.1f %0.3f %0.2f %i %i', header = 'chi2 r_equiv_primary r_equiv_secondary teff_primary teff_secondary period sma q inclination gamma v_macro t0 async_primary async_secondary pitch_primary pitch_secondary yaw_primary yaw_secondary metallicity alpha_enhancement run_id run_success')
        elif io_dict['grid_type'] == 'T':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.2f %0.2f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.1f %0.3f %0.2f %0.2f %0.1f %0.1f %0.1f %0.1f %0.3f %i %i', header = 'chi2 r_equiv_primary r_equiv_secondary teff_primary teff_secondary period sma q inclination gamma v_macro t0 async_primary async_secondary pitch_primary pitch_secondary yaw_primary yaw_secondary metallicity run_id run_success')
    if io_dict['object_type'] == 'single':
        if io_dict['grid_type'] in ['FW', 'FWNN']:
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %d %0.1f %0.1f %0.4f %0.1f %0.1f %0.2f %0.1f %0.1f %0.1f %0.3f %0.3f %0.2f %0.2f %0.2f %i %i', header = 'chi2 teff vsini rotation_rate v_crit_frac mass r r_pole inclination gamma v_macro t0 he c n o run_id run_success')
        elif io_dict['grid_type'] == 'K':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %d %0.1f %0.1f %0.4f %0.1f %0.1f %0.2f %0.1f %0.1f %0.1f %0.3f %0.3f %0.2f %i %i', header = 'chi2 teff vsini rotation_rate v_crit_frac mass r r_pole inclination gamma v_macro t0 metallicity alpha_enhancement run_id run_success')
        elif io_dict['grid_type'] == 'T':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %d %0.1f %0.1f %0.4f %0.1f %0.1f %0.2f %0.1f %0.1f %0.1f %0.3f %0.3f %i %i', header = 'chi2 teff vsini rotation_rate v_crit_frac mass r r_pole inclination gamma v_macro t0 metallicity run_id run_success')

    run_successes = [i[-1]==1 for i in chi_full_array]
    if np.all(np.array(run_successes)):
        print('SPAMMS run finished. \nAll models ran successfully!')
    else:
        print('SPAMMS run finished. \nSome models failed to run.')

py_ver = sys.version_info[0]
try:
    phoebe_ver = float(phoebe.__version__[:3])
except:
    phoebe_ver = 2.3

if __name__ == "__main__":
    main()

# python spamms.py -n 2
