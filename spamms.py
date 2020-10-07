import numpy as np
import glob
import os
import shutil
from schwimmbad import MPIPool
import sys
import itertools
import time
import functools
from scipy.interpolate import splrep, splev
from scipy import stats
import phoebe
from phoebe import u,c
import math
from tqdm import tqdm, trange
import settings
from astropy.constants import R_sun, M_sun, G
import getopt




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


    fit_params = ['fillout_factor', 'teff_primary', 'teff_secondary', 'period', 'sma', 'inclination', 'q', 't0', 'async_primary', 'async_secondary', 'gamma']
    abundance_params = ['he_abundances', 'cno_abundances']

    fit_param_values = {}
    abund_param_values = {}
    io_dict = {'object_type':object_type, 'path_to_obs_spectra':path_to_obs_spectra, 'output_directory':output_directory, 'path_to_grid':path_to_grid, 'input_file':input_file}
    try:
        times_ind = [i for i in range(len(lines)) if lines[i].startswith('times')][0]
        times = lines[times_ind].split('=')[1].strip()
        io_dict['times'] = arg_parse(times)
    except:
        pass

    for param in fit_params:
        arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
        fit_param_values[param] = arg_parse(arg)

    abund_param_values['he_abundances'] = [0.06, 0.1, 0.15, 0.2]
    abund_param_values['cno_abundances'] = [6.5, 7.0, 7.5, 8.0, 8.5]
    abund_param_values['lp_bins'] = 161
    for param in abundance_params:
        arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
        abund_param_values[param] = arg_parse(arg)

    if type(abund_param_values['lp_bins']) is list:
        abund_param_values['lp_bins'] = int(abund_param_values['lp_bins'][0])

    abund_param_values['interpolate_abundances'] = False
    # interp_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('interpolate_abundances')][0]].split('=')[1].strip()
    # if interp_arg == 'False' or interp_arg == '0':
    #     abund_param_values['interpolate_abundances'] = False
    # else:
    #     abund_param_values['interpolate_abundances'] = True

    line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_line_list =')][0]].split('=')[1].strip()
    line_list = parse_line_list(line_list_arg)

    return fit_param_values, abund_param_values, line_list, io_dict


def read_b_input_file(input_file):
    lines = tuple(open(input_file, 'r'))

    object_type_ind = [i for i in range(len(lines)) if lines[i].startswith('object_type')][0]
    object_type = lines[object_type_ind].split('=')[1].strip()

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


    fit_params = ['r_equiv_primary', 'r_equiv_secondary', 'teff_primary', 'teff_secondary', 'period', 'sma', 'inclination', 'q', 't0', 'async_primary', 'async_secondary', 'pitch_primary', 'pitch_secondary', 'yaw_primary', 'yaw_secondary', 'gamma']
    abundance_params = ['he_abundances', 'cno_abundances']

    fit_param_values = {}
    abund_param_values = {}
    io_dict = {'object_type':object_type, 'path_to_obs_spectra':path_to_obs_spectra, 'output_directory':output_directory, 'path_to_grid':path_to_grid, 'input_file':input_file}
    try:
        times_ind = [i for i in range(len(lines)) if lines[i].startswith('times')][0]
        times = lines[times_ind].split('=')[1].strip()
        io_dict['times'] = arg_parse(times)
    except:
        pass

    for param in fit_params:
        arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
        fit_param_values[param] = arg_parse(arg)

    abund_param_values['he_abundances'] = [0.06, 0.1, 0.15, 0.2]
    abund_param_values['cno_abundances'] = [6.5, 7.0, 7.5, 8.0, 8.5]
    abund_param_values['lp_bins'] = 161
    for param in abundance_params:
        arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
        abund_param_values[param] = arg_parse(arg)

    if type(abund_param_values['lp_bins']) is list:
        abund_param_values['lp_bins'] = int(abund_param_values['lp_bins'][0])

    abund_param_values['interpolate_abundances'] = False
    # interp_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('interpolate_abundances')][0]].split('=')[1].strip()
    # if interp_arg == 'False' or interp_arg == '0':
    #     abund_param_values['interpolate_abundances'] = False
    # else:
    #     abund_param_values['interpolate_abundances'] = True

    line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_line_list =')][0]].split('=')[1].strip()
    line_list = parse_line_list(line_list_arg)

    return fit_param_values, abund_param_values, line_list, io_dict


def read_s_input_file(input_file):
    lines = tuple(open(input_file, 'r'))

    object_type_ind = [i for i in range(len(lines)) if lines[i].startswith('object_type')][0]
    object_type = lines[object_type_ind].split('=')[1].strip()

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


    fit_params = ['teff', 'rotation_rate', 'requiv', 'inclination', 'mass', 't0', 'gamma']
    fit_params_alt = ['teff', 'vsini', 'rotation_rate', 'requiv', 'inclination', 'mass', 't0', 'gamma']
    abundance_params = ['he_abundances', 'cno_abundances']

    fit_param_values = {}
    abund_param_values = {}
    io_dict = {'object_type':object_type, 'path_to_obs_spectra':path_to_obs_spectra, 'output_directory':output_directory, 'path_to_grid':path_to_grid, 'input_file':input_file}
    try:
        times_ind = [i for i in range(len(lines)) if lines[i].startswith('times')][0]
        times = lines[times_ind].split('=')[1].strip()
        io_dict['times'] = arg_parse(times)
    except:
        pass

    for param in fit_params_alt:
        try:
            arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
            fit_param_values[param] = arg_parse(arg)
        except:
            if param == 'vsini':
                fit_param_values[param] = [-1.0]
    for param in abundance_params:
        arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
        abund_param_values[param] = arg_parse(arg)

    abund_param_values['he_abundances'] = [0.06, 0.1, 0.15, 0.2]
    abund_param_values['cno_abundances'] = [6.5, 7.0, 7.5, 8.0, 8.5]
    abund_param_values['lp_bins'] = 161
    for param in abundance_params:
        arg = lines[[i for i in range(len(lines)) if lines[i].startswith(param)][0]].split('=')[1].strip()
        abund_param_values[param] = arg_parse(arg)

    if type(abund_param_values['lp_bins']) is list:
        abund_param_values['lp_bins'] = int(abund_param_values['lp_bins'][0])

    abund_param_values['interpolate_abundances'] = False
    # interp_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('interpolate_abundances')][0]].split('=')[1].strip()
    # if interp_arg == 'False' or interp_arg == '0':
    #     abund_param_values['interpolate_abundances'] = False
    # else:
    #     abund_param_values['interpolate_abundances'] = True

    line_list_arg = lines[[i for i in range(len(lines)) if lines[i].startswith('selected_line_list =')][0]].split('=')[1].strip()
    line_list = parse_line_list(line_list_arg)

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

    print('Checks Complete')


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
                obs_specs[str(times_temp[i]).ljust(13, '0')] = dic
        return np.array(times), obs_specs



def create_runs_and_ids(fit_param_values):
    keys = []
    values = []
    for key, value in fit_param_values.iteritems():
        keys.append(key)
        values.append(value)

    runs = list(itertools.product(*values))
    run_dictionaries = [dict(zip(keys,i)) for i in runs]
    for i, j in enumerate(run_dictionaries):
        j['run_id'] = i

    run_ids = range(len(run_dictionaries))

    dictionary_of_run_dicts = dict(zip(run_ids, run_dictionaries))

    return run_dictionaries


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
        cb['ntriangles'].set_value_all(value = 5000)
    else:
        cb['ntriangles'].set_value(value = 5000)
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

    cb['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
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
    if phoebe_ver < 2.2:
        b['ntriangles'].set_value_all(value = 5000)
    else:
        b['ntriangles'].set_value(value = 5000)
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
    b['ld_func'].set_value_all(value = 'logarithmic')
    b['atm'].set_value_all(value='blackbody')
    b['include_times'] = 't0_ref@binary'
    b.flip_constraint('t0_ref@binary', 't0_supconj')
    b['t0_ref@binary@component'].set_value(value = run_dictionary['t0'])

    b['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
    b.run_compute()

    execution_time = time.time() - start_time_prog_1
    # print execution_time
    return b


def run_s_phoebe_model(times, abund_param_values, io_dict, run_dictionary):
    start_time_prog_1 = time.time()
    logger = phoebe.logger(clevel='ERROR')

    s = phoebe.default_star()
    # s['distortion_method'].set_value(value='sphere')
    s['teff@component'].set_value(value = run_dictionary['teff'])
    s['gravb_bol'].set_value(value = 1.0)
    s['irrad_frac_refl_bol'].set_value(value = 1.0)
    s['mass@component'].set_value(value = run_dictionary['mass'])
    s['requiv@component'].set_value(value = run_dictionary['requiv'])
    if run_dictionary['rotation_rate'] == 0:
        s['distortion_method'].set_value('sphere')
    if run_dictionary['rotation_rate'] == -1:
        period = rotation_rate_to_period(run_dictionary['vsini'] / (np.sin(run_dictionary['inclination'] * np.pi/180.)), run_dictionary['requiv'])
    else:
        period = rotation_rate_to_period(run_dictionary['rotation_rate'], run_dictionary['requiv'])
    s['period@component'].set_value(value = period)
    if run_dictionary['inclination'] == -1:
        s['incl@component'].set_value(value = np.arcsin(run_dictionary['vsini'] / run_dictionary['rotation_rate']) * 180./np.pi)
    else:
        s['incl@component'].set_value(value = run_dictionary['inclination'])

    s['ntriangles'].set_value(value = 5000)

    t = list(times)

    s.add_dataset('lc', times=t, dataset='lc01')
    s.add_dataset('rv', times=t, dataset='rv01')
    s.add_dataset('mesh', times=t, dataset='mesh01')
    s.add_dataset('orb', times=t, dataset='orb01')
    s['ld_func'].set_value_all(value = 'logarithmic')
    s['atm'].set_value(value='blackbody')
    s['include_times'] = 't0@system'
    s['t0@system'].set_value(value = run_dictionary['t0'])

    s['columns'] = ['*@lc01', '*@rv01', 'us', 'vs', 'ws', 'vus', 'vvs', 'vws', 'loggs', 'teffs', 'mus', 'visibilities', 'rs', 'areas']
    s.run_compute()

    execution_time = time.time() - start_time_prog_1
    # print execution_time
    return s


def update_output_directories(times, abund_param_values, io_dict, run_dictionary):
    model_path = io_dict['output_directory'] + 'Model_' + str(run_dictionary['run_id']).zfill(4)
    os.mkdir(model_path)
    if abund_param_values['interpolate_abundances']:
        print('abundance interpolation is not supported yet.')
    he_abundances = [i for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
    cno_abundances = [j for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
    # he_abundances = [0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2]
    # cno_abundances = [6.5, 6.5, 6.5, 6.5, 7.0, 7.0, 7.0, 7.0, 7.5, 7.5, 7.5, 7.5, 8.0, 8.0, 8.0, 8.0, 8.5, 8.5, 8.5, 8.5]
    for i in range(len(he_abundances)):
        os.mkdir(model_path + '/He' + str(he_abundances[i]) + '_CNO' + str(cno_abundances[i]))
    return model_path


def calc_spec_by_phase(mesh_vals, hjd, model_path, lines, abund_param_values, lines_dic, io_dict):
    nw = []
    nf = []

    # print 'assigning spectra'
    for line in lines:
        assign_and_calc_abundance(mesh_vals, hjd, model_path, abund_param_values, lines_dic, io_dict, line)


def assign_and_calc_abundance(mesh_vals, hjd, model_path, abund_param_values, lines_dic, io_dict, line):
    start_time = time.time()
    he_abundances = [i for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
    cno_abundances = [j for j in abund_param_values['cno_abundances'] for i in abund_param_values['he_abundances']]
    # he_abundances = [0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2, 0.06, 0.1, 0.15, 0.2]
    # cno_abundances = [6.5, 6.5, 6.5, 6.5, 7.0, 7.0, 7.0, 7.0, 7.5, 7.5, 7.5, 7.5, 8.0, 8.0, 8.0, 8.0, 8.5, 8.5, 8.5, 8.5]
    ws, star_profs, wind_profs = assign_spectra_interp(mesh_vals, line, lines_dic, io_dict, abund_param_values)
    # if 'upper' in lines_dic.keys():
    #     ws, star_profs, wind_profs = assign_spectra_interp(mesh_vals, line, lines_dic)
    # else:
    #     ws, star_profs, wind_profs = assign_spectra(mesh_vals, line, lines_dic)
    waves = []
    phots = []
    lp_bins = abund_param_values['lp_bins']
    for i in range(len(ws[0])/lp_bins):
        wavg_single, phot_avg_single = calc_flux_optimize(ws[:,lp_bins*i:lp_bins*(i+1)], ws, star_profs[:,lp_bins*i:lp_bins*(i+1)], wind_profs[:,lp_bins*i:lp_bins*(i+1)], mesh_vals)
        waves.append(wavg_single)
        phots.append(phot_avg_single/phot_avg_single[-1])
        np.savetxt(model_path + '/He' + str(he_abundances[i]) + '_CNO' + str(cno_abundances[i]) + '/hjd' + str(hjd).ljust(13, '0') + '_' + line + '.txt', np.array([wavg_single, phot_avg_single/phot_avg_single[-1]]).T)

    # phot_avg = interp_abundances(phots, He_abundance, CNO_abundance)
    # wavg = waves[0]
    # return np.array(wavg), np.array(phot_avg)/phot_avg[0]
    # wave, phots = calc_flux_bulk(ws, star_profs, wind_profs, mesh_vals)
    # for i in range(20):
    #     np.savetxt(model_path + '/He' + str(he_abundances[i]) + '_CNO' + str(cno_abundances[i]) + '/hjd' + str(hjd).ljust(13, '0') + '_' + line + '.txt', np.array([wave, phots[i]/phots[i][0]]).T)
    print(time.time() - start_time)


def assign_spectra(mesh_vals, line, lines_dic, io_dict):
    ts = np.around(mesh_vals['teffs'] / 1000.0) * 1000.0
    lgs = np.around(mesh_vals['loggs']*10.) / 10.
    rads = np.around(mesh_vals['rs'] * 4.0) / 4.0
    if io_dict['rad_bound']:
        rads = rads * (rads <= 9.) + 9.*(rads > 9.)
        rads = rads * (rads >= 6.5) + 6.5*(rads < 6.5)

    # lgs = lgs * (lgs >= 3.) + 3.*(lgs < 3.)


    ws = []
    star_profs = []
    wind_profs = []
    start_time = time.time()
    for i in tqdm(range(len(ts))):
        w, st, wi = lookup_line_profs_from_dic(ts[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic)
        ws.append(w)
        star_profs.append(st)
        wind_profs.append(np.array(wi))
    elapsed_time = time.time() - start_time
    # print 'Average iterations per second: ' + str(len(ts) / elapsed_time)
    ws = dopler_shift(np.array(ws), np.array([mesh_vals['rvs']]*len(ws[0])).T)
    ws=np.array(ws, dtype='float')
    return np.array(ws), np.array(star_profs), np.array(wind_profs)


def assign_spectra_interp(mesh_vals, line, lines_dic, io_dict, abund_param_values):
    ts = mesh_vals['teffs']
    tls = np.floor(mesh_vals['teffs'] / 1000.0) * 1000.0
    tus = np.ceil(mesh_vals['teffs'] / 1000.0) * 1000.0
    w1s = (tus - ts)/1000.0
    w2s = (ts - tls)/1000.0
    lgs = np.around(mesh_vals['loggs']*10.) / 10.
    rads = np.around(mesh_vals['rs'] * 4.0) / 4.0
    if io_dict['rad_bound']:
        rads = rads * (rads <= 9.) + 9.*(rads > 9.)
        rads = rads * (rads >= 6.5) + 6.5*(rads < 6.5)
    # lgs = lgs * (lgs >= 3.) + 3.*(lgs < 3.)

    ws = []
    star_low_profs = []
    wind_low_profs = []
    star_high_profs = []
    wind_high_profs = []
    start_time = time.time()
    for i in tqdm(range(len(ts))):
        # w, stl, wil = lookup_line_profs_from_dic(tls[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic['lower'])
        # wu, stu, wiu = lookup_line_profs_from_dic(tus[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic['upper'])
        w, stl, wil = lookup_line_profs_from_dic(tls[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic)
        wu, stu, wiu = lookup_line_profs_from_dic(tus[i], lgs[i], rads[i], mesh_vals['mus'][i], mesh_vals['viss'][i], line, lines_dic)
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
    elapsed_time = time.time() - start_time
    # print 'Average iterations per second: ' + str(len(ts) / elapsed_time)
    # print mesh_vals['rvs']
    ws = dopler_shift(np.array(ws), np.array([mesh_vals['rvs']]*len(ws[0])).T)
    ws=np.array(ws, dtype='float')
    return np.array(ws), np.array(star_profs), np.array(wind_profs)


def dopler_shift(w, rv):
    c = 299792.458
    return w*c/(c-rv)


def lookup_line_profs_from_dic(t, g, r, m, v, line, lines_dic):
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


def calc_flux_optimize(ws, ws_all, star_profs, wind_profs, mesh_vals):
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

    wave = np.arange(w_min, w_max, 0.1)
    I_star = []
    for i in range(len(ws)):
        I_star.append(np.interp(wave, ws[i], star_profs[i]))

    I_star = np.array(I_star)
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

    wave = np.arange(w_min, w_max, 0.1)
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


def spec_by_phase_cb(cb, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = cb['times@dataset@lc'].value
    interp = True

    combs, mode_combs = determine_tgr_combinations(cb, io_dict)
    grid = glob.glob(io_dict['path_to_grid'] + 'T*')
    grid_entries = [i.split('/')[-1] for i in grid]
    missing_combs = [i for i in combs if i not in grid_entries]
    if len(missing_combs) > 0:
        print('WARNING: input grid entries missing.')
        print(missing_combs)
    #print combs
    lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)
    # lines_dic = line_dictionary_structure(combs, line_list, io_dict)
    # if interp:
    #     lines_dic = interp_line_dictionary_structure(lines_dic)

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

        mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict)
        # print time.time() - start_time


def spec_by_phase_b(b, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = b['times@dataset@lc'].value
    interp = True

    combs, mode_combs = determine_tgr_combinations(b, io_dict)
    lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)
    # lines_dic = line_dictionary_structure(combs, line_list, io_dict)
    # if interp:
    #     lines_dic = interp_line_dictionary_structure(lines_dic)

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

        mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict)
        # print time.time() - start_time


def spec_by_phase_s(s, line_list, abund_param_values, io_dict, run_dictionary, model_path):
    times = s['times@dataset@lc'].value

    combs, mode_combs = determine_tgr_combinations(s, io_dict)
    #print(combs)
    lines_dic = interp_line_dictionary_structure_new(combs, line_list, io_dict, mode_combs, abund_param_values)

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

        mesh_vals = {'teffs':teffs, 'loggs':loggs, 'rs':rs, 'mus':mus, 'rvs':rvs, 'viss':viss, 'abs_intens':abs_intens, 'areas':areas, 'ldints':ldints, 'rs_sol':rs_sol}

        calc_spec_by_phase(mesh_vals, hjd, model_path, line_list, abund_param_values, lines_dic, io_dict)
        # print time.time() - start_time



def determine_tgr_combinations(cb, io_dict):
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

    ts = np.around(np.array(teffs) / 1000.0) * 1000.0
    tls = np.floor(np.array(teffs) / 1000.0) * 1000.0
    tus = np.ceil(np.array(teffs) / 1000.0) * 1000.0
    lgs = np.around(np.array(loggs)*10.) / 10.
    rads = np.around(np.array(rs) * 4.0) / 4.0
    if io_dict['rad_bound']:
        rads = rads * (rads <= 9.) + 9.*(rads > 9.)
        rads = rads * (rads >= 6.5) + 6.5*(rads < 6.5)
    # lgs = lgs * (lgs >= 3.) + 3.*(lgs < 3.)
    # combinations = ['T' + str(int(ts[i])) + '_G' + str(lgs[i]) + '_R' + format(rads[i], '.2f') for i in range(len(ts))]
    combinations = ['T' + str(int(tls[i])) + '_G' + str(lgs[i]) + '_R' + format(rads[i], '.2f') for i in range(len(ts))]
    combinations.extend(['T' + str(int(tus[i])) + '_G' + str(lgs[i]) + '_R' + format(rads[i], '.2f') for i in range(len(ts))])
    return list(set(combinations)), stats.mode(np.array(combinations))[0][0]


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
                        run_id = run_dictionary['run_id']
                        chi2_info = [chi2, fillout_factor, teff_primary, teff_secondary, period, sma, q, inclination, gamma, t0, async_primary, async_secondary, float(he), float(c), float(n), float(o), run_id]

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
                        run_id = run_dictionary['run_id']
                        chi2_info = [chi2, r_equiv_primary, r_equiv_secondary, teff_primary, teff_secondary, period, sma, q, inclination, gamma, t0, async_primary, async_secondary, pitch_primary, pitch_secondary, yaw_primary, yaw_secondary, float(he), float(c), float(n), float(o), run_id]

                    elif io_dict['object_type'] == 'single':
                        teff = run_dictionary['teff']
                        if run_dictionary['rotation_rate'] == -1:
                            rotation_rate = run_dictionary['vsini'] / np.sin(run_dictionary['inclination'] * np.pi/180.)
                        else:
                            rotation_rate = run_dictionary['rotation_rate']
                        mass = run_dictionary['mass']
                        r = run_dictionary['requiv']
                        if run_dictionary['inclination'] == -1:
                            inclination = np.arcsin(run_dictionary['vsini'] / run_dictionary['rotation_rate']) * 180/np.pi
                        else:
                            inclination = run_dictionary['inclination']
                        t0 = run_dictionary['t0']
                        gamma = run_dictionary['gamma']
                        run_id = run_dictionary['run_id']
                        if run_dictionary['vsini'] == -1:
                            vsini = run_dictionary['rotation_rate'] * np.sin(run_dictionary['inclination'] * np.pi/180.)
                        else:
                            vsini = run_dictionary['vsini']
                        chi2_info = [chi2, teff, vsini, rotation_rate, mass, r, inclination, gamma, t0, float(he), float(c), float(n), float(o), run_id]

                    chi_array.append(chi2_info)
    return chi_array


def PFGS(times, abund_param_values, line_list, io_dict, obs_specs, run_dictionary):
    print('starting...' + str(run_dictionary['run_id']))
    model_path = update_output_directories(times, abund_param_values, io_dict, run_dictionary)
    if io_dict['object_type'] == 'contact_binary':
        cb = run_cb_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
        spec_by_phase_cb(cb, line_list, abund_param_values, io_dict, run_dictionary, model_path)
        try:
            if obs_specs == None:
                chi_array = [0]
            else:
                chi_array = calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
        except:
            chi_array = [[9999, run_dictionary['fillout_factor'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], -1, -1, -1, -1, run_dictionary['run_id']]]
    elif io_dict['object_type'] == 'binary':
        try:
            b = run_b_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
            spec_by_phase_b(b, line_list, abund_param_values, io_dict, run_dictionary, model_path)
            if obs_specs == None:
                chi_array = [0]
            else:
                chi_array = calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
        except:
            chi_array = [[9999, run_dictionary['r_equiv_primary'], run_dictionary['r_equiv_secondary'], run_dictionary['teff_primary'], run_dictionary['teff_secondary'], run_dictionary['period'], run_dictionary['sma'], run_dictionary['q'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['t0'], run_dictionary['async_primary'], run_dictionary['async_secondary'], run_dictionary['pitch_primary'], run_dictionary['pitch_secondary'], run_dictionary['yaw_primary'], run_dictionary['yaw_secondary'], -1, -1, -1, -1, run_dictionary['run_id']]]
    elif io_dict['object_type'] == 'single':
        try:
            s = run_s_phoebe_model(times, abund_param_values, io_dict, run_dictionary)
            spec_by_phase_s(s, line_list, abund_param_values, io_dict, run_dictionary, model_path)
            if obs_specs == None:
                chi_array = [0]
            else:
                chi_array = calc_chi2_per_model_new(line_list, abund_param_values, obs_specs, run_dictionary, io_dict, model_path)
        except:
            chi_array = [[9999, run_dictionary['teff'], run_dictionary['vsini'], run_dictionary['rotation_rate'], run_dictionary['mass'], run_dictionary['requiv'], run_dictionary['inclination'], run_dictionary['gamma'], run_dictionary['t0'], -1, -1, -1, -1, run_dictionary['run_id']]]

    return chi_array


def main():
    phoebe.mpi.off()

    try:
        pool = MPIPool()
    except:
        pass

    if 'pool' in locals():
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        MPI = True
    else:
        MPI = False

    input_file = 'input.txt'
    rad_bound = False
    opts, args = getopt.getopt(sys.argv[1:], 'i:b', ['input=', 'bound'])
    for opt, arg in opts:
        if opt in ('-i', '--input'):
            input_file = str(arg)
        if opt in ('-b', '--bound'):
            rad_bound = True



    fit_param_values, abund_param_values, line_list, io_dict = read_input_file(input_file)
    setup_output_directory(io_dict)
    check_input_spectra(io_dict)
    times, obs_specs = get_obs_spec_and_times(io_dict)
    io_dict['rad_bound'] = rad_bound

    run_dictionaries = create_runs_and_ids(fit_param_values)
    # run_dictionary = run_dictionaries[0]
    # chi2 = run_phoebe_model(times, abund_param_values, io_dict, run_dictionary)


    print('hello')

    if MPI:
        chi2 = pool.map(functools.partial(PFGS, times, abund_param_values, line_list, io_dict, obs_specs), run_dictionaries)
        pool.close()
    else:
        chi2 = map(functools.partial(PFGS, times, abund_param_values, line_list, io_dict, obs_specs), run_dictionaries)

    if obs_specs != None:
        chi_full_array = []
        for i in chi2:
            chi_full_array.extend(i)
        # print len(chi_full_array)

        if io_dict['object_type'] == 'contact_binary':
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.3f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.3f %0.2f %0.2f %0.3f %0.2f %0.2f %0.2f %s', header = 'chi2 fillout_factor teff_primary teff_secondary period sma q inclination gamma t0 async_primary async_secondary he c n o run_id')
        if io_dict['object_type'] == 'binary':
            [chi2, r_equiv_primary, r_equiv_secondary, teff_primary, teff_secondary, period, sma, q, inclination, gamma, t0, async_primary, async_secondary, pitch_primary, pitch_secondary, yaw_primary, yaw_secondary, he, cno, run_id]
            np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %0.2f %0.2f %d %d %f %0.2f %0.2f %0.1f %0.1f %0.3f %0.2f %0.2f %0.1f %0.1f %0.1f %0.1f %0.2f %0.2f %s', header = 'chi2 r_equiv_primary r_equiv_secondary teff_primary teff_secondary period sma q inclination gamma t0 async_primary async_secondary pitch_primary pitch_secondary yaw_primary yaw_secondary he cno run_id')
        elif io_dict['object_type'] == 'single':
            try:
                np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %d %0.1f %0.1f %0.1f %0.2f %0.1f %0.1f %0.3f %0.3f %0.2f %0.2f %0.2f %s', header = 'chi2 teff vsini rotation_rate mass r inclination gamma t0 he c n o run_id')
            except:
                np.savetxt(io_dict['output_directory'] + 'chi_square_summary.txt', np.array(chi_full_array), fmt='%f %d %0.1f %0.1f %0.2f %0.1f %0.1f %0.3f %0.3f %0.2f %0.2f %0.2f %s', header = 'chi2 teff rotation_rate mass r inclination gamma t0 he c n o run_id')

py_ver = sys.version_info[0]
phoebe_ver = float(phoebe.__version__[:3])

if __name__ == "__main__":
    main()

# mpiexec -n 2 python PFGS.py
