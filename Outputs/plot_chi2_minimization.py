import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scInterp
import sys
import glob
from scipy.stats import stats
sys.path.append('..')
import settings
from spamms import read_input_file

def update_annot(ind):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))), " ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
    annot.get_bbox_patch().set_alpha(0.4)

def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            update_annot(ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()


# Provide the information needed to get the data
run_location = sys.argv[1]

filename = run_location+ '/chi_square_summary.txt'

# Read in the data
# chi2, pot, Ta, Tb, period, sma, q, t0, [He], [CNO]
data = np.genfromtxt(filename, dtype=None, delimiter=' ', names=True)
ind = data['chi2'] < 9000
data = data[ind]

# params = ['pot', 'Ta', 'Tb', 'period', 'sma', 'q', 't0', 'he', 'cno']
#
# fitted_params = {}
# for i in range(len(params)):
#     param = np.unique(data[:, i + 1])
#     if len(param) > 1:
#         fitted_params[params[i]] = np.unique(data[:, i + 1])

columns = list(data.dtype.names)
columns.remove('chi2')
columns.remove('run_id')

fit_params = []
fit_params_unique = {}
for param in columns:
    if len(np.unique(data[param])) > 1:
        fit_params.append(param)
        fit_params_unique[param] = np.unique(data[param])


temp1, temp2, line_list, temp3 = read_input_file(run_location + '/input.txt')
spec_file = glob.glob(run_location + '/input_spectra/*spec.txt')
x = np.loadtxt(spec_file[0]).T

line_bounds = settings.line_bounds()
bounds = [line_bounds[line] + np.array([-5, 5]) for line in line_list]
obs_wave = x[0]
inds = [i for b in bounds for i, val in  enumerate(obs_wave) if val >= b[0] and val <= b[1]]
inds = list(set(inds))
inds.sort()
final_wave = obs_wave[inds]
n_dof = len(final_wave)*len(x[1:]) - len(fit_params)


fit_params_chi = {}
fit_param_interp = {}
fit_param_range = {}
fit_param_min = {}
fit_param_sig = {}
# fit_param_sig_ga = {}
min_chi = 9999999999999
for param in fit_params:
    print(param)
    fit_params_chi[param] = [np.min(data['chi2'][data[param] == unique]) for unique in fit_params_unique[param]]
    try:
        fit_param_interp[param] = scInterp.interp1d(fit_params_unique[param], fit_params_chi[param], kind='cubic')
    except:
        fit_param_interp[param] = scInterp.interp1d(fit_params_unique[param], fit_params_chi[param], kind='quadratic')
    fit_param_range[param] = np.linspace(np.min(fit_params_unique[param]), np.max(fit_params_unique[param]), 500)

    temp_pr = np.linspace(np.min(fit_params_unique[param]), np.max(fit_params_unique[param]), 2000)
    temp_chi = fit_param_interp[param](temp_pr)
    ind_min = np.argmin(temp_chi)
    fit_param_min[param] = temp_pr[ind_min]
    min_chi = min(min_chi, min(temp_chi))

    chi_sig = temp_chi[ind_min]*(1. + np.sqrt(2./n_dof))
    fit_param_sig[param] = chi_sig
    sig_range_ind = (temp_chi < chi_sig)
    sig_range = [min(temp_pr[sig_range_ind]), max(temp_pr[sig_range_ind])]
    print(temp_pr[ind_min])
    print(sig_range)
    # # print chi_sig
    # Ps = [stats.distributions.chi2.sf(i*n_dof/temp_chi[ind_min], n_dof) for i in temp_chi]
    # inds = (np.array(Ps) >= 0.05)
    # chi_sig2 = max(temp_chi[inds])
    # fit_param_sig_ga[param] = chi_sig2
    # print chi_sig2

chi_sig = min_chi*(1. + np.sqrt(2./n_dof))


i = np.argmin(data['chi2'])
try:
    min_mod = 'Model' + str(int(data['run_id'][i])).zfill(4) + ' He' + str(float(data['he'][i])) + '_C' + str(float(data['c'][i])) + '_N' + str(float(data['n'][i])) + '_O' + str(float(data['o'][i]))
except:
    min_mod = 'Model' + str(int(data['run_id'][i])).zfill(4) + ' He' + str(float(data['he'][i])) + '_CNO' + str(float(data['cno'][i]))
print(min_mod)


'''
for param in fit_params:
    try:
        labels = ['Model' + str(int(data['run_id'][i])).zfill(4) + ' He' + str(float(data['he'][i])) + ' C' + str(float(data['c'][i])) + ' N' + str(float(data['n'][i])) + ' O' + str(float(data['o'][i])) for i in range(len(data['run_id']))]
    except:
        labels = ['Model' + str(int(data['run_id'][i])).zfill(4) + ' He' + str(float(data['he'][i])) + ' CNO' + str(float(data['cno'][i])) for i in range(len(data['run_id']))]

    # print labels
    fig, ax = plt.subplots()
    scatter = ax.scatter(data[param], data['chi2'], c = 'black', alpha=0.5, s = 4)
    ax.set_ylabel('Chi square',  fontsize=12)
    ax.set_xlabel(param,  fontsize=12)
    ax.plot(fit_param_range[param], fit_param_interp[param](fit_param_range[param]), 'k:', alpha=0.8, lw=1)
    ax.plot(fit_param_range[param], np.ones_like(fit_param_range[param])*fit_param_sig[param], 'r')
    # ax.plot(fit_param_range[param], np.ones_like(fit_param_range[param])*fit_param_sig_ga[param], c='purple')

    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = scatter.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}".format(" ".join([labels[n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = scatter.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()
    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.show()
'''

cols = np.ceil(np.sqrt(len(fit_params)))
rows = np.ceil(len(fit_params) / cols)
fig, axs = plt.subplots(int(cols), int(rows), sharey=True, figsize = (8, 6))
for ind, param in enumerate(fit_params):
    c = int(ind%cols)
    r = int(np.floor(ind/cols))
    axs[r][c].scatter(data[param], data['chi2'], c = 'black', alpha=0.5, s = 4)
    if r == 0:
        axs[c][r].set_ylabel('Chi square')
    axs[r][c].set_xlabel(param)
    axs[r][c].plot(fit_param_range[param], fit_param_interp[param](fit_param_range[param]), 'k:', alpha=0.8, lw=1)
    axs[r][c].plot(fit_param_range[param], np.ones_like(fit_param_range[param])*fit_param_sig[param], 'r')

for ind in range(len(fit_params), int(rows*cols)):
    c = int(ind%cols)
    r = int(np.floor(ind/cols))
    axs[r][c].axis('off')

plt.subplots_adjust(top=0.935, bottom=0.095, left=0.11, right=0.9, hspace=0.22, wspace=0.2)
plt.show()





chi_sig = min_chi*(1. + np.sqrt(2./n_dof))
inds = data['chi2'] < chi_sig

mods = data[inds]

i = np.argmin(data['chi2'])
min_mod = data[i]



f_specs = glob.glob(run_location + '/input_spectra/*_spec.txt')
f_hjds = glob.glob(run_location + '/input_spectra/*_hjd.txt')
x = np.loadtxt(f_specs[0]).T
hjds = np.loadtxt(f_hjds[0], ndmin=1)

temp1, temp2, line_list, io_dict = read_input_file(run_location + '/input.txt')
try:
    times = io_dict['times']
except:
    times = hjds

line_bounds = settings.line_bounds()

fig, axs = plt.subplots(1, len(line_list))

hjd = times[0]
for ind, line in enumerate(line_list):
    y = np.loadtxt(run_location + '/Model_' + str(int(min_mod['run_id'])).zfill(4) + '/He' + str(float(min_mod['he'])) + '_CNO7.5/hjd' + str(hjd).ljust(13, '0') + '_' + line + '.txt').T
    bounds = line_bounds[line] + np.array([-10, 10])

    inds_obs = (x[0] >= bounds[0]) * (x[0] <= bounds[1])
    axs[ind].plot(x[0][inds_obs], x[1][inds_obs], 'k')

    inds_mod = (y[0] >= bounds[0]) * (y[0] <= bounds[1])
    axs[ind].plot(y[0][inds_mod], y[1][inds_mod], 'r')

    axs[ind].set_xlim(bounds)
    axs[ind].set_xlabel('Wavelength ($\AA$)')

    w = y[0][inds_mod]

    chi_mods = []
    for mod in mods:
        z = np.loadtxt(run_location + '/Model_' + str(int(mod['run_id'])).zfill(4) + '/He' + str(float(mod['he'])) + '_CNO7.5/hjd' + str(hjd).ljust(13, '0') + '_' + line + '.txt').T
        chi_mods.append(np.interp(w, z[0], z[1]))

    upper = np.max(chi_mods, axis=0)
    lower = np.min(chi_mods, axis=0)

    axs[ind].fill_between(w, lower, upper, color='red', alpha = 0.3)

axs[0].set_ylabel('Scaled Flux')
plt.show()
