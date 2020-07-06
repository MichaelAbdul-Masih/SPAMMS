import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scInterp
import sys
import glob
from scipy.stats import stats
import sys
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
for param in fit_params:
    print param
    fit_params_chi[param] = [np.min(data['chi2'][data[param] == unique]) for unique in fit_params_unique[param]]
    try:
        fit_param_interp[param] = scInterp.interp1d(fit_params_unique[param], fit_params_chi[param], kind='cubic')
    except:
        fit_param_interp[param] = scInterp.interp1d(fit_params_unique[param], fit_params_chi[param], kind='quadratic')
    fit_param_range[param] = np.linspace(np.min(fit_params_unique[param]), np.max(fit_params_unique[param]), 50.)

    temp_pr = np.linspace(np.min(fit_params_unique[param]), np.max(fit_params_unique[param]), 200.)
    temp_chi = fit_param_interp[param](temp_pr)
    ind_min = np.argmin(temp_chi)
    fit_param_min[param] = temp_pr[ind_min]

    chi_sig = temp_chi[ind_min]*(1. + np.sqrt(2./n_dof))
    fit_param_sig[param] = chi_sig
    sig_range_ind = (temp_chi < chi_sig)
    sig_range = [min(temp_pr[sig_range_ind]), max(temp_pr[sig_range_ind])]
    print sig_range
    # # print chi_sig
    # Ps = [stats.distributions.chi2.sf(i*n_dof/temp_chi[ind_min], n_dof) for i in temp_chi]
    # inds = (np.array(Ps) >= 0.05)
    # chi_sig2 = max(temp_chi[inds])
    # fit_param_sig_ga[param] = chi_sig2
    # print chi_sig2



i = np.argmin(data['chi2'])
try:
    min_mod = 'Model' + str(int(data['run_id'][i])).zfill(4) + ' He' + str(float(data['he'][i])) + '_C' + str(float(data['c'][i])) + '_N' + str(float(data['n'][i])) + '_O' + str(float(data['o'][i]))
except:
    min_mod = 'Model' + str(int(data['run_id'][i])).zfill(4) + ' He' + str(float(data['he'][i])) + '_CNO' + str(float(data['cno'][i]))
print min_mod

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
pot = np.unique(data[:, 1])
Ta = np.unique(data[:, 2])
Tb = np.unique(data[:, 3])
period = np.unique(data[:, 4])
sma = np.unique(data[:, 5])
inc = np.unique(data[:, 6])
q = np.unique(data[:, 7])
t0 = np.unique(data[:, 8])
he = np.unique(data[:, 9])
cno = np.unique(data[:, 10])

pot_chi = np.zeros_like(pot)
Ta_chi = np.zeros_like(Ta)
Tb_chi = np.zeros_like(Tb)
period_chi = np.zeros_like(period)
sma_chi = np.zeros_like(sma)
q_chi = np.zeros_like(q)
t0_chi = np.zeros_like(t0)
he_chi = np.zeros_like(he)
cno_chi = np.zeros_like(cno)

for i in range(len(pot)):
    pot_chi[i] = np.min(data[data[:, 1] == pot[i]][:, 0])
for i in range(len(Ta)):
    Ta_chi[i] = np.min(data[data[:, 2] == Ta[i]][:, 0])
for i in range(len(Tb)):
    Tb_chi[i] = np.min(data[data[:, 3] == Tb[i]][:, 0])
for i in range(len(period)):
    period_chi[i] = np.min(data[data[:, 4] == period[i]][:, 0])
for i in range(len(sma)):
    sma_chi[i] = np.min(data[data[:, 5] == sma[i]][:, 0])
for i in range(len(q)):
    q_chi[i] = np.min(data[data[:, 7] == q[i]][:, 0])
for i in range(len(t0)):
    t0_chi[i] = np.min(data[data[:, 8] == t0[i]][:, 0])
for i in range(len(he)):
    he_chi[i] = np.min(data[data[:, 9] == he[i]][:, 0])
for i in range(len(cno)):
    cno_chi[i] = np.min(data[data[:, 10] == cno[i]][:, 0])

print len(he), len(he_chi)


# potINT = scInterp.interp1d(pot, pot_chi, kind='quadratic')
TaINT = scInterp.interp1d(Ta, Ta_chi, kind='cubic')
TbINT = scInterp.interp1d(Tb, Tb_chi, kind='cubic')
# periodINT = scInterp.interp1d(period, period_chi, kind='cubic')
# smaINT = scInterp.interp1d(sma, sma_chi, kind='cubic')
# qINT = scInterp.interp1d(q, q_chi, kind='cubic')
# t0INT = scInterp.interp1d(t0, t0_chi, kind='cubic')
heINT = scInterp.interp1d(he, he_chi, kind='cubic')
cnoINT = scInterp.interp1d(cno, cno_chi, kind='cubic')

# potNEW = np.linspace(np.min(pot), np.max(pot), 50.)
TaNEW = np.linspace(np.min(Ta), np.max(Ta), 50.)
TbNEW = np.linspace(np.min(Tb), np.max(Tb), 50.)
periodNEW = np.linspace(np.min(period), np.max(period), 50.)
smaNEW = np.linspace(np.min(sma), np.max(sma), 50.)
qNEW = np.linspace(np.min(q), np.max(q), 50.)
t0NEW = np.linspace(np.min(t0), np.max(t0), 50.)
heNEW = np.linspace(np.min(he), np.max(he), 50.)
cnoNEW = np.linspace(np.min(cno), np.max(cno), 50.)


fig = plt.figure(figsize=(14, 8))
axTa = fig.add_subplot(211)
axTb = fig.add_subplot(212, sharey=axTa)

# axpot.plot(data[:, 1], data[:, 0], 'k.', alpha=0.5)
# axpot.set_xlabel('Potential')
axTa.plot(data[:, 2], data[:, 0], 'k.', alpha=0.5)
axTa.set_ylabel('Chi square',  fontsize=15)
axTa.set_xlabel('Teff (K)',  fontsize=15)
axTb.plot(data[:, 3], data[:, 0], 'k.', alpha=0.5)
axTb.set_ylabel('Chi square',  fontsize=15)
axTb.set_xlabel('Teff (K)',  fontsize=15)

# axpot.plot(pot, pot_chi, 'midnightblue', lw=2, alpha=.75)
# axpot.plot(potNEW, potINT(potNEW), 'k:', alpha=0.8, lw=1)
# axTa.plot(Ta, Ta_chi, 'indigo', lw=2, alpha=.75)
axTa.plot(TaNEW, TaINT(TaNEW), 'k:', alpha=0.8, lw=1)
# axTb.plot(Tb, Tb_chi, 'purple', lw=2, alpha=.75)
axTb.plot(TbNEW, TbINT(TbNEW), 'k:', alpha=0.8, lw=1)

plt.show()

fig = plt.figure(figsize=(14, 8))
axhe = fig.add_subplot(211)
axcno = fig.add_subplot(212, sharey=axhe)

# axpot.plot(data[:, 1], data[:, 0], 'k.', alpha=0.5)
# axpot.set_xlabel('Potential')
axhe.plot(data[:, 9], data[:, 0], 'k.', alpha=0.5)
axhe.set_ylabel('Chi square',  fontsize=15)
axhe.set_xlabel('He',  fontsize=15)
axcno.plot(data[:, 10], data[:, 0], 'k.', alpha=0.5)
axcno.set_ylabel('Chi square',  fontsize=15)
axcno.set_xlabel('CNO',  fontsize=15)

# axpot.plot(pot, pot_chi, 'midnightblue', lw=2, alpha=.75)
# axpot.plot(potNEW, potINT(potNEW), 'k:', alpha=0.8, lw=1)
# axhe.plot(he, he_chi, 'indigo', lw=2, alpha=.75)
axhe.plot(heNEW, heINT(heNEW), 'k:', alpha=0.8, lw=1)
# axcno.plot(cno, cno_chi, 'purple', lw=2, alpha=.75)
axcno.plot(cnoNEW, cnoINT(cnoNEW), 'k:', alpha=0.8, lw=1)

plt.show()
'''
