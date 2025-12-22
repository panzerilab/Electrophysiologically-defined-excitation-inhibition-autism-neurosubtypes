# -*- coding: utf-8 -*-
"""
This is a collection of analysis tools, both for individual cells and populations

Created on Mon Nov  7 11:08:13 2022

@author: glorenz
"""

import os, sys
import numpy as np
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import h5py
import scipy
from itertools import product
from scipy.stats import rankdata
from scipy import integrate
from scipy import signal
from src import tools
from fooof import FOOOF
import pickle
import math
#import scipy.signal


# Spike data analysis
def cut_warmup_time(spikes, warmup_time):
    # Removes initial warmup time from recorded spikes
    spikes_cut = spikes[spikes['times'] > warmup_time]
#    spikes_cut['times'] = spikes['times'][
#        spikes['times'] > warmup_time]
    return spikes_cut


def compute_rate(spikes, n_rec, sim_time):
    # Computes average rate in sp/s from recorded spikes
    if n_rec > 0:
        return 1. * len(spikes['times']) / n_rec / sim_time * 1e3
    else:
        return 0.


def sort_spikes(spikes):
    # Sorts recorded spikes by node ID
    unique_node_ids = sorted(np.unique(spikes['senders']))
    spiketrains = []
    for node_id in unique_node_ids:
        spiketrains.append(spikes['times'][spikes['senders'] == node_id])
    return unique_node_ids, spiketrains


def compute_isidist(spiketrains, model=True):
    # Computes interspike interval distribution from
    if len(spiketrains) == 0:
        return 0., 0.
    if model:
        isi_dist = np.hstack([np.diff(st) if len(st) >= 3 else 0 for st in spiketrains])
    else:
        isi_dist = np.hstack([np.diff(st.T) if len(st.T) >= 3 else 0 for st in spiketrains])
    return isi_dist


def compute_cv2(isi_dist):
    # Computes coefficient of variation from sorted spikes
    if len(isi_dist) > 1:
        misi = np.mean(isi_dist)
        return np.std(isi_dist) / misi, misi
    else:
        return 0., 0.


def compute_cv(spiketrains, model=True):
    # Computes coefficient of variation from sorted spikes
    if len(spiketrains) == 0:
        return 0., 0.
    if model:
        isi_dist = np.hstack([np.diff(st) if len(st) >= 3 else 0 for st in spiketrains])
    else:
        isi_dist = np.hstack([np.diff(st.T) if len(st.T) >= 3 else 0 for st in spiketrains])
    if len(isi_dist) > 1:
        misi = np.mean(isi_dist)
        return np.std(isi_dist) / misi, misi
    else:
        return 0., 0.


def compute_lv(isi_dist):
    diff_isi = np.diff(isi_dist, axis=-1)
    summ_isi = isi_dist[:, 1:]+isi_dist[:, -1]
    quotient = np.flatten(diff_isi**2/summ_isi**2)
    return 3/(quotient.shape[0]-1)*np.sum(quotient)


def bin_spiketrains(spiketrains, t_min, t_max, t_bin):
    # Bins sorted spikes
    bins = np.arange(t_min, t_max, t_bin)
    return bins, np.array([np.histogram(s, bins=bins)[0] for s in spiketrains])


def compute_correlations(binned_spiketrains):
    # Computes correlations from binned spiketrains
    n = len(binned_spiketrains)
    if n > 1:
        cc = np.corrcoef(binned_spiketrains)
        return 1. / (n * (n - 1.)) * (np.sum(cc) - n)
    else:
        return 0.


def spike_statistics(espikes, ispikes, rec_time, warmup_t=500., N_e=4000, N_i=1000):
    # Computes population-averaged rates coefficients of variation and
    # correlations from recorded spikes of excitatory and inhibitory
    # populations

    if warmup_t != 0:
        espikes = cut_warmup_time(espikes, warmup_t)
        ispikes = cut_warmup_time(ispikes, warmup_t)

    eactivcells = espikes["senders"].nunique()
    iactivcells = ispikes["senders"].nunique()
    print(f"{eactivcells} exc cell spiked")
    print(f"{iactivcells} inh cell spiked")
    erate = compute_rate(espikes, eactivcells, rec_time-warmup_t)
    irate = compute_rate(espikes, iactivcells, rec_time-warmup_t)
    print(f"Rate [Hz]: r_ex = {erate}, r_in={irate}")

    enode_ids, espiketrains = sort_spikes(espikes)
    inode_ids, ispiketrains = sort_spikes(ispikes)

    ecv, eisi = compute_cv(espiketrains)
    icv, iisi = compute_cv(ispiketrains)
    print(f"Coefficient of variation: cv_ex = {ecv}, cv_in={icv}")
    print(f"Mean ISI [ms]: ISI_ex={eisi}, ISI_in={iisi}")
    N_corr = 50
    if len(espiketrains) > N_corr:
        rnd_e = np.random.choice(len(espiketrains), N_corr, replace=False)
        e_selected = [espiketrains[i] for i in rnd_e]
    else:
        e_selected = espiketrains

    if len(ispiketrains) > N_corr:
        rnd_i = np.random.choice(len(ispiketrains), N_corr, replace=False)
        i_selected = [ispiketrains[i] for i in rnd_i]
    else:
        i_selected = espiketrains

    ecorr = compute_correlations(
        bin_spiketrains(e_selected, warmup_t, rec_time, 1.)[1])
    icorr = compute_correlations(
        bin_spiketrains(i_selected, warmup_t, rec_time, 1.)[1])
    print(f"Correlation matrices: corr_ex = {ecorr}, corr_in={icorr}")

    return erate, irate, ecv, icv, eisi, iisi, ecorr, icorr


def spike_statistics_1pop(spikes, rec_time, ncells, warmup_t=500.):

    spikes = cut_warmup_time(spikes, warmup_t)
    rate = compute_rate(spikes, ncells, rec_time-warmup_t)
    print(f"Rate [Hz]: {rate} sp/s")

    node_ids, spiketrains = sort_spikes(spikes)

    cv, isi = compute_cv(spiketrains)
    print(f"Coefficient of variation: cv = {cv}")
    print(f"Mean ISI [ms]: ISI={isi}")
    N_corr = 50
    if len(spiketrains) > N_corr:
        rnd_c = np.random.choice(len(spiketrains), N_corr, replace=False)
        c_selected = [spiketrains[i] for i in rnd_c]
    else:
        c_selected = spiketrains

    corr = compute_correlations(
        bin_spiketrains(c_selected, warmup_t, rec_time, 1.)[1])

    print(f"Correlation matrix: corr_ex = {corr}")

    return [rate, cv, isi, corr]


def plot_raster_and_pop_activity(espikes, ispikes, ne_max=100, ni_max=100,
                                 t_min=0, t_max=700, binwidth=1, lbl="", title=""):
    colors = [c["color"] for c in plt.rcParams["axes.prop_cycle"]]

    espikes = espikes.loc[(t_min < espikes.times) &
                          (espikes.times < t_max)]
    ispikes = ispikes.loc[(t_min < ispikes.times) &
                          (ispikes.times < t_max)]
    espikes_small = espikes.loc[(espikes.senders < ne_max)]
    ispikes_small = ispikes.loc[(ispikes.senders < ispikes.senders.min() + ni_max)]
    fig = plt.figure(figsize=(10, 4), dpi=300)

    xlabel = "Time (ms)"
    ylabel = "Neuron ID"
    ax1 = fig.add_axes([0.1, 0.3, 0.85, 0.6])
    ax2 = fig.add_axes([0.1, 0.1, 0.85, 0.17])
    esenders_rank = rankdata(espikes_small.senders, method='dense')
    isenders_rank = rankdata(ispikes_small.senders, method='dense')
    ax1.plot(espikes_small.times, esenders_rank, ".", color=colors[0])
    ax1.plot(ispikes_small.times, isenders_rank + ne_max, ".", color=colors[1])
    # ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    ax1.set_xticks([])
    # ax1.set_yticks([])

    num_neurons = len(np.unique(espikes.senders)) + len(np.unique(ispikes.senders))
    hist_times = pd.concat([espikes.times, ispikes.times], ignore_index=True)
    weights = 1000 / (binwidth * num_neurons) * np.ones(hist_times.shape[0])
    sns.histplot(x=hist_times,
                 weights=weights,
                 binwidth=binwidth,
                 color='k', ax=ax2)
    ax2.set_ylabel("Rate (Hz)")
    ax2.set_xlabel(xlabel)
    fig.savefig(f"{lbl}.png")


def plot_raster_and_pop_activity_var_npop(spikes_list, pop_names, ncells=100,
                                          t_min=0, t_max=700, binwidth=1,
                                          lbl="", title=""):
    colors = [c["color"] for c in plt.rcParams["axes.prop_cycle"]]
    spikes_selection = []
    for pop in spikes_list:
        spk_twindow = pop.loc[(t_min < pop.times) &
                              (pop.times < t_max)]
        spikes_selection.append(spk_twindow.loc[spk_twindow.senders < (spk_twindow.senders.min() + ncells)])
    fig = plt.figure(figsize=(10, 4), dpi=300)

    xlabel = "Time (ms)"
    ylabel = "Neuron ID"
    spikes_ranked = [rankdata(pop.senders, method='dense') for pop in spikes_selection]
    ax1 = fig.add_axes([0.1, 0.3, 0.85, 0.6])
    ax2 = fig.add_axes([0.1, 0.1, 0.85, 0.17])
    for i, [pop_sel, pop_rank] in enumerate(zip(spikes_selection, spikes_ranked)):
        ax1.plot(pop_sel.times, pop_rank + i * ncells, ".", markersize=1, color=colors[i], label=pop_names[i])
    # ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    ax1.set_xticks([])
    # ax1.set_yticks([])

    num_neurons = np.sum([pop.senders.nunique() for pop in spikes_selection])
    hist_times = pd.concat([pop.times for pop in spikes_selection], ignore_index=True)
    weights = 1000 / (binwidth * num_neurons) * np.ones(hist_times.shape[0])
    sns.histplot(x=hist_times,
                 weights=weights,
                 binwidth=binwidth,
                 color='k', ax=ax2)
    ax2.set_ylabel("Rate (Hz)")
    ax2.set_xlabel(xlabel)
    fig.savefig(f"{lbl}.png")


def plot_raster_and_pop_activity_all_pops(filename, spikes_list, pop_names,
                                          t_min=0, t_max=700, binwidth=1,
                                          title=""):
    """
    Plot a rastergram of the entire neuronal population and a bar plot of the population frequency
    :param filename: name of the file to be saved
    :param spikes_list: list object of spike recordings per population
    :param pop_names: list of names for the legend
    :param t_min: initial time with the beginning of the simulation as zero
    :param t_max: final time with the beginning of the simulation as zero
    :param binwidth: time resolution to calculate the population firing
    :param title: figure title of both subplots
    :return: saves a figure and returns None
    """
    colors = [c["color"] for c in plt.rcParams["axes.prop_cycle"]]
    spikes_selection = []
    for pop in spikes_list:
        spk_twindow = pop.loc[(t_min < pop.times) &
                              (pop.times < t_max)]
        spikes_selection.append(spk_twindow)
    fig = plt.figure(figsize=(10, 4), dpi=300)

    xlabel = "Time (ms)"
    ylabel = "Neuron ID"
    ax1 = fig.add_axes([0.1, 0.3, 0.85, 0.6])
    ax2 = fig.add_axes([0.1, 0.1, 0.85, 0.17])
    for i, pop_sel in enumerate(spikes_selection):
        ax1.plot(pop_sel.times, pop_sel.senders, ".", markersize=.5, color=colors[i], label=pop_names[i])
    # ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    ax1.set_xticks([])
    # ax1.set_yticks([])

    num_neurons = np.sum([pop.senders.nunique() for pop in spikes_selection])
    hist_times = pd.concat([pop.times for pop in spikes_selection], ignore_index=True)
    weights = 1000 / (binwidth * num_neurons) * np.ones(hist_times.shape[0])
    sns.histplot(x=hist_times,
                 weights=weights,
                 binwidth=binwidth,
                 color='k', ax=ax2)
    ax2.set_ylabel("Rate (Hz)")
    ax2.set_xlabel(xlabel)
    fig.savefig(f"{filename}.png")
    return None


# TODO
def spike_synchrony(espikes, ispikes, binwidth=1.):
    num_neurons = len(np.unique(espikes.senders) +
                      np.unique(ispikes.senders))
    hist_times = pd.concat([espikes.times, ispikes.times], ignore_index=True)
    bins = np.arange(hist_times.min(), hist_times.max(), binwidth)
    hist = np.histogram(hist_times, bins=bins)[0]
    time_interval = hist_times.max() - hist_times.min()
    rate = 1000. * hist_times.shape[0] / (time_interval * num_neurons)  # rate in Hz averaged through cells
    ac = np.correlate(hist, hist, "same")
    zero_lag_index = int((ac.shape[0] - 1) / 2)
    ac = ac[zero_lag_index - 25:zero_lag_index + 25]
    ac = ac / ac.shape[0]
    print(ac.shape[0])
    print(zero_lag_index)
    print(type(zero_lag_index))
    zero_lag_index = int((ac.shape[0] - 1) / 2)
    return ac, ac[zero_lag_index] / rate ** 2


# Analysis of analog signals
def lfp(multimeter_pd, E_ex, E_in, warmup_time):
    multimeter_pd = multimeter_pd.loc[warmup_time < multimeter_pd.times]
    multimeter_pd["AMPA_curr"] = - multimeter_pd['g_ex'] * (multimeter_pd['V_m'] - E_ex)
    multimeter_pd["GABA_curr"] = multimeter_pd['g_in'] * (multimeter_pd['V_m'] - E_in)
    abs_AMPA = np.abs(multimeter_pd.groupby(["times"])["AMPA_curr"].sum())
    abs_GABA = np.abs(multimeter_pd.groupby(["times"])["GABA_curr"].sum())
    return abs_AMPA + abs_GABA


# TODO
def membrane_synchrony(mm, num_cells=100):
    cells = np.random.choice(mm.senders, size=num_cells, replace=False)
    mm_small = mm.loc[mm.senders.isin(cells)]
    mm_small.sort_values(by=['times'], inplace=True)
    ms = []
    lags = np.arange(-int(25 / 0.05) + 1, int(25 / 0.05))
    for cell1, cell2 in itertools.combinations(cells, 2):
        x = mm_small['V_m'][mm_small['senders'] == cell1]
        y = mm_small['V_m'][mm_small['senders'] == cell2]
        correlation = np.correlate(x, y, "full")
    return


def powerspectrum(signal, simtime, simstep):
    """
    Estimation of power spectral density using Welch’s method. We divided the data into
    8 overlapping segments, with 50% overlap (default)
    :param signal:
    :param simtime:
    :param simstep:
    :return:
    """
    return scipy.signal.welch(signal, fs=1000. / simstep, nperseg=(simtime / (8.0 * simstep)))


def estimatesynapcurr(mm_pd, rev_p_ampa=0., rev_p_gaba=-80., warmup_time=500, decim_factor=10):
    mm_pd = cut_warmup_time(mm_pd, warmup_time)
    mm_pd["AMPA_curr"] = mm_pd['g_ex'] * (mm_pd['V_m'] - rev_p_ampa)
    mm_pd["GABA_curr"] = mm_pd['g_in'] * (mm_pd['V_m'] - rev_p_gaba)
    I_ampa = np.array(mm_pd.groupby(["times"])["AMPA_curr"].sum())
    I_gaba = np.array(mm_pd.groupby(["times"])["GABA_curr"].sum())
    if decim_factor == 1:
        return I_ampa, I_gaba
    else:
        return signal.decimate(I_ampa, decim_factor), signal.decimate(I_gaba, decim_factor)


def estimatesynapcurrPSP(mm_pd, rev_p_ampa=0., rev_p_gaba=-80., warmup_time=500, decim_factor=10, n_dummy=10):
    mm_pd = cut_warmup_time(mm_pd, warmup_time)
    mm_pd["AMPA_curr"] = mm_pd['g_ex'] * (mm_pd['V_m'] - rev_p_ampa)
    mm_pd["GABA_curr"] = mm_pd['g_in'] * (mm_pd['V_m'] - rev_p_gaba)
    pops_names = ['pyr', 'pv', 'som', 'vip']
    mm_pd["Pop"] = pd.cut(mm_pd['senders'], 4, labels=pops_names)
    mean_signal = mm_pd.groupby(["Pop", "times"]).mean()

    I_ampa = np.array(mm_pd.groupby(["Pop", "times"])["AMPA_curr"].mean())
    I_gaba = np.array(mm_pd.groupby(["Pop", "times"])["GABA_curr"].mean())
    return signal.decimate(I_ampa, decim_factor), signal.decimate(I_gaba, decim_factor)


def estimateei(mm_pd, rev_p_ampa=0., rev_p_gaba=-80., warmup_time=500, decim_factor=10):
    mm_pd = cut_warmup_time(mm_pd, warmup_time)
    mm_pd["AMPA_curr"] = mm_pd['g_ex'] * (mm_pd['V_m'] - rev_p_ampa)
    mm_pd["GABA_curr"] = mm_pd['g_in'] * (mm_pd['V_m'] - rev_p_gaba)
    ei_raw = np.abs(mm_pd.groupby(["times"])["AMPA_curr"].sum())/np.abs(mm_pd.groupby(["times"])["GABA_curr"].sum())
    return signal.decimate(ei_raw, decim_factor)


def powerSpectrum(lfp, rec_step):
    return signal.welch(lfp, fs=1000./rec_step)


def fit_slope(freq, psd, start_freq=50, end_freq=100, plot=False, plot_name=""):
    x = freq[(freq >= start_freq) & (freq <= end_freq)]
    y = psd[(freq >= start_freq) & (freq <= end_freq)]
    x_log = np.log(x)
    y_log = np.log(y)
    a, b = np.polyfit(x_log, y_log, 1)
    if plot:
        plt.figure()
        plt.plot(x_log, y_log, label="Data")
        plt.plot(x_log, b + a * x_log, label="Fit")
        plt.xlabel("log(Frequency)")
        plt.ylabel("log(PSD)")
        plt.savefig(f"{plot_name}.png")
    return a


def postSynaptic(simtime, simstep, senders_v, data_v, pop_ex, pop_cloned_AMPA, pop_cloned_GABA, E_ex, E_in):
    PSP_AMPA = np.zeros(int(100.0/simstep))
    PSP_GABA = np.zeros(int(100.0/simstep))
    PSC_AMPA = np.zeros(int(100.0/simstep))
    PSC_GABA = np.zeros(int(100.0/simstep))

    for n in range(len(pop_cloned_AMPA)):
        pos_exc = np.where(senders_v[0]==pop_ex[n])
        pos_cloned_AMPA = np.where(senders_v[2]==pop_cloned_AMPA[n])
        pos_cloned_GABA = np.where(senders_v[3]==pop_cloned_GABA[n])

        V_dif_AMPA = (data_v[2][0]['V_m'])[pos_cloned_AMPA]-\
        (data_v[0][0]['V_m'])[pos_exc]
        V_dif_GABA = (data_v[3][0]['V_m'])[pos_cloned_GABA]-\
        (data_v[0][0]['V_m'])[pos_exc]

        current_ex_AMPA = -(data_v[0][0]['g_ex'])[pos_exc] *\
        ((data_v[0][0]['V_m'])[pos_exc] - E_ex)

        current_ex_GABA = -(data_v[0][0]['g_in'])[pos_exc] *\
        ((data_v[0][0]['V_m'])[pos_exc] - E_in)

        current_cloned_AMPA = -(data_v[2][0]['g_ex'])[pos_cloned_AMPA] *\
        ((data_v[2][0]['V_m'])[pos_cloned_AMPA] - E_ex)

        current_cloned_GABA = -(data_v[3][0]['g_in'])[pos_cloned_GABA] *\
        ((data_v[3][0]['V_m'])[pos_cloned_GABA] - E_in)

        current_dif_AMPA = current_ex_AMPA - current_cloned_AMPA
        current_dif_GABA = current_ex_GABA - current_cloned_GABA

        i=0
        # first 500 ms of the simulations are not included
        for j in np.arange(500.0,simtime-100.0,100.0):
            PSP_AMPA += V_dif_AMPA[int(j/simstep):
            int(j/simstep)+len(PSP_AMPA)]
            PSP_GABA += V_dif_GABA[int(j/simstep):
            int(j/simstep)+len(PSP_GABA)]
            PSC_AMPA += current_dif_AMPA[int(j/simstep):
            int(j/simstep)+len(PSC_AMPA)]
            PSC_GABA += current_dif_GABA[int(j/simstep):
            int(j/simstep)+len(PSC_GABA)]
            i+=1

    PSP_AMPA/=(i*len(pop_cloned_AMPA))
    PSP_GABA/=(i*len(pop_cloned_GABA))
    PSC_AMPA/=(i*len(pop_cloned_AMPA))
    PSC_GABA/=(i*len(pop_cloned_GABA))

    return [PSP_AMPA,PSP_GABA,PSC_AMPA,PSC_GABA]


# TODO
def runsaveplot(networkarch, base_params, study_params, parameter_db_path, stats_db_path, results_folder):
    """

    :param networkarch:
    :param base_params: dictionary of the values that are supposed to be the default ones in this situation
    :param study_params: dictionary of the values to vary. the keys should match the ones in base_params
    :param parameter_db_path:
    :param stats_db_path:
    :param results_folder:
    :return:
    """
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    # let's translate the study_params dict into two lists
    var_names = list(study_params.keys())
    var_values = [study_params[k] for k in var_names]
    resulting_ids = []
    for var_combination in product(*var_values):
        new_params = tools.dictupdate(base_params, {vname: var_combination[i] for i, vname in enumerate(var_names)})
        params_flat = new_params.copy()
        for key, value in new_params.items():
            if type(value) is list:
                flat_val = tools.flatten_matrix(value)
                params_flat = tools.dictupdate(params_flat, {key: flat_val})
            elif type(value) is np.ndarray:
                flat_val = tools.flatten_matrix(value.tolist())
                params_flat = tools.dictupdate(params_flat, {key: flat_val})

        new_params_pd = pd.DataFrame(params_flat, index=[0])

        # construct new network
        network = networkarch(params=new_params)
        entry_id = 0
        stats = []
        freq = 0
        psd = []
        lfp = []
        eiratio = []
        param_psd =[]
        for t in range(base_params["trials"]):
            print(f"\nTrial {t + 1}")
            print("\nBuilding network")
            print(20 * "#")
            network.create_network()

            print("\nSimulating network")
            print(20 * "#")
            mm_trial, espikes_trial, ispikes_pv_trial, ispikes_som_trial, ispikes_vip_trial = network.simulate_network()
            spikes = [espikes_trial, ispikes_pv_trial, ispikes_som_trial, ispikes_vip_trial]
            ncells_var = ["N_pyr", "N_pv", "N_som", "N_vip"]

            print("\nStatistical characterization")
            print(20 * "#")
            stats_names = ["rate", "cv", "isi", "corr"]
            full_stats_names = ["pyr" + sub for sub in stats_names] + ["pv" + sub for sub in stats_names] \
                               + ["som" + sub for sub in stats_names] + ["vip" + sub for sub in stats_names]

            trial_stats = [
                spike_statistics_1pop(spk, new_params["simtime"], new_params[popname], warmup_t=500.)
                for spk, popname in zip(spikes, ncells_var)]
            trial_stats = sum(trial_stats, [])
            # stats_pd = pd.DataFrame(columns=full_stats_names, data=stats_values)
            stats.append(trial_stats)  # Measure statistics
            if t == 0:
                print("\nCreating parameter database entry")
                print(20 * "#")
                entry = tools.create_database_entry(params_flat, parameter_db_path)
                entry_id = entry["ID"][0]
            if t == (new_params["trials"] - 1):
                print("\nSaving parameter database entry")
                print(20 * "#")
                append_stat = tools.append_database_entry(entry, parameter_db_path)

                print("\nPloting population activity")
                print(20 * "#")
                plot_raster_and_pop_activity_all_pops(spikes, ["Pyr", "PV", "SOM", "VIP"],
                                                      t_min=500, t_max=new_params["simtime"], binwidth=1,
                                                      filename="SOMFFPyr/" + str(entry_id),
                                                      title=f"Base network: {new_params_pd['label'][0]}")

            print("\nCalculating LFP and its spectrum")
            print(20 * "#")
            ## Delete warmup, decimate by factor 10
            mm_trial = cut_warmup_time(mm_trial, 500.)
            mm_trial["AMPA_curr"] = - mm_trial['g_1'] * (mm_trial['V_m'] - new_params["exc_E_rev"][0])
            mm_trial["GABA_curr"] = (mm_trial['g_2'] + mm_trial['g_3'] + mm_trial['g_4']) * \
                                    (mm_trial['V_m'] - new_params["pv_E_rev"][1])
            I_ampa = np.array(mm_trial.groupby(["times"])["AMPA_curr"].sum())
            I_gaba = np.array(mm_trial.groupby(["times"])["GABA_curr"].sum())
            I_ampa = signal.decimate(I_ampa, 10)
            I_gaba = signal.decimate(I_gaba, 10)

            trial_lfp = np.abs(I_ampa) + np.abs(I_gaba)
            trial_eiratio = np.abs(I_ampa) / np.abs(I_gaba)
            eiratio.append(trial_eiratio)
            lfp.append(trial_lfp)

            trial_freq, trial_psd = powerspectrum(trial_lfp, simtime=new_params["simtime"] - 500.,
                                                  simstep=10 * new_params["recstep"])
            freq = trial_freq
            psd.append(trial_psd)

            print("\nSaving neural activity")
            print(20 * "#")
            espikes_trial.to_hdf(os.path.join(results_folder, f"ID{entry_id}-espikes.hf5"), key=f'T{t}', mode="a")
            ispikes_pv_trial.to_hdf(os.path.join(results_folder, f"ID{entry_id}-ispikes_pv.hf5"), key=f'T{t}', mode="a")
            ispikes_som_trial.to_hdf(os.path.join(results_folder, f"ID{entry_id}-ispikes_som.hf5"), key=f'T{t}',
                                     mode="a")
            ispikes_vip_trial.to_hdf(os.path.join(results_folder, f"ID{entry_id}-ispikes_vip.hf5"), key=f'T{t}',
                                     mode="a")

        np.save(os.path.join(results_folder, f"ID{entry_id}-lfp"), np.array(lfp))
        np.save(os.path.join(results_folder, f"ID{entry_id}-ei"), np.array(eiratio))
        print("\nCalculating peak and slope with FOOOF")
        print(20 * "#")
        param_psd.append(np.median(psd, axis=0))
        fm = FOOOF(peak_width_limits=[3, 15], min_peak_height=0.05, max_n_peaks=1)
        fm.fit(freq, param_psd[-1], [3, 300])
        # exps = fm.peak_params_[-1][1]
        try:
            cfs = fm.peak_params_[-1][0]
            print(f"Peak frequency: {cfs[-1]}")
        except:
            cfs = -1
            print("No peak found")
        # Extract aperiodic parameters
        exps = fm.get_params('aperiodic_params', 'exponent')
        # Extract peak frequencies
        # cfs = fm.get_params('peak_params', 'CF')
        # Extract goodness-of-fit metrics
        r2s = fm.get_params('r_squared')
        fm.plot(save_fig=True, file_name=f"SOMFFPyr/ID{entry_id}.png")
        stats_pd = pd.DataFrame(columns=full_stats_names, data=stats)
        stats_summ = stats_pd.mean()
        stats_sem = stats_pd.std() / np.sqrt(base_params["trials"])
        stats_summ = stats_summ.append(stats_sem.add_suffix("_sem"))
        stats_summ["slope_exponent"] = exps
        stats_summ["peak_frequency"] = cfs
        stats_summ["r_squared"] = r2s
        stats_entry = tools.create_database_entry(pd.DataFrame([stats_summ]),
                                                  stats_db_path, entryid=entry_id)
        stats_append = tools.append_database_entry(stats_entry, stats_db_path)
        resulting_ids.append(entry_id)
    print(f"Created IDs: {resulting_ids}")
    return None


def rerunanalysis(parameter_db_path, stats_db_path, results_folder, id_list=[]):
    param_db = pd.read_csv(parameter_db_path)
    stats_db = pd.read_csv(stats_db_path)

    if len(id_list)==0:
        selected_sims = param_db
    else:
        selected_sims = param_db[param_db["ID"].isin(id_list)]

    stats_names = ["erate", "irate", "ecv", "icv",
                   "eisi", "iisi", "ecorr", "icorr"]

    for index, row in selected_sims.iterrows():
        espikes = h5py.File(os.path.join(results_folder,f"ID{row['ID']}-espikes.hf5"), 'r')
        # espikes.visit(printname)
        print(type(espikes['T1']))
        print(espikes['T1'].keys())
        print(espikes.name)
        print(espikes['T1'].name)
        print(espikes['T1']['axis0'])
        print(espikes['T1']['block0_values'][:10])
        print(espikes['T1']['block1_values'][:10])
        ispikes = h5py.File(os.path.join(results_folder, f"ID{row['ID']}-ispikes.hf5"))
        lfp = np.load(os.path.join(results_folder, f"ID{row['ID']}-lfp.npy"))
        stats = []
        psd=[]
        for trial in range(row["trials"]):
            esp = pd.DataFrame.from_dict({'senders': np.array(espikes[f'T{trial}']['block0_values']).flatten(),
                                          'times': np.array(espikes[f'T{trial}']['block1_values']).flatten()})
            isp = pd.DataFrame.from_dict({'senders': np.array(ispikes[f'T{trial}']['block0_values']).flatten(),
                                          'times': np.array(ispikes[f'T{trial}']['block1_values']).flatten()})
            stats.append(spike_statistics(esp, isp,
                                          row["simtime"], warmup_t=500.,
                                          N_e=4000, N_i=1000))  # Measure statistics

            decim_factor = 10
            trial_freq, trial_psd = powerspectrum(lfp[trial], simtime=row["simtime"] - 500.,
                                                  simstep=decim_factor * row["recstep"])
            freq = trial_freq
            psd.append(trial_psd)
        print("\nCalculating peak and slope with FOOOF")
        print(20 * "#")
        psd_mean = np.mean(psd, axis=0)
        fm = FOOOF(peak_width_limits=[3, 15], min_peak_height=0.05, max_n_peaks=1)
        fm.fit(freq, psd_mean, [5, 200])
        # exps = fm.peak_params_[-1][1]
        try:
            cfs = fm.peak_params_[-1][0]
            print(f"Peak frequency: {cfs}")
        except:
            cfs = -1
            print("No peak found")
        # Extract aperiodic parameters
        exps = fm.get_params('aperiodic_params', 'exponent')
        # Extract goodness-of-fit metrics
        r2s = fm.get_params('r_squared')
        fm.plot(save_fig=True, file_name=f"ID{row['ID']}.png")
        stats_pd = pd.DataFrame(columns=stats_names, data=stats)
        stats_summ = stats_pd.mean()
        stats_sem = stats_pd.std() / np.sqrt(row["trials"])
        stats_summ = stats_summ.append(stats_sem.add_suffix("_sem"))
        stats_summ["slope_exponent"] = exps
        stats_summ["peak_frequency"] = cfs
        stats_summ["r_squared"] = r2s

        stats_entry = tools.create_database_entry(pd.DataFrame([stats_summ]),
                                                  stats_db_path, entryid=row['ID'])
        print(stats_entry)
        stats_append = tools.append_database_entry(stats_entry, stats_db_path)
        print(stats_append)
        file_stats = os.stat(stats_db_path)
        print(file_stats)
        print(f'File Size in Bytes is {file_stats.st_size}')
        print(f'File Size in MegaBytes is {file_stats.st_size / (1024 * 1024)}')
    return


def ballon_solver(t, Y, Z):
    z = Z[int(t)]
    s, f, v, q = Y
    kappa = 0.65
    gamma = 0.41
    tau = 0.98
    alpha = .32
    rho = 0.34

    ds_dt = z - kappa*s - gamma * (f - 1)
    df_dt = s
    dv_dt = 1/tau*(f-np.power(v, 1/alpha))
    dq_dt = f/rho * (1-np.power(1-rho,1/f)) - q*np.power(v, 1/alpha-1)
    return [ds_dt, df_dt, dv_dt, dq_dt]


def compute_bold_balloon(Z):
    rho = 0.34
    V0 = 0.02
    k1 = 7. * rho
    k2 = 2.
    k3 = 2. * rho - 0.2
    print(len(Z))
    T = np.arange(len(Z))
    s0 = 0.
    f0 = 1.
    v0 = 1.
    q0 = 1.
    Y0 = [s0, f0, v0, q0]
    sol = scipy.integrate.solve_ivp(ballon_solver, [np.min(T), np.max(T)], y0=Y0, max_step=len(Z), args=[Z])
    s, f, v, q = sol.y
    t = sol.t
    print(s.shape)
    bold = V0 * (k1 * (1-q) + k2 * (1-q/v) + k3 * (1-v))
    return t, bold


def compute_bold_kernel(Z):
    rho = 0.34
    V0 = 0.02
    k1 = 7. * rho
    k2 = 2.
    k3 = 2. * rho - 0.2
    print(len(Z))
    T = np.arange(len(Z))
    s0 = 0.
    f0 = 1.
    v0 = 1.
    q0 = 1.
    Y0 = [s0, f0, v0, q0]
    Y = scipy.integrate.solve_ivp(ballon_solver, [np.min(T), np.max(T)], y0=Y0, max_step=len(Z), args=[Z])
    s, f, v, q = Y

    bold = V0 * (k1 * (1-q) + k2 * (1-q/v) + k3 * (1-v))
    return bold


def bold_hrf(lfp_signal):
    HRF = pickle.load(open(os.path.join('data', 'HRF'), "rb"), encoding='latin1')
    # Downsample the HRF. Resulting rate: 2ms
    HRF = signal.decimate(HRF, 10)
    HRF = signal.decimate(HRF, 2)
    return np.convolve(lfp_signal, HRF)


def lfp_proxies(ampa_current, gaba_current, rec_step=2.):
    LFP = []
    # Sum of AMPA PSCs
    LFP.append(ampa_current)
    # Sum of GABA PSCs (change sign)
    LFP.append(-gaba_current)
    # Sum of AMPA and GABA PSCs (change sign)
    I1 = ampa_current+gaba_current
    LFP.append(-I1)
    # Sum of absolute values of AMPA and GABA PSCs
    I2 = np.abs(ampa_current) + np.abs(gaba_current)
    LFP.append(I2)
    # Reference Weighted Sum (RWS)
    AMPA_delayed = np.roll(ampa_current, int(6.0 / rec_step))
    I3 = AMPA_delayed[:len(gaba_current)]-1.65*gaba_current
    LFP.append(I3)
    return LFP


def erws1(ampa_current, gaba_current, params, rec_step=2.):
    ampa_delayed = np.roll(ampa_current, int(params[0]/rec_step))
    gaba_delayed = np.roll(gaba_current, int(params[1]/rec_step))
    return ampa_delayed - params[2]*gaba_delayed


def erws2(ampa_current, gaba_current, params, v_0=0, rec_step=2.):
    tau_AMPA = params[0] * np.power(v_0,-params[1]) + params[2]
    tau_GABA = params[3] * np.power(v_0,-params[4]) + params[5]
    alpha = params[6] * np.power(v_0,-params[7]) + params[8]
    ampa_delayed = np.roll(ampa_current, int(tau_AMPA/rec_step))
    gaba_delayed = np.roll(gaba_current, int(tau_GABA/rec_step))
    return ampa_delayed - alpha*gaba_delayed


def eeg_proxies(ampa_current, gaba_current, v_0, rec_step=2.):
    # Parameters of the ERWS-1 and ERWS-2 proxies
    ERWS1_params_causal = [0., 3.1, 0.1]
    ERWS2_params_causal = [0.,  0.,  0., -1.5,  0.2,  4.,  0.5,  0.5,  0.]
    ERWS1_params_noncausal = [-0.9, 2.3, 0.3]
    ERWS2_params_noncausal = [-0.6, 0.1, -0.4, -1.9, 0.6, 3.0, 1.4, 1.7, 0.2]
    EEG = []
    ERWS1_causal = erws1(ampa_current, gaba_current, ERWS1_params_causal, rec_step=rec_step)
    EEG.append(ERWS1_causal)
    ERWS2_causal = erws2(ampa_current, gaba_current, ERWS2_params_causal, v_0=v_0, rec_step=rec_step)
    EEG.append(ERWS2_causal)
    ERWS1_noncausal = erws1(ampa_current, gaba_current, ERWS1_params_noncausal, rec_step=rec_step)
    EEG.append(ERWS1_noncausal)
    ERWS2_noncausal = erws2(ampa_current, gaba_current, ERWS2_params_noncausal, v_0=v_0, rec_step=rec_step)
    EEG.append(ERWS2_noncausal)
    return EEG

def fmri_proxies():
    return
import cmath
from scipy import optimize
def piecewise_linear(x, y1,y2, k1, k2):
    fx0 = np.log10(30.0)
    return np.piecewise(x, [x < fx0], [lambda x:k1*x + y1, lambda x:k2*x + y2])


def bold_hrf_conv(lfp_signal, hrf, simstep=2):
    # HRF = pickle.load(open(os.path.join('data', 'HRF'), "rb"), encoding='latin1')
    # # Downsample the HRF. Resulting rate: 2ms
    # HRF = signal.decimate(HRF, 10)
    # HRF = signal.decimate(HRF, 2)
    fft_LFP = np.fft.rfft(lfp_signal)
    freq_fft_LFP = np.fft.rfftfreq(len(lfp_signal), d=simstep / 1000.)

    # Correct LFP spectra for low freq.
    # Piecewise linear regression
    pwlf_range = [np.where(freq_fft_LFP >= 1.)[0][0], np.where(freq_fft_LFP >= 99)[0][0]]
    pp, ee = optimize.curve_fit(piecewise_linear, np.log10(freq_fft_LFP[pwlf_range[0]:pwlf_range[1]]),
                                np.log10(np.abs(fft_LFP[pwlf_range[0]:pwlf_range[1]])))

    surrogate_slopes = [pp[-2], pp[-1]]
    new_fft_LFP = fft_LFP.copy()

    # Minimum freq. considered valid
    min_freq_fft_LFP = freq_fft_LFP[np.where(freq_fft_LFP >= 1.)[0][0]]
    min_fft_LFP = fft_LFP[np.where(freq_fft_LFP >= 1.)[0][0]]

    # Extrapolate low freq.
    for ffi, pbpf in enumerate(freq_fft_LFP):
        if pbpf < min_freq_fft_LFP:
            ex_slope = surrogate_slopes[0]
            # # Random phase
            # fft_phase = np.exp(1j*np.random.uniform(0, 2*np.pi, (1)))
            # Same phase of the old FFT
            fft_phase = np.exp(1j * cmath.phase(fft_LFP[ffi]))

            # Freq. at 0 Hz
            if pbpf == 0:
                pbpf = freq_fft_LFP[1]

            fft_LFP_log_real = np.log10(np.abs(np.real(min_fft_LFP))) - (np.log10(min_freq_fft_LFP) -
                                                                         np.log10(pbpf)) * ex_slope
            fft_LFP_log_imag = np.log10(np.abs(np.imag(min_fft_LFP))) - (np.log10(min_freq_fft_LFP) -
                                                                         np.log10(pbpf)) * ex_slope

            new_fft_LFP[ffi] = np.abs(np.power(10, fft_LFP_log_real) + 1j * np.power(10, fft_LFP_log_imag)) * fft_phase

    # Invert LFP spectrum
    #I1 = np.fft.irfft(new_fft_LFP, len(HRF))
    fft_LFP = new_fft_LFP

    # High-pass filter (HPF)
    slope1 = 2
    slope2 = 1.5
    cutoff_freq1 = 0.03
    cutoff_freq2 = 20

    fft_HPF = np.ones(len(freq_fft_LFP)) * (0 + 0j)
    cutoff_pos1 = np.where(freq_fft_LFP >= cutoff_freq1)[0][0]
    cutoff_pos2 = np.where(freq_fft_LFP >= cutoff_freq2)[0][0]

    # First segment
    for pbpf in np.arange(0, cutoff_pos1):
        # fft_phase = np.exp(1j*0.5*np.pi)
        fft_HPF[pbpf] = freq_fft_LFP[pbpf] * slope1

    # Second segment
    for pbpf in np.arange(cutoff_pos1 - 1, cutoff_pos2):
        fft_HPF[pbpf] = np.power(freq_fft_LFP[pbpf], slope2) + fft_HPF[cutoff_pos1 - 1] - \
                        np.power(freq_fft_LFP[cutoff_pos1 - 1], slope2)

    # Plateau
    for pbpf in np.arange(cutoff_pos2, len(freq_fft_LFP)):
        fft_HPF[pbpf] = fft_HPF[cutoff_pos2 - 1]

    # Offset of the HPF
    fft_HPF *= 10 ** (-2)
    # HRF spectrum and offset
    fft_HRF = np.fft.rfft(hrf)
    # fft_HRF *= 10 ** (-3)

    # White noise
    # fft_noise = 10 ** (-5) * np.ones(len(fft_LFP))
    fft_noise = 10 ** (-2) * np.ones(len(fft_LFP))

    # BOLD signal
    print(f"fft_LFP shape {fft_LFP.shape}")
    print(f"fft_noise shape {fft_noise.shape}")
    print(f"fft_HPF shape {fft_HPF.shape}")
    print(f"fft_HPF shape {fft_HRF.shape}")
    fft_BOLD = fft_LFP * (fft_HRF + fft_noise) * fft_HPF
    # Offset
    # fft_BOLD *= 10 ** 3

    # Time-series of the BOLD signal
    recovered_BOLD = np.fft.irfft(fft_BOLD, len(hrf))

    # Downsample BOLD signal
    # BOLD_decimate = 500
    # recovered_BOLD = tools.decimate(recovered_BOLD, BOLD_decimate)
    print("len(BOLD) = %s" % len(recovered_BOLD))
    return recovered_BOLD


def get_s_peak(tau_m, tau_r, tau_d):
    return (tau_m/tau_d) * (tau_r/tau_d)**(tau_r/(tau_d - tau_r))

