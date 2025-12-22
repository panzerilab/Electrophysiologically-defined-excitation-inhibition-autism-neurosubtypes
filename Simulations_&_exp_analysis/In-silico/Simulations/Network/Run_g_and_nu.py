"""
72022-11-23
In this script we run simulation of E-I network with different values of external input and synaptic conductance ratios
"""
import os, sys
import numpy as np
import pandas as pd
from itertools import product
import matplotlib.pyplot as plt
from scipy import signal
from fooof import FOOOF
import seaborn as sns

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from src import networks
from src.networks import ExcInhNetworkrecordrecurrent
from src import analysis
from src import tools

results_folder = os.path.join("..", "data")
parameter_db_path = os.path.join(results_folder, "metadata", "parameter_db.csv")
stats_db_path = os.path.join(results_folder, "metadata", "stats_db.csv")

Server = [["MatiRovereto", 60],
          ["GabRovereto", 60],
          ["MatiZMNH", 100]]

base_params = networks.params_pablo.copy()
base_params = tools.dictupdate(base_params, {"trials": 2})
base_params = tools.dictupdate(base_params, {"label": 'Proxy to isolate recurrent contribution'})
base_params = tools.dictupdate(base_params, {"simtime": 40500.})
base_params = tools.dictupdate(base_params, {"recstep": 1.})
base_params = tools.dictupdate(base_params, {"server": Server[0][0],
                                             "num_threads": Server[0][1]})
base_params = tools.dictupdate(base_params, {"autapse": False})
base_params = tools.dictupdate(base_params, {"multapse": False})
base_params = tools.dictupdate(base_params, {"v_0": 2.})


resulting_ids = []
freq = []
param_psd = []
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

default_ei = np.abs(base_params['exc_exc_recurrent']/base_params['inh_exc_recurrent'])

modulation_factor = [0.05]
#resting_potential = np.arange(-70., -75, step=-0.5)
#resting_potential = np.arange(-70., -65.+.01, step=0.5)
resting_potential = np.linspace(-70., -72., 50)
resting_potential_weight = resting_potential/-70

for m, el in product(modulation_factor, resting_potential):
    # Build new parameter set
    d_el = el - (-70)
    new_params = tools.dictupdate(base_params, {
       
         "exc_E_L": el,
         "inh_E_L": el,
        
        "inh_inh_recurrent": base_params["inh_inh_recurrent"] * (1. + d_el * m),
        "inh_exc_recurrent": base_params["inh_exc_recurrent"] * (1. + d_el * m),
        "exc_inh_recurrent": base_params["exc_inh_recurrent"] * (1. + d_el * m),
        "exc_exc_recurrent": base_params["exc_exc_recurrent"] * (1. + d_el * m),

        "label": f'E_l inh and modulation factor {m} for fixed gratio changes, g decreases'
    })
    ei_ratio = np.abs(new_params['inh_exc_recurrent']/new_params['exc_exc_recurrent'])
    print(f"\nModulation {m}")
    print(f"\nResting potential {el}")
    print(f"\nEI ratio {ei_ratio}")

    keys_flat = []
    values_flat = []
    params_flat = new_params.copy()
    for key, value in new_params.items():
        if type(value) is list:
            flat_val = tools.flatten_matrix(value)
            params_flat = tools.dictupdate(params_flat, {key: flat_val})
        elif type(value) is np.ndarray:
            flat_val = tools.flatten_matrix(value.tolist())
            params_flat = tools.dictupdate(params_flat, {key: flat_val})

    new_params_pd = pd.DataFrame(params_flat, index=[0])
    # Check that simulation is new
    is_there, dupentry_id = tools.check_database_entry(new_params_pd, parameter_db_path)
    existing_trials =0
    if is_there is True:
        print(f"Already simulated, ID {dupentry_id}")
        resulting_ids.append(dupentry_id)
        continue
    # construct new network
    network = ExcInhNetworkrecordrecurrent(params=new_params)
    entry = 0
    entry_id = 0
    stats = []
    freq = 0
    psd = []
    lfp = []
    I_gaba = []
    I_ampa = []
    g_exc = []
    g_inh = []
    v_m = []
    g_excrecur = []
    g_inhrecur = []
    v_mrecur = []
    eiratio = []
    for t in range(new_params["trials"]):
        print(f"\nTrial {t + 1}")
        print("\nBuilding network")
        print(20 * "#")
        network.create_network()

        print("\nSimulating network")
        print(20 * "#")
        mm_trial, mm_recurrent_trial, espikes_trial, ispikes_trial, recurspikes = network.simulate_network()
        spikes = [espikes_trial, ispikes_trial]
        ncells_var = ["N_exc", "N_inh"]

        print("\nStatistical characterization")
        print(20 * "#")
        stats_names = ["rate", "cv", "isi", "corr"]
        full_stats_names = ["exc" + sub for sub in stats_names]   + ["inh" + sub for sub in stats_names] + ["Trial"]

        trial_stats = [analysis.spike_statistics_1pop(spk, new_params["simtime"], new_params[popname],
                                                      warmup_t=new_params["warmup"])
                       for spk, popname in zip(spikes, ncells_var)]
        trial_stats = sum(trial_stats, [])
        trial_stats = trial_stats + [t+1]
        #stats_pd = pd.DataFrame(columns=full_stats_names, data=stats_values)
        stats.append(trial_stats)  # Measure statistics
        if t == 0:
            print("\nCreating parameter database entry")
            print(20 * "#")
            entry = tools.create_database_entry(params_flat, parameter_db_path)
            entry_id = entry["ID"][0]
            param_append = tools.append_database_entry(entry, parameter_db_path)
        print("\n Entry ID")
        print(entry_id)
        #print(trial_stats)
        stats_pd = pd.DataFrame(columns=full_stats_names, data=[trial_stats])
        stats_entry = tools.create_database_entry(stats_pd,
                                                  stats_db_path,
                                                  entryid=entry_id)
        stats_append = tools.append_database_entry(stats_entry, stats_db_path)


        print("\nCalculating LFP and its spectrum")
        print(20 * "#")
        # Delete warmup, decimate by factor 10
        I_ampa_trial, I_gaba_trial = analysis.estimatesynapcurr(mm_trial, decim_factor=new_params['decim_factor'])
        print(f'IAMPA shape {I_ampa_trial.shape}')
        print(f'IGABA shape {I_gaba_trial.shape}')
        I_ampa.append(I_ampa_trial)
        I_gaba.append(I_gaba_trial)
        trial_lfp = np.abs(I_ampa_trial) + np.abs(I_gaba_trial)
        print(f'LFP shape {trial_lfp.shape}')
        trial_eiratio = np.abs(I_ampa_trial) / np.abs(I_gaba_trial)
        print(f'Mean EI ratio {np.mean(trial_eiratio)}')

        eiratio.append(trial_eiratio)
        lfp.append(trial_lfp)


        mm_trial = analysis.cut_warmup_time(mm_trial, 500)
        gex = np.array(mm_trial.groupby(["times"])["g_ex"].sum())
        gin = np.array(mm_trial.groupby(["times"])["g_in"].sum())
        vm = np.array(mm_trial.groupby(["times"])["V_m"].sum())

        g_exc.append(signal.decimate(gex, new_params['decim_factor']))
        g_inh.append(signal.decimate(gin, new_params['decim_factor']))
        v_m.append(signal.decimate(vm, new_params['decim_factor']))

        mmrecur_trial = analysis.cut_warmup_time(mm_recurrent_trial, 500)
        gexrecur = np.array(mmrecur_trial.groupby(["times"])["g_ex"].sum())
        ginrecur = np.array(mmrecur_trial.groupby(["times"])["g_in"].sum())
        vmrecur = np.array(mmrecur_trial.groupby(["times"])["V_m"].sum())

        g_excrecur.append(signal.decimate(gexrecur, new_params['decim_factor']))
        g_inhrecur.append(signal.decimate(ginrecur, new_params['decim_factor']))
        v_mrecur.append(signal.decimate(vmrecur, new_params['decim_factor']))

    print("\nSaving neural activity")
    print(20 * "#")
    np.save(os.path.join(results_folder, f"ID{entry_id}"), np.array([lfp, I_ampa, I_gaba, g_exc, g_inh, v_m,
                                                                     g_excrecur, g_inhrecur, v_mrecur]))

    resulting_ids.append(entry_id)
    print(entry_id)



print(f"Created IDs: {resulting_ids}")
