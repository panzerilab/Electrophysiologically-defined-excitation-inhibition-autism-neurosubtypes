"""
72022-11-23
In this script we run simulation of E-I networks with camk manipulation
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

# ------------------------------------------------------------------
# ------------- CHANGE ONLY THIS PART (like HM4Di) -----------------
# ------------------------------------------------------------------

results_folder = os.path.join("..", "Simulation results", "camk")
metadata_folder = os.path.join(results_folder, "metadata")
data_folder     = os.path.join(results_folder, "data")

# Create folders if missing
os.makedirs(metadata_folder, exist_ok=True)
os.makedirs(data_folder, exist_ok=True)

parameter_db_path = os.path.join(metadata_folder, "parameter_db.csv")
stats_db_path     = os.path.join(metadata_folder, "stats_db.csv")

# ------------------------------------------------------------------

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
base_params = tools.dictupdate(base_params, {"v_0": 3.})

resulting_ids = []
freq = []
param_psd = []

default_ei = np.abs(base_params['exc_exc_recurrent']/base_params['inh_exc_recurrent'])

modulation_factor = [0.05]
resting_potential = np.linspace(-70., -65., 50)
resting_potential_weight = resting_potential / -70

for m, el in product(modulation_factor, resting_potential):

    d_el = el - (-70)
    new_params = tools.dictupdate(base_params, {

        "exc_E_L": el

    })

    ei_ratio = np.abs(new_params['inh_exc_recurrent'] / new_params['exc_exc_recurrent'])
    print(f"\nModulation {m}")
    print(f"Resting potential {el}")
    print(f"EI ratio {ei_ratio}")

    # Flatten lists/arrays for DB
    params_flat = new_params.copy()
    for key, value in new_params.items():
        if type(value) is list:
            params_flat[key] = tools.flatten_matrix(value)
        elif type(value) is np.ndarray:
            params_flat[key] = tools.flatten_matrix(value.tolist())

    new_params_pd = pd.DataFrame(params_flat, index=[0])

    # Check simulation existence
    is_there, dupentry_id = tools.check_database_entry(new_params_pd, parameter_db_path)
    if is_there:
        print(f"Already simulated, ID {dupentry_id}")
        resulting_ids.append(dupentry_id)
        continue

    # Simulate
    network = ExcInhNetworkrecordrecurrent(params=new_params)

    entry_id = 0
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
        network.create_network()

        mm_trial, mm_recurrent_trial, espikes_trial, ispikes_trial, recurspikes = network.simulate_network()
        spikes = [espikes_trial, ispikes_trial]
        ncells_var = ["N_exc", "N_inh"]

        stats_names = ["rate", "cv", "isi", "corr"]
        full_stats_names = ["exc" + s for s in stats_names] + \
                           ["inh" + s for s in stats_names] + ["Trial"]

        # Compute statistics
        trial_stats = [
            analysis.spike_statistics_1pop(spk, new_params["simtime"], new_params[pop],
                                           warmup_t=new_params["warmup"])
            for spk, pop in zip(spikes, ncells_var)
        ]
        trial_stats = sum(trial_stats, [])
        trial_stats.append(t + 1)

        if t == 0:
            entry = tools.create_database_entry(params_flat, parameter_db_path)
            entry_id = entry["ID"][0]
            tools.append_database_entry(entry, parameter_db_path)

        stats_pd = pd.DataFrame(columns=full_stats_names, data=[trial_stats])
        stats_entry = tools.create_database_entry(stats_pd, stats_db_path, entryid=entry_id)
        tools.append_database_entry(stats_entry, stats_db_path)

        # LFP + signals
        I_ampa_trial, I_gaba_trial = analysis.estimatesynapcurr(
            mm_trial, decim_factor=new_params['decim_factor'])

        I_ampa.append(I_ampa_trial)
        I_gaba.append(I_gaba_trial)

        trial_lfp = np.abs(I_ampa_trial) + np.abs(I_gaba_trial)
        lfp.append(trial_lfp)

        trial_eiratio = np.abs(I_ampa_trial) / np.abs(I_gaba_trial)
        eiratio.append(trial_eiratio)

        # Cut warm startup
        mm_trial = analysis.cut_warmup_time(mm_trial, 500)
        mmrecur_trial = analysis.cut_warmup_time(mm_recurrent_trial, 500)

        g_exc.append(signal.decimate(np.array(mm_trial.groupby("times")["g_ex"].sum()),
                                     new_params['decim_factor']))
        g_inh.append(signal.decimate(np.array(mm_trial.groupby("times")["g_in"].sum()),
                                     new_params['decim_factor']))
        v_m.append(signal.decimate(np.array(mm_trial.groupby("times")["V_m"].sum()),
                                   new_params['decim_factor']))

        g_excrecur.append(signal.decimate(np.array(mmrecur_trial.groupby("times")["g_ex"].sum()),
                                          new_params['decim_factor']))
        g_inhrecur.append(signal.decimate(np.array(mmrecur_trial.groupby("times")["g_in"].sum()),
                                          new_params['decim_factor']))
        v_mrecur.append(signal.decimate(np.array(mmrecur_trial.groupby("times")["V_m"].sum()),
                                        new_params['decim_factor']))

    # ------------------------------------------------------------------
    # -------- SAVE DATA IN /Simulation results/camk/data --------------
    # ------------------------------------------------------------------
    print("\nSaving neural activity")
    np.save(os.path.join(data_folder, f"ID{entry_id}.npy"),
            np.array([lfp, I_ampa, I_gaba, g_exc, g_inh, v_m,
                      g_excrecur, g_inhrecur, v_mrecur]))

    resulting_ids.append(entry_id)
    print(entry_id)

print(f"Created IDs: {resulting_ids}")
