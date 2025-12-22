# -*- coding: utf-8 -*-
"""
All the main models and architectures used as classes. This code is planned for NEST 3.2

"""
import os
import numpy as np
import scipy as sc
from scipy import signal
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

# Properties of plots
DPI = 300  # journals typically requires a minimum of 300 dpi
font_1 = {'fontname': 'Arial', 'fontsize': 8}
font_2 = {'fontname': 'Arial', 'fontsize': 10}
font_3 = {'fontname': 'Arial', 'fontsize': 12}
FIGDIM = (16, 7)


def dictupdate(dict1, dict2):
    """
    Update (safely) dictionaries

    :param dict1:
    :param dict2:
    :return:
    """
    assert (isinstance(dict1, dict))
    assert (isinstance(dict2, dict))

    tmp = dict1.copy()
    for key, value in dict2.items():
        if key in dict1:
            tmp[key] = value
        else:
            print(f"Check the varable name '{key}'")
    return tmp


def flatten_matrix(mult_list):
    if type(mult_list[0]) is list:
        mult_list = [item for sublist in mult_list for item in sublist] #sum(mult_list, [])
    return ",".join(map(str, mult_list))


def unflatten_matrix(string, npop):
    values = np.asarray([float(v) for v in string.split(",")])
    if values.shape[0] == npop:
        return values
    else:
        return values.reshape((int(values.shape[0]/npop), npop))


def ounoise(tsteps, simstep, sigma, tau):
    # Ornstein-Uhlenbeck process (solved by Euler method)
    # Define renormalized variables (to avoid recomputing these constants
    # at every time step)
    ou_sigma_bis = sigma * np.sqrt(2. / tau)
    ou_sqrtdt = np.sqrt(simstep)
    ou_x = np.zeros(tsteps)  # OU output
    for i in range(tsteps - 1):
        ou_x[i + 1] = ou_x[i] + simstep * (-ou_x[i] / tau) + \
                      ou_sigma_bis * ou_sqrtdt * np.random.randn()
    return ou_x


def loadData(experiment_id, filename, extension):
    data = pickle.load(open('../results/' + experiment_id + '/' + filename + extension, "rb"))
    return data


def saveData(experiment_id, filename, extension, data):
    pickle.dump(data, open('../results/' + experiment_id + '/' + filename + extension, "wb"))

'''
def decimate(x, q=10, n=4, k=0.8, filterfun=signal.cheby1):
    """
    scipy.signal.decimate like downsampling using filtfilt instead of lfilter,
    and filter coeffs from butterworth or chebyshev type 1.

    Parameters
    ----------
    x : ndarray
        Array to be downsampled along last axis.
    q : int
        Downsampling factor.
    n : int
        Filter order.
    k : float
        Aliasing filter critical frequency Wn will be set as Wn=k/q.
    filterfun : function
        `scipy.signal.filter_design.cheby1` or
        `scipy.signal.filter_design.butter` function

    Returns
    -------
    ndarray
        Downsampled signal.

    """
    if not isinstance(q, int):
        raise TypeError("q must be an integer")

    if n is None:
        n = 1

    if filterfun == signal.butter:
        b, a = filterfun(n, k / q)
    elif filterfun == signal.cheby1:
        b, a = filterfun(n, 0.05, k / q)
    else:
        raise Exception('only scipy.signal.butter or scipy.signal.cheby1 supported')

    try:
        y = signal.filtfilt(b, a, x)
    except:  # Multidim array can only be processed at once for scipy >= 0.9.0
        y = []
        for data in x:
            y.append(signal.filtfilt(b, a, data))
        y = np.array(y)

    try:
        return y[:, ::q]
    except:
        return y[::q]
'''


def flt2str(number, precision=1):
    number = str(round(number, precision)).split(".")
    if len(number) == 2:
        integer, decimal = number
    else:
        integer = number[0]
        decimal = "0"
    if precision <= len(decimal):
        return integer + 'c' + decimal[:precision]
    else:
        return integer + 'c' + decimal + (precision - len(decimal)) * '0'


def str2flt(string):
    number_parts = str(string).split("c")
    number = float(number_parts[0]) + float(number_parts[1]) * 10 ** (-len(number_parts[1]))
    return round(number, len(number_parts[1]))


def cleandir():
    dir_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
    for f in os.listdir(dir_path):
        if 'dat' in f or 'gdf' in f:
            os.remove(os.path.join(dir_path, f))


def make_filename(parameter_dictionary, indexes, precisions):
    filename = ''
    for key, index, prec in zip(parameter_dictionary.keys(), indexes, precisions):
        filename += key + flt2str(parameter_dictionary[key][index], prec)
    return filename


def create_database_entry(entry, database_filepath, entryid=None):
    if not isinstance(entry, pd.DataFrame):
        entry = pd.DataFrame(entry, index=[0])
    if os.path.exists(database_filepath):
        old_database = pd.read_csv(database_filepath)
        if "ID" not in old_database.columns:
            old_database["ID"] = pd.Series()
            last_exp_index = 0
        else:
            last_exp_index = old_database.ID.max()
        if entryid is None:
            entry["ID"] = [last_exp_index + 1]
        else:
            entry["ID"] = [entryid]
        # new_database = pd.concat([old_database, entry], ignore_index=True)
        # new_database.to_csv(database_filepath, index=False)
        return entry#["ID"][0]
    else:
        entry["ID"] = [1]
        entry.to_csv(database_filepath, index=False)
        return entry#["ID"][0]


def append_database_entry(entry, database_filepath):
    if not isinstance(entry, pd.DataFrame):
        print("This is not a Pandas Dataframe")
        return -1
    if os.path.exists(database_filepath):
        old_database = pd.read_csv(database_filepath)
        new_database = pd.concat([old_database, entry], ignore_index=True)
        new_database.to_csv(database_filepath, index=False)
        return 0


def check_database_entry(entry, database_filepath):
    try:
        database = pd.read_csv(database_filepath)
        database = database.round(6)
    except FileNotFoundError:
        print("No database with that name")
        return False, -1
    #database = database.loc[:, database.columns != 'ID']
    entry_pd = pd.DataFrame(entry, index=[0])
    entry_pd = entry_pd.round(6)
    #values = np.array([v for k, v in params.items()])

    check_duplicates = pd.merge(database, entry_pd,
                                on=[i for i in database.columns if ((i != 'ID') &
                                                                    (i != 'server') &
                                                                    (i != 'num_threads') &
                                                                    (i != 'trials')
                                                                    )])
    if check_duplicates.shape[0] == 1:
        print(f"Parameters already simulated in entry {check_duplicates.iloc[0]['ID']}")
        return True, check_duplicates.iloc[0]['ID']#, check_duplicates.iloc[0]['trials']
    elif check_duplicates.shape[0] > 1:
        print("Something wrong, there are multiple entries with the same parameters")
        return True, list(check_duplicates['ID'].to_numpy())#, check_duplicates['trials']
    else:
        print("How original, new parameters to see")
        return False, -1#, -1
# params_pd = pd.DataFrame(params_pablo, index=[0])


def add_trials(entry_id, new_trials, results_filepath, new_recording, parameter_filepath):
    param_db = pd.read_csv(parameter_filepath)
    existing_trials = param_db[param_db.ID == entry_id]["trials"]
    print(type(existing_trials.iloc[0]))
    print(f"now there are {existing_trials.iloc[0]+new_trials}trials ")
    param_db["trials"][param_db.ID == entry_id] = existing_trials.iloc[0] + new_trials
    param_db.to_csv(parameter_filepath, index=False)
    # stats_db = pd.read_csv(stats_filepath)
    old_recording = np.load(os.path.join(results_filepath, f"ID{entry_id}.npy"), allow_pickle=True)

    # old_stats = stats_db['ID' == ID]
    # old_params = param_db['ID' == ID]
    # old_ntrials = old_params['ntrials'][0]

    #new_trials = new_params['ntrials'][0]

    # tot_trials = old_ntrials + ntrials
    # tot_stats = old_ntrials * old_stats + ntrials * new_stats
    # tot_stats /= tot_trials

    tot_recording = np.concatenate([old_recording, new_recording], axis=1)
    print(tot_recording.shape)
    np.save(os.path.join(results_filepath, f"ID{entry_id}.npy"), tot_recording)
    return


def decimate(x, q=10, n=4, k=0.8, filterfun=signal.cheby1):
    """
    scipy.signal.decimate like downsampling using filtfilt instead of lfilter,
    and filter coeffs from butterworth or chebyshev type 1.

    Parameters
    ----------
    x : ndarray
        Array to be downsampled along last axis.
    q : int
        Downsampling factor.
    n : int
        Filter order.
    k : float
        Aliasing filter critical frequency Wn will be set as Wn=k/q.
    filterfun : function
        `scipy.signal.filter_design.cheby1` or
        `scipy.signal.filter_design.butter` function

    Returns
    -------
    ndarray
        Downsampled signal.

    """
    if not isinstance(q, int):
        raise TypeError("q must be an integer")

    if n is None:
        n = 1

    if filterfun == signal.butter:
        b, a = filterfun(n, k / q)
    elif filterfun == signal.cheby1:
        b, a = filterfun(n, 0.05, k / q)
    else:
        raise Exception('only scipy.signal.butter or scipy.signal.cheby1 supported')

    try:
        y = signal.filtfilt(b, a, x)
    except: # Multidim array can only be processed at once for scipy >= 0.9.0
        y = []
        for data in x:
            y.append(signal.filtfilt(b, a, data))
        y = np.array(y)

    try:
        return y[:, ::q]
    except:
        return y[::q]

def normalize(signal):
    return (signal - np.mean(signal, axis=-1, keepdims=True))/np.std(signal, axis=-1, keepdims=True)
