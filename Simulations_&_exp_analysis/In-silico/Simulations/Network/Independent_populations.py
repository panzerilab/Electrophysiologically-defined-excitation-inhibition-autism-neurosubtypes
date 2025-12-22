import time
import os
import statistics as sta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.stats
import scipy.signal
import pandas as pd
from numpy import trapz
from scipy import integrate
from scipy.interpolate import interp1d
from fooof import FOOOF, FOOOFGroup
import nest
from scipy.signal import welch
import nest.raster_plot
from scipy.stats import spearmanr
from fooof.plts.spectra import plot_spectra
from fooof.plts.annotate import plot_annotated_model
from fooof.plts.aperiodic import plot_aperiodic_params
from fooof.sim.params import Stepper, param_iter
from fooof.sim import gen_power_spectrum, gen_group_power_spectra
from fooof.utils.params import compute_time_constant, compute_knee_frequency

def power_spectrum(signal, sim_time, sim_step):
    return scipy.signal.welch(signal, fs=1000./sim_step, nperseg=(sim_time/(8.0 * sim_step)))
    
def remove_mean(signal):
    mean_value = sum(signal) / len(signal)
    signal_without_mean = [x - mean_value for x in signal]
    return signal_without_mean

def compute_spectral_slope(signal, fs):
    # signals: Matrix of input signals (each signal is a column)
    # fs: Sampling frequency
    # f_range: Frequency range for analysis
    f, S = welch(signal, fs=fs)
    psd = S
    fx = f
    idx = (fx > 30) & (fx < 70)
    fx = fx[idx]
    psd = psd[idx]
    # Take the logarithm of frequency and PSD
    log_f = np.log10(fx)
    log_psd = np.log10(psd)

    # Fit a straight line to the log-log plot using polyfit
    p = np.polyfit(log_f, log_psd, 1)

    # The slope of the line corresponds to the spectral slope
    spectral_slope = p[0]

    return spectral_slope


nest.Install('mymodule3')


NE = 1
NI = 1
simtime = 5000
epsilon = 0.2

count=0

#scaling =  np.linspace(1., 3., 15)
#scaling =  np.linspace(0.1, 1.2, 15)
scaling =  np.linspace(0.5, 3., 5)

#exc_exc_rec_values = np.linspace(0.1, 0.8, 10)

n_trials = scaling.size
interval=int(simtime)
LFP_signal = np.zeros(((n_trials)*15,interval-1))


param = []
slope = []
slope_knee = []
knee = []
offset = []
PSD_plot=[]
AMPA_plot=[]
GABA_plot =[]

for s in scaling:
    for trial in range(15):
    
            nest.ResetKernel()
       
            
            #nest.local_num_threads = 60
            
            np.random.seed(int(time.time()))
            nest.rng_seed= int(time.time())
                            
            simstep = 1.0
            fs = simtime
            alpha_e = 0.

            # Firing rate of the populations
            nu_ex = 2.0 # E
            nu_in = 5.0 # I

            # Uncomment for having same global FR but different proportions E/I
#           nu_in = 2.6/((2/5)*0.8*s + 0.2)
#           nu_ex = s*(2/5)*nu_in

            # Synaptic conductance values
            exc_exc_rec = s*1
            inh_exc_rec = -6
            
            exc_exc_recurrent = exc_exc_rec
            inh_exc_recurrent = inh_exc_rec
            
            delay = 1.0

            excitatory_cell_params_out = {
                "V_th": -2.0,
                "V_reset": -59.0,
                "t_ref": 2.0,
                "g_L": 25.0,
                "C_m": 500.0,
                "E_ex": 0.0,
                "E_in": -80.0,
                "E_L": -70.0,
                "tau_rise_AMPA": 0.1,
                "tau_decay_AMPA": 2.0,
                "tau_rise_GABA_A": 0.5,
                "tau_decay_GABA_A": 10.0,
                "tau_m": 20.0,
                "I_e": 0.0
            }

            nest.CopyModel("iaf_bw_2003_NMDA", "exc_cell_out", excitatory_cell_params_out)
            nest.CopyModel("inhomogeneous_poisson_generator", "cortical_input")

            nest.CopyModel("static_synapse", "exc_out", {"weight": exc_exc_recurrent, "delay": delay})
            nest.CopyModel("static_synapse", "inh_out", {"weight": inh_exc_recurrent, "delay": delay})

            output_cell = nest.Create("exc_cell_out", 1)
            pr_ex = nest.Create("parrot_neuron", 4000)
            pr_in = nest.Create("parrot_neuron", 1000)
            poisson_ex = nest.Create("poisson_generator", params={'rate': nu_ex})
            poisson_in = nest.Create("poisson_generator", params={'rate': nu_in})

            conn_params_out = {'rule': 'all_to_all'}
            multimeter = nest.Create('multimeter')
            multimeter.set(record_from=["V_m", "g_ex", "g_in", "g_NMDA"])
            spikerecorder_ex = nest.Create("spike_recorder")
            spikerecorder_in = nest.Create("spike_recorder")
            spikerecorder_out = nest.Create("spike_recorder")

            nest.Connect(poisson_ex, pr_ex)
            nest.Connect(poisson_in, pr_in)
            nest.Connect(pr_ex, output_cell, conn_params_out, syn_spec="exc_out")
            nest.Connect(pr_in, output_cell, conn_params_out, syn_spec="inh_out")
            nest.Connect(multimeter, output_cell)
            nest.Connect(pr_ex, spikerecorder_ex)
            nest.Connect(pr_in, spikerecorder_in)
            nest.Connect(output_cell, spikerecorder_out)

            nest.Simulate(simtime)

            m_out = multimeter.get('events')
            g_ex_out = m_out['g_ex']
            g_in_out = m_out['g_in']
            g_NMDA_out = m_out['g_NMDA']
            V_mm_out = m_out['V_m']

            AMPA = (g_ex_out * 65)
            GABA = (g_in_out * 15)
            
            f_AMPA, psd_AMPA = power_spectrum(AMPA, fs, 1)
            f_GABA, psd_GABA = power_spectrum(GABA, fs, 1)
                    
            LFP_trace = np.abs(GABA) + np.abs(AMPA)

            f_LFP, psd_LFP = power_spectrum(LFP_trace, fs, 1)

            PSD_plot.append(psd_LFP)
            AMPA_plot.append(psd_AMPA)
            GABA_plot.append(psd_GABA)
            

            sl = compute_spectral_slope(LFP_trace, 1000)
            slope.append(sl)
            
            param.append(exc_exc_rec/inh_exc_rec)
            #param.append(nu_ex/nu_in)
            #param.append(0.8*nu_ex+0.2*nu_in)
            
            LFP_signal [count] = LFP_trace
            count+=1;

# Convert lists to numpy arrays for easier manipulation
param = np.array(param)
slope = np.array(slope)
knee = np.array(knee)
offset = np.array(offset)

# Scatter plot for Slope
fig, axs = plt.subplots(1, n_trials, figsize=(15, 5))
fig.suptitle('Slope Analysis')

for i in range(n_trials):
    axs[i].loglog(f_LFP, AMPA_plot[i], label='AMPA')
    axs[i].loglog(f_LFP, GABA_plot[i], label='GABA')
    axs[i].loglog(f_LFP, PSD_plot[i], label='PSD')
    axs[i].set_xlabel('Excitatory Rate')
    axs[i].set_ylabel('Power Spectral Density')
    axs[i].legend()

plt.show()


out_dir = os.path.join("..", "Simulation results", "Independent_populations")
os.makedirs(out_dir, exist_ok=True)


# Concatenate arrays for saving
results = np.column_stack((param, slope))

np.savetxt(
    os.path.join(out_dir, "results_g.csv"),
    results,
    delimiter=',',
    header='Param, Slope',
    comments=''
)

data = pd.DataFrame(LFP_signal)
data.to_csv(
    os.path.join(out_dir, "LFP_g.csv"),
    index=False
)
