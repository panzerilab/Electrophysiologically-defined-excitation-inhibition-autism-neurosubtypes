# Sim_ephys_analysis_Bertelsen

This repository contains all the code required to **reproduce the analyses and figures** presented in:

**Bertelsen et al., 2023**  
*Electrophysiologically-defined excitation–inhibition autism neurosubtypes*

The codebase includes both **in-silico network simulations** and **in-vivo electrophysiological data analysis**, covering feature extraction, statistical analysis, and figure generation.

---

## Repository Structure

```

├── in-silico/
│   ├── Simulations/
│   └── Analysis & plot/
│
├── in-vivo/
│   ├── Scripts/
│   ├── DATA/
│   ├── results/
│   └── Figures/
│
└── functions/

```

---

## 1. In-silico Analysis

The `in-silico` folder contains everything needed to generate network simulations and analyze simulated LFP activity corresponding to:

- **Figure 1**
- **Supplementary Figures 1, 2, and 3**
- Associated statistical tables

---

### 1.1 Simulations

Path:
```

in-silico/Simulations/

````

This folder contains:
- Network model definitions
- Simulation scripts
- Output directories for simulation results

#### Network Model

- Simulations are based on **NEST Simulator**
- Custom neuron model:
  - Conductance-based LIAF neuron
  - AMPA and GABA synapses
  - Synaptic conductances modeled as a **difference of two exponentials**
- Model implementation follows:
  - *Cavallari et al., 2015*

---

### 1.2 Building the Neuron Model (NEST Extension)

The custom neuron model must be compiled before running simulations.  
The procedure follows the NEST tutorial *“Writing an extension module”*.

#### Step 1 — Set the NEST installation directory

```bash
export NEST_INSTALL_DIR=/Users/gabriele/NEST/ins
````

#### Step 2 — Create and enter the build directory

```bash
cd neuron_model
mkdir build
cd build
```

#### Step 3 — Configure the extension module

If `nest-config` is not in your `PATH`, specify it explicitly:

```bash
cmake -Dwith-nest=${NEST_INSTALL_DIR}/bin/nest-config ..
```

#### Step 4 — Compile and install

```bash
make
make install
```

You may need to update the library path:

```bash
export LD_LIBRARY_PATH=${NEST_INSTALL_DIR}/lib/python2.7/site-packages/nest:$LD_LIBRARY_PATH
```

---

### 1.3 Running Simulations

Once the neuron model is installed, simulations can be run from:

```
in-silico/Simulations/Network/
```

* Each script runs a different set of simulations used in the paper
* The number of CPU cores can be specified in the main script to speed up simulations
* Simulation outputs are saved in:

```
in-silico/Simulations/Simulation results/
```

---

### 1.4 Simulation Analysis and Figures

Path:

```
in-silico/Analysis & plot/
```

The analysis pipeline is divided into:

1. Feature extraction
2. Statistical analysis and figure generation

#### Feature Extraction

Run the scripts in:

```
in-silico/Analysis & plot/feature_extraction/
```

This step extracts:

* Periodic and aperiodic components of simulated LFPs
* All associated parameters and statistics

#### Figure Generation

To reproduce:

* **Figure 1**
* **Supplementary Figures 1, 2, and 3**

Run the scripts named after each figure.
Individual figure panels are saved in:

```
in-silico/Figures/
```

---

## 2. In-vivo Analysis

The `in-vivo` folder contains all code required to analyze experimental electrophysiological data.

---

### 2.1 Preprocessing Raw Data

Raw data must first be processed and split into trials (chunks) using external code available at:

**[LINK TO RAW DATA PROCESSING REPOSITORY]**

---

### 2.2 Data Organization

After preprocessing, place the chunked data in:

```
in-vivo/DATA/
```

Organized as follows:

```
DATA/
├── LFP/
│   ├── manipulation_1/
│   ├── manipulation_2/
│   └── ...
└── MUA/
    ├── manipulation_1/
    ├── manipulation_2/
    └── ...
```

Each subfolder corresponds to a different **chemogenetic manipulation**.

---

### 2.3 Requirements (MATLAB)

The following MATLAB toolboxes are required:

* `chronux_2_12`
* `fooof_mat-main`
* `gramm-master`
* `nonfractal-master`
* `wmtsa-matlab-0.2.6`

Custom helper functions are located in:

```
/functions
```

> The path to `/functions` is automatically added at the beginning of each script.

---

### 2.4 Feature Extraction (In-vivo)

To extract features from LFP and MUA signals, run:

```
in-vivo/Scripts/Feature_extraction.m
```

This script:

* Extracts periodic and aperiodic components
* Computes all relevant statistics
* Saves results to:

```
in-vivo/results/
```

---

### 2.5 In-vivo Figures and Statistics

After feature extraction, run the scripts corresponding to each figure to reproduce:

* **Figure 2**
* **Supplementary Figures 4 and 5**

Outputs:

* Figures are saved in:

  ```
  in-vivo/Figures/
  ```
* Statistical results are saved in the same directory



## Contact

For questions or issues related to the code, please contact:

**Gabriele Mancini**
gabriele.mancini@iit.it


# Bertelsen_et_al
