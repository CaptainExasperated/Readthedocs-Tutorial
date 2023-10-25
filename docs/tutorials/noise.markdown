---
layout: page
title: Noise
permalink: /Noise/
---
<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

The HopsTrajectory object creates and manages the correlated stochastic process known as noise, creating "noise trajectories" $$z_{n,t}$$ for each of its independent baths $$n$$. These trajectories are made up of noise samples at timepoints $$0, \tau, 2\tau...T_{len}$$ on the noise time axis with $$Nt$$ time point.

<h2>Noise Initialization</h2>

To create this noise, HopsTrajectory takes in following user-specified input parameters in the noise_param input dictionary:
- "seed": defines noise trajectories for all modes
- "model": describes noise generation
- "tlen" ($$T_{len}$$): the end time
- "tau" ($$\tau$$): the time step

Other optional parameters include:
- "rand_model": method of generating random complex numbers
- "store_raw_noise": determines if the uncorrelated "raw" complex Gaussian random noise trajectory should be stored in Hops.noise["Z_UNCORRELATED"]

Note that tlen should exceed the duration of the calculation and thus be extended by the reorganization timescale, which is typically around 1000 fs. The number of sampled time points in each bath’s noise trajectory, $$N_t$$, will be $$\frac{T_{len}}{\tau}$$.

The noise also uses input parameters from the sys_param input dictionary:
- PARAM_NOISE1: list of these parameters for each correlation function mode
- L_NOISE1: matches each correlation function mode with the site-projection operator $$\hat{L}_n$$ of its associated bath
- ALPHA_NOISE1: gives the functional form of the correlation function modes, taking in parameters $$g_{j_n}$$ and $$w_{j_n}$$ to define exponential mode $$C_{j_n} = g_{j_n}e^{-w_{j_n}t/\hslash}$$.

While the user may enter any correlation function as a sum of exponentials, several helper functions are defined in "bath_corr_functions.py" to produce exact exponential decompositions of common correlation functions. Options include an overdamped Drude-Lorentz spectral density (with optional Matsubara modes), an underdamped Drude-Lorentz spectral density (without Matsubara modes), and a Drude-Lorentz-like overdamped spectral density (with optional Matsubara modes).

Both dictionaries will be used when creating the HopsTrajectory object. In the quickstart tutorial, they are defined here.

```
sys_param = {
    "HAMILTONIAN": np.array(hs, dtype=np.complex128),  
    "GW_SYSBATH": gw_sysbath,  
    "L_HIER": lop_list,  
    "L_NOISE1": lop_list,
    "ALPHA_NOISE1": bcf_exp,
    "PARAM_NOISE1": gw_sysbath,
}
noise_param = {
    "SEED": trajectory_index,
    "MODEL": "FFT_FILTER",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

<h2>Models and Seeding</h2>
While the noise is a stochastic process, it is necessary to be able to accurately reproduce the noise of a given trajectory for the purposes of reproducible data. As such, there are three general cases of noise trajectory seeding used to generate reproducible  noise trajectories, defined by the "MODEL" parameter of the "noise_param" dictionary. The role of the "SEED" parameter will vary based on the value of "MODEL". In this section, we will describe these roles along with a model of what the code should look like.

<h3>MODEL=ZERO</h3>
The noise trajectory is automatically set to 0 at all time points, rendering the SEED parameter irrelevant. 

Example:
```
noise_param = {
    "SEED": None,
    "MODEL": "ZERO",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

<h3>MODEL=PRE-CALCULATED</h3>
The noise trajectory is directly determined by the "SEED" parameter of the noise_param input dictionary. 

- The seed may be a two-dimensional iterable consisting of "N_L2" rows and "Nt" columns of complex numbers, where each row represents the noise trajectory of a given bath, $$z_{n,t}$$. In this example, the values are generated randomly.

Example: 
```
N_L2 = 5
Nt = 3

seed = np.array([[complex(random.random(), random.random()) for _ in range(Nt)] for _ in range(N_L2)])

noise_param = {
    "SEED": seed,
    "MODEL": "PRE_CALCULATED",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

- The seed may be a string that is a valid address for a ".npy" file. HopsNoise then loads the file at the specified address, which must be a two-dimensional iterable with "N_L2" rows and "Nt" columns of complex numbers, such that each row is the noise trajectory of a given bath, $$z_{n,t}$$. In this example, we are assuming that the variable "seed" already contains that iterable and we need to store it in a .npy file.

Example:
```
np.save("noise_seed.npy", seed)

noise_param = {
    "SEED": "/path/to/noise_seed.npy",
    "MODEL": "PRE_CALCULATED",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

<h3>MODEL=FFT_FILTER</h3>

The HopsTrajectory object generates a sample of the stochastic noise on its own, specifically a stochastic Gaussian white noise of an appropriate length for each independent bath, $$Y_{n,t}$$ where the 'uncorrelated noise' is given by $$4Nt − 4$$ uniform real Gaussian random variables, used to produce $$2Nt − 2$$ uniform complex Gaussian random variables. 

- The seed may be an integer index. This is used to seed a np.random.RandomState object which then produces the uncorrelated noise.

Example:
```
noise_param = {
    "SEED": 12345,
    "MODEL": "FFT_FILTER",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

- The seed may be a two-dimensional iterable consisting of "N_L2" rows and $$2Nt − 2$$ columns of complex numbers, where each row represents the complex uncorrelated noise of a given bath, $$Y_{n,t}$$. In this example, the values are generated randomly.

Example:
```
N_L2 = 5
Nt = 3

seed = np.array([[complex(random.random(), random.random()) for _ in range(2Nt − 2)] for _ in range(N_L2)])

noise_param = {
    "SEED": None,
    "MODEL": "FFT_FILTER",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

- The seed may be a string that is a valid address for a ".npy" file. HopsNoise then loads the file at the specified address, which must be a two-dimensional iterable with "N_L2" rows and $$2Nt − 2$$ columns of complex numbers, such that each row is the uncorrelated noise trajectory of a given bath, $$Y_{n,t}$$. This will be similar to the pre-calculated model example, where seed is a variable containing that iterable and it is converted into a .npy file.

Example:
```
np.save("noise_seed.npy", seed)

noise_param = {
    "SEED": None,
    "MODEL": "FFT_FILTER",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

- The seed may be a NoneType object. Similarly to the integer case, a np.random.RandomState produces the uncorrelated noise, but in this case the np.random.RandomState instance is not seeded and the noise is not reproducible.

Example:
```
noise_param = {
    "SEED": None,
    "MODEL": "FFT_FILTER",
    "TLEN": 500.0,
    "TAU": 1.0,
}
```

With the uncorrelated noise prepared, a circulant embedding process generates correlated noise, which is detailed in the supplementary information of this paper: [https://pubs.acs.org/doi/10.1021/acs.jpclett.3c00086](https://pubs.acs.org/doi/10.1021/acs.jpclett.3c00086)
