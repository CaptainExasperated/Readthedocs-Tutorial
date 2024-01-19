---
layout: page
title: Storage
permalink: /Storage/
---
<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

The HopsTrajectory object saves information about the trajectory after each time step. At each time point, it takes in input parameters consisting of the full wave function ($$\psi_{t}$$), the auxiliary basis, the state basis, the current time point $$t$$, and the set of noise memory drift terms $$\xi_{n,t}$$, then saves a user-selected set of output information for later analysis.

<h2>Saving and Accessing Information</h2>

By default, the HopsTrajectory object stores the system wave function $$\vert \psi^{\vec{0}}_{t} \rangle$$ and time $$t$$ when the calculation is not adaptive. When the calculation is adaptive (see [Adaptivity](/Readthedocs-Tutorial/Adaptivity/) for more details), the HopsTrajectory object also stores the outputs of list of indexing vectors of each member of $$A_t$$, the list of the integer indices of each state in $$S_t$$,  number of states in $$S_t$$, and number of auxiliary wave functions in $$A$$.

The information is stored in a list dictionary which can be accessed through "hopsobj.storage" (where hopsobj is the name of the HopsTrajectory created). The information can be accessed by referencing the information’s key name.

We already do this in the base tutorial:

```
hops.propagate(t_max, t_step)
psi_traj = hops.storage['psi_traj']
```

The HopsTrajectory object automatically saves the system wave function when hops.propagate is called. We then access the system wave function at all time points through the storage dictionary (hops.storage) via its key name ('psi_traj').

<h2>Specifying Saved Information</h2>

We can specify which information we choose to save by adding keys to a new dictionary (storage_param) when creating the HopsTrajectory object. This dictionary takes in the key name and the function used to store that information. There are eight built-in and preconstructed saving functions in mesoHOPS that can be found and modified in the file storage_functions.py. Note that the key name for each saving function is the saving function's name excluding the inital "save_".  

| Saving Function    | Output |
| -------- | ------- |
| save_psi_traj  | system wave function $$\vert \psi^{\vec{0}}_{t} \rangle$$   |
| save_phi_traj | full wave function $$\psi_{t}$$    |
| save_phi_norm    | norm of full wave function $$\rVert \psi_{t} \rVert_{2}$$   |
| save_t_axis  | time $$t$$   |
| save_aux_list | list of indexing vectors of each member of $$A_t$$    |
| save_list_nhier    | number of auxiliary wave functions in $$A_t$$   |
| save_state_list  | list of the integer indices of each state in $$S_t$$   |
| save_list_nstate | number of states in $$S_t$$    |

For example, if we wanted our tutorial code to additionally save the full wave function as well as its default system wave function and time, we would first create a dictionary

```
storage_param = {
    "psi_traj": save_psi_traj,
    "t_axis": save_t_axis,
    "phi_traj": save_phi_traj
}
```

And then store it in our hops object

```
hops = HOPS(
        sys_param,
      noise_param=noise_param,
      hierarchy_param=hierarchy_param,
      eom_param=eom_param,
      integration_param=integration_param,
      storage_param=storage_param
    )
```

The "data" dictionary will then contain the "psi_traj," "t_axis" and "phi_traj" for each time point after propagate() is called. We can access the full wave function as such:

```
phi_traj = hops.storage['phi_traj']
```

Afterwards, we can optionally save this information into its own file using a function in numpy called np.save(), which takes in the name of the file and the information saved

```
np.save("full_wave_function_file", phi_traj)
```

<h2>User-Defined Saving Functions</h2>
A user may create and modify the key names and their saving functions easily to tailor the saved data to a particular analysis.

To modify the key names, we simply change the key name stored in storage_param. For, example, we can rename t_axis to time so that we store it with 

```
storage_param = {
    "psi_traj": save_psi_traj,
    "time": save_t_axis,
    "phi_traj": save_phi_traj
}
```

And access it with

```
t = hops.storage["time"]
```

When creating functions, we need to account for five parameters:
- "phi_new": the full wavefunction ($$\psi_{t}$$).
- "aux_list": the list of AuxilaryVector objects in the hierarchy ($$A_{t}$$).
- "state_list": the list of states in the state basis ($$\S_{t}$$).
- "t_new": the time t.
- "z_mem_new": the set of noise memory drift terms $$\xi_{n,t}$$ for each independent bath $$n$$.

Saving functions are constructed in the following format: 
```
def save_basis_size (aux_list, state_list, **kwargs):
"""
returns the integer size of the full Hilbert space at a given time point.
"""
    n_aux = len(aux_list)
    n_site = len(state_list)
    return n_aux*n_site
```

where the parameters that will be used to calculate the output are defined as inputs and any other input parameters are allowed but unused via the "**kwargs" parameter.

We can then store that information by modifying storage_param
```
storage_param = {
    "psi_traj": save_psi_traj,
    "time": save_t_axis,
    "basis_size": save_basis_size
}
```

Here is a modified version of the example code that creates a saving function, specifies for HopsTrajectory to store that information at each timestep, and then accesses that information through a variable.

```
import os
import numpy as np
import scipy as sp
from scipy import sparse
from pyhops.dynamics.hops_trajectory import HopsTrajectory as HOPS
from pyhops.dynamics.eom_hops_ksuper import _permute_aux_by_matrix
from pyhops.dynamics.bath_corr_functions import bcf_exp, bcf_convert_sdl_to_exp
import matplotlib.pyplot as plt


def save_basis_size (aux_list, state_list, **kwargs):
    """
    Returns the integer size of the full Hilbert space at a given time point.
    """


    n_aux = len(aux_list)
    n_site = len(state_list)
    return n_aux * n_site


wf_list = []
ntraj = 100


for trajectory_index in range(ntraj):


    nsite = 2
    hs = np.zeros([nsite, nsite])
    hs[0, 0] = 100
    hs[0, 1] = 50
    hs[1, 0] = 50
    hs[1, 1] = 0


    e_lambda = 60.0  
    gamma = 60.0  
    temp = 300.0  
    (g_0, w_0) = bcf_convert_sdl_to_exp(e_lambda, gamma, 0.0, temp)  


    loperator = np.zeros([nsite, nsite, nsite], dtype=np.float64)  
    gw_sysbath = []  
    lop_list = []  
    for i in range(nsite):  
        loperator[i, i, i] = 1.0  
        gw_sysbath.append([g_0, w_0])  
        lop_list.append(sp.sparse.coo_matrix(loperator[i]))  
        gw_sysbath.append([-1j * np.imag(g_0), 500.0])  
        lop_list.append(loperator[i])


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


    eom_param = {"EQUATION_OF_MOTION": "NORMALIZED NONLINEAR"}
    integration_param = {"INTEGRATOR": "RUNGE_KUTTA",
                         "EARLY_ADAPTIVE_INTEGRATOR": "INCH_WORM",
                         "EARLY_INTEGRATOR_STEPS": 5,
                         "INCHWORM_CAP": 5}  
    hierarchy_param = {"MAXHIER": 4}


    storage_param = {
        "psi_traj" : save_psi_traj,
        "time" : save_t_list,
        "basis_size" : save_basis_size
    }


    hops = HOPS(
        sys_param,
      noise_param=noise_param,
      hierarchy_param=hierarchy_param,
      eom_param=eom_param,
      integration_param=integration_param,
      storage_param=storage_param
    )


    psi_0 = np.array([0.0] * nsite, dtype=np.complex)
    psi_0[0] = 1.0


    psi_0 = psi_0 / np.linalg.norm(psi_0)


    t_max = 200.0  
    t_step = 4.0  
   
    hops.initialize(psi_0)
    hops.propagate(t_max, t_step)


    psi_traj = hops.storage['psi_traj']
    wf_list.append(psi_traj)
    np.save("trajectory" + str(trajectory_index), psi_traj)


    b_size = hops.storage["basis_size"]
    np.save("basis_size" + str(trajectory_index), b_size)


mean_pop = np.mean(np.abs(np.array(wf_list))**2, axis=0)
pop1 = mean_pop[:,0]
pop2 = mean_pop[:,1]


np.save("mean_pop_depth_{}_datatype_{}_traj_{}_eq".format("nparray", str(ntraj), "linear"), mean_pop)


t_axis = np.arange(0, t_max+t_step, t_step)
plt.xlabel("Time (fs)")
plt.ylabel("Population")
plt.plot(t_axis,pop1)
plt.plot(t_axis,pop2)
plt.show()
```