---
layout: page
title: The Low-Temperature Correction and Effective Integration of the Noise
permalink: /LowTempCorr/
---
[comment]: <> (This allows latex formatting in the markdown files)
<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

In many systems, the correlation function describing the thermal environment contains modes that are significantly faster than the others: these ultrafast modes lead to time-sensitive dynamics that can be difficult to converge. Thankfully, we can use the HOPS low-temperature correction in conjunction with an effective integration of the noise to approximate the effect of ultrafast modes with high accuracy and low computational cost.

## What is the Low-Temperature Correction?
In HOPS, the correlation function of each independent bath is comprised of a sum of exponentials. When an exponential term decays sufficiently quickly, the following delta function approximation becomes appropriate:

$$ge^{-\gamma t} \approx \frac{g}{\gamma}\delta(t).$$
  
We refer to modes that can be approximated as delta functions as "ultrafast." The delta function approximation turns out to significantly simplify the HOPS equation of motion. Better yet, because all delta functions decay on the same (infinitely-fast) timescale, we find that for two ultrafast modes in the same bath:

  

$$g_1e^{-\gamma_1 t} + g_2e^{-\gamma_2 t} \approx \frac{g_1}{\gamma_1}\delta(t) + \frac{g_2}{\gamma_2}\delta(t) = (\frac{g_1}{\gamma_1}\ + \frac{g_2}{\gamma_2})\delta(t).$$

  

Thus, each independent bath has a single low-temperature correction factor, given by the sum of correlation function mode magnitude (g) divided by correlation function mode decay frequency (w) for each ultrafast mode in the bath correlation function: these modes will be included in the list of correlation function modes for the noise, but not in the list of correlation function modes for the hierarchy. The proper way of introducing low-temperature-corrected ultrafast modes is displayed below:

  ```
e_lambda = 60.0
gamma = 60.0
temp = 300.0

(g_0, w_0) = bcf_convert_sdl_to_exp(e_lambda, gamma, 0.0, temp)
# Arbitrary ultrafast modes to illustrate the example.
(g_1, w_1) = (np.real(g_0)/5, w_0*5)
(g_2, w_2) = (np.real(g_0)/10, w_0*10)

# Defines the set of modes that go into the noise and hierarchy
# correlation functions of each bath. The modes that are in the 
# noise correlation function but not the hierarchy correlation
# function are treated by the low-temperature correction.
modes_noise = [g_0, w_0, g_1, w_1, g_2, w_2, -1j * np.imag(g_0), 500.0]
modes_hier = [g_0, w_0, -1j * np.imag(g_0), 500.0]
ltc_coeff = g_1/w_1 + g_2/w_2
loperator = np.zeros([nsite, nsite, nsite], dtype=np.float64)

gw_noise = []
lop_noise = []
gw_sysbath = []
lop_sysbath = []
lt_corr_factor = []
lop_lt_corr = []

for i in range(nsite):
loperator[i, i, i] = 1.0

# Assign modes of the noise correlation function to the appropriate L-operator.
gw_noise.append(modes_noise)
lop_noise.append([sp.sparse.coo_matrix(loperator[i])]*len(modes_noise))

# Assign modes of the hierarchy correlation function to the appropriate L-operator.
gw_sysbath.append(modes_hier)
lop_sysbath.append([sp.sparse.coo_matrix(loperator[i])]*len(modes_hier))

# Assign low-temperature correction facrtors to the appropriate L-operator.
lt_corr_factor.append(ltc_coeff)
lop_lt_corr.append(sp.sparse.coo_matrix(loperator[i]))

sys_param = {
"HAMILTONIAN": np.array(hs, dtype=np.complex128),
"GW_SYSBATH": gw_sysbath,
"L_HIER": lop_sysbath,
"L_NOISE1": lop_noise,
"ALPHA_NOISE1": bcf_exp,
"PARAM_NOISE1": gw_noise,
"PARAM_LT_CORR": lt_corr_factor,
"L_LT_CORR": lop_lt_corr
}
```
  
  
  

## What is the Effective Integration of the Noise?

MesoHOPS works by using the HOPS equation of motion works to find the derivative of each auxiliary wave function and then integrating with a fourth-order Runge-Kutta integrator. If a single time step of integration propagates the dynamics from time $$t$$ to time $$t + \Delta t$$, the fourth-order Runge-Kutta integrator samples the noise at times $$t$$, $$t + \frac{1}{2}\Delta t$$, and $$t + \Delta t$$.

  

In many calculations - especially when the bath correlation function has ultrafast modes - the noise has dynamics on a faster timescale than anything else in the equation of motion. Therefore, one of the limiting factors in the accuracy of a simulation is whether half the time step of integration $$\Delta t/2$$ is fine-grained enough to capture short-lived fluctuations in the noise that we refer to as ultrafast noise dynamics.

  

Because a larger time step of integration results in faster calculations, we use an "effective integration of the noise" to treat the noise with a smaller effective time step. The effective integration works by taking a moving average over all noise time points between the current time and the next half time step of integration, accounting for the points on the noise time axis that are not captured by half the time step of integration.

  

The effective noise integration may be activated in the integration_param dictionary:

```
integration_param  = {"INTEGRATOR": "RUNGE_KUTTA",
					  "EARLY_ADAPTIVE_INTEGRATOR": "INCH_WORM",
					  "EARLY_INTEGRATOR_STEPS": 5,
					  "INCHWORM_CAP": 5,
					  "EFFECTIVE_NOISE_INTEGRATION": True
}

hops  =  HOPS(sys_param,
			  noise_param=noise_param,
			  hierarchy_param=hierarchy_param,
			  eom_param=eom_param,
			  integration_param=integration_param,
			  storage_param=storage_param
)
```