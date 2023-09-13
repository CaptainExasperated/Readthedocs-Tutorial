---
layout: page
title: Convergence Testing HOPS Calculations
permalink: /ConvergenceTesting/
---
[comment]: <> (This allows latex formatting in the markdown files)
<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

When simulating excited-state dynamics with mesoHOPS, it is always necessary to ensure that your simulation is accurate. We do this by convergence testing - that is, showing that systematically increasing the accuracy of a given simulation does not meaningfully alter the results.

## Convergence Parameters

Both the accuracy and the computational expense of a given HOPS calculation are governed by a set of parameters that control the calculation itself: these "convergence parameters" do not define the open quantum system that we wish to simulate. "Tightening" any of these convergence parameters results in a calculation that is both more computationally expensive and more accurate, while "loosening" them does the opposite. These convergence parameters are:

-   kmax, the maximum hierarchy depth (larger is tighter)
    
-   dt, the time step of integration (smaller is tighter)
    
-   kmats, the number of Matsubara modes in the exponential decomposition of the correlation function (larger is tighter)
    

In an adaptive HOPS calculation, there are additional convergence parameters:

-   delta_a, the auxiliary wave function basis derivative error bound (smaller is tighter)
    
-   delta_s, the state basis derivative error bound (smaller is tighter)
    
-   ut/us, the update time/update step, which is the time/number of integration time steps between updates to the adaptive basis (smaller is tighter)
    

  

Ideally, we want to find a calculation that is "converged" - that is, tightening all of the convergence parameters does not change the dynamics - with the loosest possible convergence parameters. In order to do this, we need to complete "convergence scans."

## Convergence Scans

A convergence scan is an ensemble of HOPS trajectories (usually of 1,000 trajectories) calculated at different convergence parameters but with the same noise trajectories (since Matsubara modes affect the noise, this means the same uncorrelated noise trajectories - see the [noise tutorial](https://duckduckgo.com)). This allows us to complete the "same" simulation with tighter or looser convergence parameters and meaningfully measure error introduced by insufficiently tight convergence parameters.

  

There are two major types of convergence scans: parameter-scan convergence scans begin with a set of "over-converged" parameters (determined based on other simulations or simply by selecting absurdly tight convergence parameters) and loosen a single convergence parameter at a time. While this is useful for finding the loosest set of convergence parameters for which calculations are still converged, it is also time-consuming and computationally expensive, because the bulk of the calculations will be overly expensive and each parameter must be individually surveyed.

  

On the other hand, manifold convergence scans also begin with a set of over-converged parameters, but loosen all parameters at once. This significantly reduces the number of simulations run for a convergence scan, but because we loosen all the convergence parameters at the same time, it is unlikely that we will find the loosest set of convergence parameters for which calculations are still converged. Overall, however, manifold convergence scans are often much more convenient that parameter-scan convergence scans. It can also be convenient to do manifold convergence scans and then attempt to more tightly converge a single parameter that you suspect to be making calculations unnecessarily expensive.

  

In either case, we call the most tightly-converged ensemble the "reference ensemble." In order to measure how converged each "test ensemble" (the ensembles that are less tightly-converged) is, we need to develop a measurement of error.

## Measuring Error and Determining Convergence

There are numerous ways to measure the convergence of a calculation, but we will focus on population error in this tutorial. To calculate the population error, we first need to calculate the "population vector" of each ensemble at each time point. This is simple to do in Python.

  

```
#Where ensemble_traj is a list of HopsTrajectory.storage['psi_traj'] for each trajectory

Mean_pop = np.mean([np.abs(ensemble_traj[i])**2 for i in range(len(ensemble_traj))], axis=0)
```

This gives population vector at all times $$t$$, $$\vec{P}(t)$$, comprised of the populations of the states defining the basis of the system Hamiltonian. We can find the population error between the test and reference ensemble as the norm of the difference between population vectors $$E(t) = \|\vec{P}_{ref}(t)-\vec{P}_{test}(t)\|_1$$ where $$\vec{P}_{ref}(t)$$ and $$\vec{P}_{test}(t)$$ are the population vectors at time $$t$$ for the reference and test ensembles, respectively.

  

The user may determine how to characterize error - and which levels of error are acceptable for considering a test ensemble converged - according to their own needs. Generally, defining a threshold value for either the maximum or mean of $$E(t)$$ is a good measure of convergence: we have previously used mean error thresholds of 0.02 and 0.03 to define convergence.

  

One alternative to population error is error in the expectation value of an observable. Instead of using the population vector $$\vec{P}(t)$$, it is simple to compare the expectation value of $$\hat{O}$$ over time via taking the ensemble average of the expectation value $$\langle \hat{O}\rangle_t = \langle\psi^{\vec{0}}_t\vert\hat{O}\vert\psi^{\vec{0}}_t\rangle$$. Comparing the expectation value of the observable between a reference and test ensemble is simple enough from there.

Finally, rather than using a numerical convergence threshold that depends on some calculated error measure, the user may simply determine convergence qualitatively by plotting either the populations of some or all states or some observable and deciding by eye whether or not the data appears converged.

  
  

## Statistical Convergence

It is also important to ensure that a HOPS calculation contains sufficient trajectories to accurately represent the ensemble behavior. The number of trajectories that proves statistically significant is system-dependent, but we find that generally, 5,000-10,000 trajectories are sufficient. Note that during parameter convergence scans, ensembles can consist of 1,000 trajectories without issue because the noise trajectories of the HOPS trajectories in the ensemble are matched: the object that must be statistically converged is a scalar error measure.

  

To test statistical convergence, we use a technique known as bootstrapping. We set up a number of N-trajectory ensembles, and then, for each ensemble, randomly sample (with replacement) M ensembles of N trajectories (M is usually a very large number, such as 10,000), before taking the 95% confidence interval of either the population vector or some observable expectation value across the M ensembles. This allows us to characterize the statistical significance of a given ensemble size.

  

Bootstrapping, however, is time-consuming and involved. While it provides in-depth information about statistical convergence, simply treating the number of trajectories in an ensemble as any other convergence parameter and doing a convergence scan is a perfectly reasonable approach. Statistical convergence may also be treated qualitatively: one of the main aims of statistically-converged ensembles is the production of publication-quality figures. A qualitatively statistically-converged ensemble has data that appears smooth to the eye and is not visually different from that of a larger ensemble (e.g., one with twice the number of trajectories).

  

## Notes

The noise time step, $$\tau$$, must be strategically set before beginning a convergence scan: all tested integration time steps must be an even integer multiple of $$\tau$$ due to the way the Runge-Kutta integrator works.

  

Setting the hierarchy depth too high can lead to issues with time step, as a large hierarchy depth introduces dynamics on very fast timescales into the equation of motion. As a general rule of thumb, ensure that your smallest tested integration time step is less than or equal to $$\frac{\hbar}{2\gamma_{max}k_{max}}$$, where $$\gamma_{max}$$ is the fastest exponential mode in your correlation function.

  

If the temperature is too low, ignoring Matsubara modes entirely will lead to catastrophically unconverged results. At non-cryogenic temperatures, you rarely need more than 5 Matsubara modes, although some systems are more sensitive to Matsubara modes than others. The treatment of Matsubara modes also matters: at room temperature, Matsubara modes can generally be treated with a low-temperature correction or Markovian filter, but this is not guaranteed by any means.

  

When doing convergence scans of an adaptive calculation, note that you set update step, not update time: the update time is the parameter to converge, but it is given by updates step multiplied by the integration time step. Thus, to test integration time step in isolation, you must change update step to ensure that update time is static.

  

While the other parameters may be somewhat predictable by examining the properties of the system and behavior of the analytic correlation functions of the baths, the convergence behavior of adaptive parameters has historically proven highly unpredictable between systems. A good starting value is 1E-4 for each: this is usually converged. However, you may very well find that a delta value of as large as 0.05 is comfortably converged and that using such a delta value makes your calculations much faster.

  

Small changes in convergence parameters do not often yield major changes in calculation results. I recommend testing integration time step and update time on a log-2 scale (e.g., testing time steps 1.0 fs, 2.0 fs, 4.0 fs and 8.0 fs and update times 8 fs, 16 fs, 32 fs, and 64 fs), hierarchy depth in steps of 1, 2, 5, or 10 depending on the system (e.g., 1, 2, 3, and 4 or 10, 15, 20, 25), number of Matsubara modes on the modified log scale of 0, 1, 2, 5 (in most systems at room temperature), and the adaptive error limit parameters on a log-10 scale of something like 0.003, 0.001, 0.0003, 0.0001.