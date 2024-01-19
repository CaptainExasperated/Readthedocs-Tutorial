---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
---

<img src="https://static.wixstatic.com/media/8848f1_c606880a315245bdb81af81017dd1cf2~mv2.png/v1/fill/w_298,h_260,al_c,q_85,usm_4.00_1.00_0.00,enc_auto/MesoHOPS_Logo-01.png" alt="drawing" width="150"/>

MesoHOPS is an open-source package that uses a formally exact trajectory-based approach to solve for the time evolution of open quantum systems coupled to non-Markovian thermal environments. MesoHOPS extends the Hierarchy Of Pure States (HOPS) formalism to the adaptive Hierarchy of Pure States (adHOPS), which leverages the dynamic localization imposed by system-environment interactions to construct an adaptive basis by identifying underlying sparsity in the equation-of-motion and constructing a time-evolving reduced basis. This reduced basis allows for tractable calculations in large systems and achieves a size-invariant scaling, unique among formally exact methods, in aggregates larger than the delocalization extent of excitations. The MesoHOPS library also offers a low-temperature correction and effective integration of the noise that reduce the difficulty of simulations with ultrafast vibrational relaxation. 


<h2> About Us </h2>

Led by Dr. Doran I. G. B. Raccah, the MesoScience Lab combines theoretical chemistry, biology, physics, and applied mathematics to elucidate the principles governing charge and energy transport in molecular materials. In particular, we are focusing on developing new theoretical and computational tools to simulate photosynthetic membranes and molecular semiconductors.

[Learn more on our group website!](https://www.mesosciencelab.com/)

[Get MesoHops on our GitHub!](https://github.com/MesoscienceLab/mesohops)

<h2> Installation </h2>

MesoHops requires Python >= 3.9 and is dependent on the following packages:
* [NumPy](https://numpy.org)
* [SciPy](https://scipy.org)
* [Numba](https://numba.readthedocs.io/en/stable/#)
* [pytest](https://docs.pytest.org/en/7.4.x/)
* [pytest-level](https://pypi.org/project/pytest-level/)

(Note: all packages listed will automatically be installed when installing MesoHops)

To install the MesoHops library, first use a terminal to activate the environment of your choice and enter a new directory where the library will be stored on your computer. Next, type these lines:
```
git clone https://github.com/MesoscienceLab/mesohops.git
cd mesohops
python3 -m pip install . 
``` 

Please see our Quickstart page to run your first calculation.

<h2> Citations </h2>
D. Suess, A. Eisfeld, and W. T. Strunz, [Hierarchy of Stochastic Pure States for Open Quantum System Dynamics](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.113.150403) (2014)

L. Varvelo, J. K. Lynd, and D. I. G. Bennett, [Formally exact simulations of mesoscale exciton dynamics in molecular materials](https://doi.org/10.1039/D1SC01448J) (2021)

