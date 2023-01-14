Research project, modeling astrocyte calcium dynamics driven by glutamate stimulations.

# Folders

* **older MATLAB code**: original MATLAB code from Greg prior to this project that was translated to Python. This includes code that was used to simulate GPCR model as well as classify responses from resulting calcium signals
* **figures**: collections of figures generated from Python
* **resources**: papers that this project originated from and notes from Alla of things that we originally wanted to tackle
* **xpp**: xpp files to run bifurcation and phase plane analysis from
    * G_lambda_ip3_ca_bifurcation.ode: the most up to date version of the full GPCR with lambda -> IP3 -> Ca2+ model
    * gpcr.ode: the GPCR model alone, from which we can model G* (activated G)
    * ip3_ca_bifurcation.ode: IP3 -> Ca2+ part of model. We can model the Ca2+ dynamics from fixed IP3
    * ip3_ca_bifurcation_<delay/fixed/reduced/etc.>: variations on the IP3 -> Ca2+ model to explore effects on oscillation delays
* **presentations**: progress updates, reports, and poster documents that were created
* **Python**: main Python folder including code to run simulations and generate figures
    * See "jupypter_notebook_descriptions.md" for descriptions on files in this folder

