## Beam Induced Neutrino events at a Muon Collider

Monte-Carlo for Beam Induced Neutrinos (BIN) in a Muon Collider (MC).

This code uses lattices provided by the IMCC collaboration. To implement them, we create smooth interpolators of the lattice, which are created by `create_beam_optics_files.ipynb`. The files are over 200 MB in size, they need to be created (should take no longer than a few minutes).

Several examples of how to generate events and the plots used for the paper are in `BIN_characterization.ipynb` and `BIN_forward_flux.ipynb`.

The LaTeX tables of total BIN interaction rates are generated in `BIN_rate_tables.ipynb`.

Some work on neutrino-electron scattering is in `BIN_electron_scattering_rates.ipynb`.

## Installation

Requirements:
* vegas
* pandas
* numpy
* DarkNews (pip install DarkNews)
* numba
* scipy

## Citation 

If you use this code, please cite:

@article{...}