## Beam Induced Neutrino events at a Muon Collider

Monte-Carlo for Beam Induced Neutrinos (BIN) in a Muon Collider (MC).

This code uses lattices provided by the IMCC collaboration. To implement them, we create smooth interpolators of the lattice. These are large files, so we suggest you create them by runnnig `create_beam_optics_files.ipynb` (should take no longer than a few minutes).

Several examples of how to generate events and the plots used for the paper are in `BIN_characterization.ipynb` and `BIN_forward_flux.ipynb`.

The LaTeX tables of total BIN interaction rates are generated in `BIN_characterization.ipynb`.

Some work on neutrino-electron scattering is in `BIN_electron_scattering_rates.ipynb`.

## Installation

Requirements:
* DarkNews
* vegas
* pandas
* numpy
* numba
* scipy

## Citation 

If you use this code, please cite:
```
@article{Bojorquez-Lopez:2024bsr,
    author = {Bojorquez-Lopez, Luc and Hostert, Matheus and Arg\"uelles, Carlos A. and Liu, Zhen},
    title = "{The Neutrino Slice at Muon Colliders}",
    eprint = "2412.14115",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "UMN-TH-4409/24",
    month = "12",
    year = "2024"
}
```