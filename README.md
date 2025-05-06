<h1 align="center">POVME</h1>

<h4 align="center">Detect and characterize protein pockets.</h4>

<p align="center">
    <a href="https://github.com/durrantlab/POVME/actions/workflows/tests.yml">
        <img src="https://github.com/durrantlab/POVME/actions/workflows/tests.yml/badge.svg" alt="Build Status ">
    </a>
    <a href="https://codecov.io/gh/durrantlab/POVME">
        <img src="https://codecov.io/gh/durrantlab/POVME/branch/main/graph/badge.svg" alt="codecov">
    </a>
    <a href="https://github.com/durrantlab/POVME/blob/main/LICENSE.md" target="_blank">
        <img src="https://img.shields.io/github/license/durrantlab/POVME" alt="License">
    </a>
    <a href="https://github.com/durrantlab/POVME/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/durrantlab/POVME" alt="GitHub repo size">
    </a>
</p>

POVME (**PO**cket **V**olume **ME**asurer) is an open-source tool designed for detailed analysis of ligand-binding pocket shapes and volumes in proteins.
Initially developed to support structure-based drug discovery, POVME has been widely adopted for its simplicity, speed, and flexibility.

## Why binding pocket analysis?

Binding pockets are critical for understanding protein-ligand interactions, a central theme in rational drug design.
Pocket shape and volume influence ligand binding, specificity, and the dynamics of molecular recognition.
By enabling precise volume and flexibility measurements, POVME has become an essential tool for characterizing these properties in applications such as:

-   Identifying and analyzing transient binding pockets.
-   Supporting virtual screening and druggability predictions.
-   Comparing conformational ensembles from molecular dynamics simulations.

## Development

We use [pixi](https://pixi.sh/latest/) to manage Python environments and simplify the developer workflow.
Once you have [pixi](https://pixi.sh/latest/) installed, move into `POVME` directory (e.g., `cd POVME`) and install the  environment using the command

```bash
pixi install
```

Now you can activate the new virtual environment using

```sh
pixi shell -e dev
```

## Installation

### Development

To install the latest development version of `povme`, follow these steps.

Clone the repository.

```bash
git clone https://github.com/durrantlab/POVME
```

Navigate into the cloned directory.

```bash
cd POVME
```

Install the development version using `pip`.

```bash
pip install .
```

This will install `povme` along with all required dependencies into your current Python environment.

### Tagged

To install a specific version of `povme`, such as `v2.2.2`, follow these steps.

Clone the repository and check out the desired version tag.

```bash
git clone https://github.com/durrantlab/POVME
cd POVME
git checkout v2.2.2
```

Install the tagged version using `pip`.

```bash
pip install .
```

This will install the specified version of `povme` and its dependencies into your current Python environment.

## Cite

If you use POVME in your work, please cite:

1.  Durrant, J. D., Votapka, L., Sørensen, J., & Amaro, R. E. (2014). POVME 2.0: An Enhanced Tool for Determining Pocket Shape and Volume Characteristics. *J. Chem. Theory Comput. 10*(11), 5047-5056.
2.  Durrant, J. D., de Oliveira, C. A. F., & McCammon, J. A. (2011). POVME: An algorithm for measuring binding-pocket volumes. *J. Mol. Graph. Model. 29*(5), 773-776.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).
