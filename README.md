# Overview

The purpose of this code and data is to enable reproduction
and facilitate extension of the computational
results associated with Ref. [1].


# Installation

The library dependencies can be installed with [Docker](https://www.docker.com)
using the [Dockerfile](.devcontainer/Dockerfile).
One option is to install [Docker Desktop](https://www.docker.com/products/docker-desktop)
and open the folder with [Visual Studio Code](https://code.visualstudio.com).

# Reproducing the Manuscript

After installation, Figure 1 can be reproduced via

```bash
bash calc_figure1.sh
```

Figure 2, 3, and S1 can be reproduced via

```bash
bash calc_figure2_figure3_figureS1.sh
```

Finally, Figure 4 can be reproduced via

```bash
bash calc_figure4.sh
```

# Documentation

The documentation can be compiled using [Doxygen](https://www.doxygen.nl/),

```bash
doxygen Doxyfile
```

The documentation generated can be viewed by then opening
[index.html](./doc/html/index.html) in a browser.

# Citing This Work

To cite the manuscript, use Ref. [1]. To cite the software or experimental data, use Ref. [2].

## References

  1. DeJaco, R. F.; Kearsley, A. J. Understanding Fast Adsorption in Break-through-curve Measurements, *Under Review*, 2023.

  2. De Jaco, R. F. Sofware and Data Associated with "DeJaco, R. F.; Kearsley, A. J. Understanding Fast Adsorption in Breakthrough-curve Measurements." National Institute of Standards and Technology, 2023, [doi: 10.18434/mds2-3103](https://doi.org/10.18434/mds2-3103).
