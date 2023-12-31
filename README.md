# Overview

The purpose of this code and data is to enable reproduction (see [Reproducing the manuscript](#reproducing-the-manuscript))
and facilitate extension of the computational
results associated with Ref. [1]
(see [Citing This Work](#citing-this-work)).


# Installation of dependencies

The container used in this work can be obtained using

```bash
docker pull rfd1/understanding-btcs:final
```

This requires [Docker](https://www.docker.com).

For users without previous Docker experience, we recommend the following steps

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop).

2. Install [Visual Studio Code](https://code.visualstudio.com).

3. Start the Docker Desktop application.

4. Open Visual Studio Code and install the *Dev Containers extension*.

5. With the Dev Containers extension installed, you will see the (green) new status bar at the bottom left 

![](doc/remote-status-bar.png)

6. Click on the green status bar and then open this folder in a container.

# Documentation
A pdf version of the manual can be found at [doc/manual.pdf](doc/manual.pdf).

# Citing This Work

To cite the manuscript, use Ref. [1]. To cite the software or data generated, use Ref. [2].

## References

  1. DeJaco, R. F.; Kearsley, A. J. Understanding Fast Adsorption in Single-Solute Breakthrough Curves, *Communications in Nonlinear Science and Numerical Simulation*, 2023, [doi: 10.1016/j.cnsns.2023.107794](https://doi.org/10.1016/j.cnsns.2023.107794).

  2. De Jaco, R. F. Sofware and Data Associated with "DeJaco, R. F.; Kearsley, A. J. Understanding Fast Adsorption in Breakthrough-curve Measurements." National Institute of Standards and Technology, 2023, [doi: 10.18434/mds2-3103](https://doi.org/10.18434/mds2-3103).

# Reproducing the Manuscript

## [Fig. S1](out/figureS1.png)
Performing the calculations via
```bash
bash calc_figureS1.sh
```
generates the output files

*  [spatial.csv](out/spatial.csv)
*  [temporal1.csv](out/temporal1.csv)
*  [temporal2.csv](out/temporal2.csv)
*  [temporal3.csv](out/temporal3.csv).

Plotting the results via
```bash
python3 plot_figureS1.py
```
generates the figure [out/figureS1.png](out/figureS1.png).

## [Fig. S2](out/figureS2.png)
Performing the calculations via
```bash
bash calc_figureS2.sh
```
generates the output files

*  [spatial-refinement-kappa=-0.8.csv](out/spatial-refinement-theta=0.5-kappa=-0.8-e=1.csv)
*  [spatial-refinement-kappa=-0.4.csv](out/spatial-refinement-theta=0.5-kappa=-0.4-e=1.csv)
*  [spatial-refinement-kappa=0.csv](out/spatial-refinement-theta=0.5-kappa=0-e=1.csv)
*  [spatial-refinement-kappa=1.csv](out/spatial-refinement-theta=0.5-kappa=1-e=1.csv)
*  [spatial-refinement-kappa=8.csv](out/spatial-refinement-theta=0.5-kappa=8-e=1.csv)
*  [spatial-refinement-kappa=64.csv](out/spatial-refinement-theta=0.5-kappa=64-e=1.csv)
*  [spatial-refinement-kappa=256.csv](out/spatial-refinement-theta=0.5-kappa=256-e=1.csv)
*  [temporal-refinement-kappa=-0.8.csv](out/temporal-refinement-theta=0.5-kappa=-0.8-e=1.csv)
*  [temporal-refinement-kappa=-0.4.csv](out/temporal-refinement-theta=0.5-kappa=-0.4-e=1.csv)
*  [temporal-refinement-kappa=0.csv](out/temporal-refinement-theta=0.5-kappa=0-e=1.csv)
*  [temporal-refinement-kappa=1.csv](out/temporal-refinement-theta=0.5-kappa=1-e=1.csv) 
*  [temporal-refinement-kappa=8.csv](out/temporal-refinement-theta=0.5-kappa=8-e=1.csv)
*  [temporal-refinement-kappa=64.csv](out/temporal-refinement-theta=0.5-kappa=64-e=1.csv)
*  [temporal-refinement-kappa=256.csv](out/temporal-refinement-theta=0.5-kappa=256-e=1.csv)

Plotting the results via
```bash
python3 plot_figureS2.py
```
uses these output files to generate the figure [out/figureS2.png](out/figureS2.png).

## [Fig. 3](out/figure3.png), [Fig. 4](out/figure4.png), and [Fig. S3](out/figureS3.png)

Performing the calculations via
```bash
bash calc_figure3_figure4_figureS3.sh
```
generates the following output files:
*  [rarefaction-0.02.dat](out/rarefaction-0.02.dat) 
*  [rarefaction-0.04.dat](out/rarefaction-0.04.dat) 
*  [rarefaction-0.08.dat](out/rarefaction-0.08.dat) 
*  [rarefaction-0.16.dat](out/rarefaction-0.16.dat) 
*  [rarefaction-0.32.dat](out/rarefaction-0.32.dat) 
*  [rarefaction-0.64.dat](out/rarefaction-0.64.dat) 
*  [shock-0.02.dat](out/shock-0.02.dat) 
*  [shock-0.04.dat](out/shock-0.04.dat) 
*  [shock-0.08.dat](out/shock-0.08.dat) 
*  [shock-0.16.dat](out/shock-0.16.dat) 
*  [shock-0.32.dat](out/shock-0.32.dat) 
*  [shock-0.64.dat](out/shock-0.64.dat) 
*  [rarefaction--0.1-0.02.dat](out/rarefaction--0.1-0.02.dat) 
*  [rarefaction--0.2-0.02.dat](out/rarefaction--0.2-0.02.dat) 
*  [rarefaction--0.4-0.02.dat](out/rarefaction--0.4-0.02.dat) 
*  [rarefaction--0.8-0.02.dat](out/rarefaction--0.8-0.02.dat) 
*  [rarefaction--0.1-0.04.dat](out/rarefaction--0.1-0.04.dat) 
*  [rarefaction--0.2-0.04.dat](out/rarefaction--0.2-0.04.dat) 
*  [rarefaction--0.4-0.04.dat](out/rarefaction--0.4-0.04.dat) 
*  [rarefaction--0.8-0.04.dat](out/rarefaction--0.8-0.04.dat) 
*  [rarefaction--0.1-0.08.dat](out/rarefaction--0.1-0.08.dat) 
*  [rarefaction--0.2-0.08.dat](out/rarefaction--0.2-0.08.dat) 
*  [rarefaction--0.4-0.08.dat](out/rarefaction--0.4-0.08.dat) 
*  [rarefaction--0.8-0.08.dat](out/rarefaction--0.8-0.08.dat) 
*  [rarefaction--0.1-0.16.dat](out/rarefaction--0.1-0.16.dat) 
*  [rarefaction--0.2-0.16.dat](out/rarefaction--0.2-0.16.dat) 
*  [rarefaction--0.4-0.16.dat](out/rarefaction--0.4-0.16.dat) 
*  [rarefaction--0.8-0.16.dat](out/rarefaction--0.8-0.16.dat) 
*  [rarefaction--0.1-0.32.dat](out/rarefaction--0.1-0.32.dat) 
*  [rarefaction--0.2-0.32.dat](out/rarefaction--0.2-0.32.dat) 
*  [rarefaction--0.4-0.32.dat](out/rarefaction--0.4-0.32.dat) 
*  [rarefaction--0.8-0.32.dat](out/rarefaction--0.8-0.32.dat) 
*  [rarefaction--0.1-0.64.dat](out/rarefaction--0.1-0.64.dat) 
*  [rarefaction--0.2-0.64.dat](out/rarefaction--0.2-0.64.dat) 
*  [rarefaction--0.4-0.64.dat](out/rarefaction--0.4-0.64.dat) 
*  [rarefaction--0.8-0.64.dat](out/rarefaction--0.8-0.64.dat) 

These figures can be generated by
```bash
python3 plot_figure3_figureS3.py && python3 plot_figure4.py
```
and their output ```.png``` files found in the ```out/``` folder as before.

## [Fig. 5](out/figure5.png)
Performing the calculations via
```bash
bash calc_figure5.sh
```
generates the output files that can be found in folders

* [movie-1.0](out/movie-1.0/)
* [movie-6.0](out/movie-6.0/)
* [movie-36.0](out/movie-36.0/)

The figure is generated via
```bash
python3 plot_figure5.py
```
