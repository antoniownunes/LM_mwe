# Levenberg-Marquardt algorithm for the determination of spacecraft trajectories in the Earth-Moon system: working examples

In this repository you will find working examples of the Levenberg-Marquardt (LM) algorithm for the purposes of transitioning trajectories from the Circular Restricted Three-Body Problem (CR3BP) into a high-fidelity ephemeris model (HFEM), considering the attraction forces of the Earth, Moon, and Sun. A full explanation of the LM implementation and deeper analysis of the results obtained may be found in [this paper](https://arxiv.org/abs/2510.18474).

## Running

These examples are established in the form of Jupyter Notebooks and may be run in two ways.

To run them without the need for any installation steps, you can follow [this binder link](https://mybinder.org/v2/gh/antoniownunes/LM_mwe/main), which will open a browser-based environment that automatically installs all necessary dependencies on the cloud. Then, select one of the notebooks and proceed normally. Note that building the environment for the first time may take a considerable amount of time.

Note also that some browsers have been found to have problems in establishing a connection to the python kernel. The use of Mozilla Firefox or Google Chrome is recommended. In addition, beware that the limited computational power available through this approach will make the examples provided slow to run at particular points along the code. 

If you want to run the examples locally, which will most likely significantly improve run speed, you will need to first install `tudatpy` and its required dependencies, as described [here](https://docs.tudat.space/en/latest/getting-started/installation.html), and create the `tudat-space` environment. This requires a `conda` installation. In addition, you will need to install Jupyter, if you haven't done yet:
```
conda install jupyter
````
To start, activate the tudat-space conda environment:
```
conda activate tudat-space
```
Then add the `tudat-space` environment to Jupyter:
```
python -m ipykernel install --user --name=tudat-space
```
Finally, you can create a local jupyter notebook instance:
```
jupyter notebook
```

Alternatively, if you use an IDE such as VS Code or other similar platforms, it is possible to run the notebooks directly at the program's interface through available extensions.

## Content

The examples provided are the following:
* **QPOs** : Showcase of the algorithm to determine quasi-periodics counterpart to periodic orbits from the CR3BP. *This is the starting example which should ideally be followed first*.
* **L2_L1_Transfer** : Showcase of the algorithm to determine a transfer trajectory from an L2 Halo orbit to an L1 Lyapunov orbit from the CR3BP. The full LM algorithm, with the possibility for adaptive weighting is employed. *This is a more advanced example that requires the user to be acquainted with the first example*. **It is not recommended to run this example in a browser through binder.**

Both examples may be easily adapted to change the initial CR3BP guess for additional testing.

## Distribution/Usage

Licensed under the MIT License. See [LICENSE](LICENSE) for details.
Dependencies are listed in the [environment.yml](environment.yml) file and are governed by their own licenses.

If you use this software and/or its routines in your research, please cite the most recent version of [this paper](https://arxiv.org/abs/2510.18474).
