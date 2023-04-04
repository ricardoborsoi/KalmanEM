# Kalman filtering and expectation maximization for multitemporal spectral unmixing    #

This package contains the implementation of multi temporal unmixing algorithm proposed in the paper [1].

The code implements a multi temporal hyperspectral unmixing (MTHU) algorithm using physically motivated parametric endmember representations to account for temporal end member variability. By representing the multitemporal mixing process using a state-space formulation, it exploits Bayesian filtering to estimate the endmember variability coefficients, and an efficient implementation of the expectation–maximization (EM) algorithm is used to estimate the abundances and other model parameters.


## Usage

The algorithm is implemented in Matlab, you can run it by running the following functions:  
-  `main_ex1.m` : example with synthetic data, it is the Data Sequence 1 in reference: https://arxiv.org/pdf/2303.10566. 
-  `main_ex2.m` : example with synthetic data, it is the Data Sequence 2 in reference: https://arxiv.org/pdf/2303.10566. 
-  `main_ex_Tahoe.m` : example with real data from the Lake Tahoe image, there real data example from reference https://arxiv.org/pdf/2303.10566. 


### Downloading the datasets

The datasets are available in Zenodo [here](https://zenodo.org/record/7796598#.ZCt8VC8iuEc), https://doi.org/10.5281/zenodo.7796597. Please download them and include the correct paths to them in the .json files to be able to run the examples.




## REFERENCE:
If you use this software please cite the following in any resulting publication:

    [1] [Kalman filtering and expectation maximization for multitemporal spectral unmixing](https://ieeexplore.ieee.org/abstract/document/9219120)
        R.A. Borsoi, T. Imbiriba, P Closas, J.C.M. Bermudez, C. Richard.
        IEEE Geoscience and Remote Sensing Letters, 2020.

