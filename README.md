# Polynomial preconditioners for the Helmholtz equation based on Faber polynomials

This repository contains an implementation of the polynomial preconditioner for the solution of Helmholtz linear systems resulting from Helmholtz PDE problems described in the paper [GRSN21.  It also contains standalone implementations of the multigrid method and the shifted Laplace preconditioner for Helmholtz problems. In addition to solving the problems presented in the paper and the dissertation, the code in the repository can be adapted to solve other Helmholtz problems in 1D and 2D.  


## Installation
The software has been developed in MATLAB 2017b and tested with versions up to MATLAB 2020a. To clone this repository, navigate using the terminal to your desired location and type
`git clone https://github.com/luisgarciar/faber-helmholtz.git`

Next, add the path to `faber-helmholtz` to your MATLAB path, which can be done using the graphical interface in MATLAB (File -> Set Path -> Add with Subfolders).

## Usage
The folder `faber` contains the implementations of the  
Faber preconditioner and the folder `helmholtz` the code necessary for constructing the matrices and a multigrid solver.

The folder `new_experiments` contains a variety of numerical experiments for Helmholtz problems in 1D and 2D using the GMRES method and the Faber preconditioner, including those in the paper [GRSN21]. These can be adapted to solve more general problems

## References

[GRSN21] L. García Ramos, O. Sète, [Precondition- ing the Helmholtz equation with the shifted Laplacian and Faber polynomials,](https://etna.math.kent.edu/vol.54.2021/pp534-557.dir/pp534-557.pdf) Electronic Transactions on Numerical Analysis 54, pp. 534-557, 2021.