# rkma

This repository contains the MATLAB code that generates the figures in the paper "The Randomized Kaczmarz Method with Mismatched Adjoint".

##### Authors:
- Dirk Lorenz    (<d.lorenz@tu-braunschweig.de>)
- Sean Rose    (<seanrose949@gmail.com>)
- Frank Schöpfer (<frank.schoepfer@uni-oldenburg.de>

Contents
--------

##### Drivers (run these to generate figures):
    example_convergence_perturbation.m  script to generate figure 2
	example_inconsistent.m              script to generate figure 3
	example_underdetermined.m           script to generate figure 4
	example_ct.m                        script to generate figure 5
	example_optimize_probabilities.m    script to generate figure 6

##### Routines called by the drivers:
	rkma.m                      function for the randomized Kaczmarz
	sampling.m                  function to sample according to a prob-vector

##### Dependencies
`example_optimize_probabilities.m` needs `projsplx.m` from <http://ufdc.ufl.edu/IR00000353/>.

`example_ct.m` needs the AIRtool package <https://doi.org/10.1016/j.cam.2011.09.039>

To use the writeout functionality requires `matlab2tikz` from <https://github.com/matlab2tikz/matlab2tikz>
