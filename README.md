## Description

This code implements the algorithm described in

Pelegrin, B., Fernandez, F.R. Determination of efficient solutions for point-objective locational decision problems. Ann Oper Res 18, 93–102 (1989). https://doi.org/10.1007/BF02097797

Please see Pelegrin-Fernandez or [Hamilton-Jacobi scaling limits of Pareto peeling in 2D](https://arxiv.org/abs/2110.06016) for details on the algorithm. 

Our implementation also includes a fix to Step 4 of Pelegrin-Fernandez. 

## Requirements
This code requires Julia 1.5.0-1.7.0. It might work with newer versions but has only been tested for versions 1.5.0-1.7.0. The external libraries, Plots with backend GRPlots, and LinearAlgebra are also requried. 

## Usage 
See run_gift_wrapping.jl for sample usage of this code. The file gen_figs.jl generates some of the figures in [Hamilton-Jacobi scaling limits of Pareto peeling in 2D](https://arxiv.org/abs/2110.06016).
