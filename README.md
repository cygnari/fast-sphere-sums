# fast-sphere-sums
Spherical tree code for fast computation of spherical convolutions

To compile the code, first clone the repository. `cd` into `fast-sphere-sums` and `mkdir build`. Make sure that you have Cmake, a C++ compiler and a usable MPI library, and then `cd build` and `cmake ..`. After that, `make -j`. To run the code, go into `./build/executables/` where you'll find a bunch of different executables. Modify the `namelist.txt` and then `mpirun -np XX desired_computation`. For example, to solve the Poisson equation on an icosahedral grid with 655362 points and a spherical harmonic right hand side, set `point_count=655362` and `grid=icos_grid` and `initial_condition=SH43` in the namelist. Then run `inverse_laplacian`. 
