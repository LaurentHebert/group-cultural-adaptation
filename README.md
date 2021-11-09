# LHD

Implementation of the model from Hébert-Dufresne et al. on [Source-sink cooperation dynamics constrain institutional evolution in a group-structured society](https://arxiv.org/abs/2109.08106).

### Codes

Integration of the ODEs is done in C++ using algorithms from the Gnu Scientific Library and data structures from the Boost Library.

The code tevol_source_gca.cpp takes in a series of parameters and can spit out either the final state of the system or the full time evolution.

Run the code without any argument for a list of the required parameters.

The file dyn_gca.hpp contains the actual dynamical system.