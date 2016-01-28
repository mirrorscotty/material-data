material-data
=============
Code to calculate material properties of foods as a function of temperature.

This includes a set of Choi-Okos equations (Choi and Okos 1986) as well as 
several functions related to diffusivity and drying rate of pasta. The
diffusivity and isotherm models come primarily from *The Handbook of Food
Engineering* (Second Edition), Chapter 10, and functions for pasta composition
and thermal properties are mostly for the transport model from Zhu 2012.

Contents
--------
* `composition` - Choi-Okos functions
* `diffusivity` - Effective diffusivity functions for extruded durum semolina and also for vapor diffusivity
* `fem-cheat` - Stuff to calculate pressure curves without using any finite element analysis
* `glass-transition` - Gordon-Taylor and Kwei models for glass transition temperature
* `isotherms` - Isotherm data for durum semolina and some other assorted products
* `math` - Numerical Laplace and inverse Laplace transform routines and a function to calculate L2-norm
* `mechanical` - Viscoelasticity and capillary pressure
* `pasta` - Assorted functions from Zhu 2012
* `test` - Test programs designed to output data to CSV files
* `unifac` - Beginnings of an implementation of the UNIFAC model

Building
--------
This code was designed to be compiled as a library and statically linked. To build just the library, execute `make material-data.a`. Documentation can be built using Doxygen using `make doc`.

Usage
-----
Include `material-data.h` in any C file that requires functions provided by this library, and instruct the compiler to link the executable against `material-data.a`.

Dependencies
------------
Compilation requires GCC and GNU Make. Building documentation requires Doxygen and LaTeX to generate PDF output. This code also requires the matrix library found [here](https://github.com/mirrorscotty/matrix).

This project is hosted at https://github.com/mirrorscotty/material-data, and any updates will be posted there.
