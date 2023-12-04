Ateles
======

Ateles is a discontinuous Galerkin solver operating on
octree meshes with cubical elements and a state representation
with Legendre polynomials as a basis.

It allows for arbitrary high order discretizations provides the
possibility to represent embedded boundaries.
Various equation systems are implemented, besides acoustic and
flow systems also the Maxwell equations for example.

The actual sources of the project are found in the ateles-source
repository, which is included here in the `atl` subdirectory.
This is a wrapper repository binding together all the parts
required for compilation.

Use `git clone --recurse-submodules` when cloning this repository to fetch the
gathered subdirectories from the various repositories.

Prerequisite for building the solver is an installed Python, Fortran compiler
and MPI library. For compilation you need to point `FC` to the appropiate MPI
compiler wrapper. (Usually `export FC=mpif90`).

The solver can then be built with

```
bin/waf configure build
```

To install it, run:

```
bin/waf install
```

Run `bin/waf --help` to see all options.

Documentation
-------------

See the [documentation](https://geb.inf.tu-dresden.de/doxy/ateles/index.html)
for more details.
