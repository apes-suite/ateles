Ateles
======

Ateles is a discontinuous Galerkin solver operating on
octree meshes with cubical elements and a state representation
with Legendre polynomials as a basis.

It allows for arbitrary high order discretizations provides the
possibility to represent embedded boundaries.
Various equation systems are implemented, besides acoustic and
flow systems also the Maxwell equations for example.

Documentation
-------------

See the [documentation](https://apes-suite.github.io/ateles/)
for more details.

Project organization
--------------------

The actual sources of the project are found in the ateles-source
repository, which is included here in the `atl` subdirectory.
This is a wrapper repository binding together all the parts
required for compilation.

Use `git clone --recurse-submodules` when cloning this repository to fetch the
gathered subdirectories from the various repositories.

Compiling
---------

Prerequisite for building the solver is an installed Python, a Fortran compiler
and a MPI library.
For compilation you need to point `FC` to the appropiate MPI compiler wrapper.
(Usually `export FC=mpif90`).

The solver can then be built with

```
bin/waf configure build
```

To install it, run:

```
bin/waf install
```

Run `bin/waf --help` to see all options.

Development Environments
------------------------

We support two methods to conveniently install a development environment
with the APES tools:

* via the Spack manager with the packages in [apes-spack](https://github.com/apes-suite/apes-spack)
* via a Python virtual environment as provided in [apes-pyenv](https://github.com/apes-suite/apes-pyenv)

Please see the respective READMEs for instructions on how to use one of these
methods.

Testing
-------

Unit tests are included in the build process via the waf-unit-tests.
They are always run during the compilation, unless you deactivate them
with the `--notests` option.

System tests are done with the help of [PySys-Test](https://pysys-test.github.io/pysys-test/),
which is included in the `apes-pyenv`.
The test cases to be run are found in `atl/examples`.
Running all system tests can be achieved with an activated `apes-pyenv` and either
the `ateles` executable in your `$PATH` (for example by installing it to the
`$VIRTUAL_ENV` prefix) or indicated via the `$APES_ATELES` variable, by doing:

```
    cd atl/examples
    pysys.py run
```

The test command can also be run in the subdirectories of `atl/examples` to
only run the tests within  those directories.
