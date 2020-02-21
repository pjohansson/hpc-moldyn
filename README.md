# Overview

This is a molecular dynamics integrator. To use:

```
# Print help
./md
```

The program uses GROMOS-87 formatted configuration files and only reads the positions. It is currently a bit particular about the formatting of the box size vector. See the example configuration file in `include` for specifics. An example parameter and forcefield file is also included.

# Compilation instructions

Compilation requires `cmake` and an MPI installation. Example compilation using OpenMPI:

```
cmake [SOURCE_DIR] -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++
make
```
