# icFlow4

tested with Qt 5.15.2 on Ubuntu 21.10

icFlow3 is an experimental 2D FEM code, conceived primarily for modeling fracture due to indentation. The main features are

* the use of variational integration and Incremental Potential Contact, as described in [https://ipc-sim.github.io/](https://ipc-sim.github.io/)

* fracture algorithm is inspired by [ARCSim](http://graphics.berkeley.edu/resources/ARCSim/)

* [PPR cohesive zones](https://www.sciencedirect.com/science/article/pii/S0013794412000690) are used for ductile fracture

* sparse matrix assembly algorithm is explained in this [CodeProject article](https://www.codeproject.com/Articles/5314545/Construction-of-Sparse-Matrices-for-Finite-Element)

## Required libraries:

libgomp1

libtbb-dev

libhdf5-dev

libgmsh-dev

libeigen3-dev

libboost-dev

libvtk9-dev

libvtk9-qt-dev

Additionally, [MOSEK](https://mosek.com) solver must be installed and its path should be updated manually in CMakeLists.txt.

When building, ensure that Qt libraries have the same version and location: Qt5Charts_DIR, Qt5Core_DIR, Qt5Gui_DIR, Qt5Widgets_DIR, and Qt5_DIR.

![screenshot](/s.png?raw=true)

