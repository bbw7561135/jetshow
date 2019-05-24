## jetshow
Full Stokes (_IQUV_) synchrotron radiation transfer in inhomogeneous AGN jets. The fields can be specified both
analytically as C++ classes and as the results of the simulations (i.e. on some grid).


## Installation and building

Possibly you need to install ``git``, ``cmake``, ``boost``, ``Eigen`` and ``CGAL`` (for interpolation of simulations
results).

``git clone https://ipashchenko@github.com/ipashchenko/jetshow.git``

``cd jetshow``

``mkdir build-dir; cd build_dir``

``cmake ..``

``make``

## Using

Python class ``JetModelZoom`` from ``jetmodel.py`` handles analytical model with non-uniform pixel.