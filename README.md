## jetshow
Full Stokes (_IQUV_) synchrotron radiation transfer in inhomogeneous AGN jets. The fields can be specified both
analytically as C++ classes and as the results of the simulations (i.e. on some grid).


## Installation and building

Possibly you need to install ``git``, ``cmake``, ``boost``, ``Eigen`` and ``CGAL`` (for interpolation of simulations
results).

``git clone https://ipashchenko@bitbucket.org/ipashchenko/jetshow.git``

``cd jetshow``

``mkdir build-dir; cd build_dir``

``cmake ..``

``make``

## Using

### Warning: Currently choosing a model should be done in ``main.cpp`` file. Configuration via ``json``-config file is highly likely broken now.

Configuration is managed via ``config.json`` file which is explained below.
After building and optional editing of this file we ran executable:

``./jetshow``

Result files should appear in the same directory. Results could be displayed
using ``plotting.py`` script. E.g.:

![Stokes I contours, color of frac. linear polarization and direction with
value of linear polarization.](https://www.dropbox.com/s/adr38w9f6ay2j8b/bk_jet.png)




## Parameters of ``config.json`` (outdated)
 
 ``observation`` describes observer-related parameters such as
 frequency, angle of jet axis to observer LOS, redshift
 of the source.
 
 ---
 
 ``image`` describes resulting image-specific parameters.
 
 ---
 
 ``jet`` describes parameters of the jet.
 
 * ``geometry`` with possible ``type``:
     * ``cone`` with parameters ``angle`` (half-opening angle) and ``scale_pc``.
     * (**not using for now**)``cylinder`` with parameters ``r`` and scale.
        
       The last one determines geometrical path in case of infinite
       intersections.
        
 * ``bfield`` - magnetic field with possible types:
     * ``spiral_conical`` with parameters ``b_1`` and ``pitch_angle``
     * ``radial_conical`` with parameters ``b_1`` and ``n_b``.
     * ``toroidal`` with parameters ``b_1`` and ``n_b``
        
       ``b1`` - value of magnetic field [G] on 1 pc from SBH, ``pitch_angle``
        \- ratio of toroidal to longitudinal components (in observer's frame),
        ``n_b`` - minus exponent of distance dependence.
     * ``random_fraction`` - fraction of random B-field. If ``0.0`` then no
         random field.
     * ``model_of_randomness`` - how random B-field is implemented. Possible
       values are:
         * ``points`` - in each point random B-field has random direction.
         * ``cells`` - B-field consists of cells with random but constant
           value and direction inside single cell. On cross-section of jet at
           distance 1 pc there are nearly ``n_cells_1pc`` cells.
        
 * ``vfield`` - velocity field with possible types:
     * ``const_flat`` with parameter ``gamma`` (lorentz factor of bulk
       motion). This corresponds to flat velocity field.
     * ``const_central`` with parameter ``gamma`` (lorentz factor of bulk
       motion). This corresponds to radial velocity field.
        
 * ``nfield`` - density field with possible types:
     * ``bk`` - power law dependence on distance from SBH with parameters
       ``n_1`` and ``n_n``.
        
       ``n1`` - density [cm^(-3)] on 1 pc from SBH, `n_n`` - minus exponent of
       distance dependence.
  
 ---
  
 ``integration`` describes parameters of numerical solving of ODEs.

 * ``step_type`` - Step type with possible values:
     * ``constant`` with parameters ``tau_max``, ``n``, ``tau_n_min``,
       ``n_tau_max``
     * ``adaptive`` with parameters ``dl_max_pc``, ``n``, ``abs`` and
       ``rel`` errors.
 * ``tau_max`` - maximum optical depth to reach. All part of jet further
       than region with this optical depth will be discarded in transport.
 * ``dl_max_pc`` - maximum length of adaptive step [pc].
 * ``log10_tau_min`` - log10 of minimum optical depth to calculate flux.
 * ``n`` - preliminary number of steps. The real number will be different
   for adaptive steps.
 * ``tau_n_min`` - border value of the optical depth: if depth is less then
   this value - ``n`` steps will be used with ``constant`` step type when
   integrating Stokes parameters. If depth is higher - value higher then
   ``n`` will be used.
 * ``n_tau_max`` - number of steps to take while integrating Stokes
   parameters when optical depth is higher or equal to ``tau_max``.
 * ``..._abs/rel_error...`` - absolute or relative error for adaptive type.
    
 ---
    
 ``output`` describes output files.
 
 * ``file_tau`` - name of the file with image of optical depth.
 * ``file_length`` - name of the file with image of geometrical depth [pc].
 * ``file_i`` - name of the file with image of Stokes *I* intensity
   [Jy/pixel]. The same for each Stokes parameter (*Q*, *U* and *V*).
 
 ---
    
 ``calculate`` chooses what to calculate. Possible values are:
 
 * ``I`` - calculate only ``I`` Stokes parameter
 * ``full`` - calculate all Stokes parameters and optical depth
 * ``tau`` - calculate only stripe of the optical depth along the jet, which
   is used while iteratively finding the best image parameters.
    
    
## Warnings

 
* Using of ``adaptive`` integration type slows down the running time in case of
``PointsRandomVectorBField``.