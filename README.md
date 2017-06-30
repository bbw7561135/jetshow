## jetshow
Synchrotron radiation transfer in inhomogeneous AGN jet model (a.k.a. Blandford & Königl jet model).


## Installation and building

Possibly you need to install ``git``, ``cmake``, ``boost`` and ``Eigen`` (C++
libraries).

``git clone https://ipashchenko@bitbucket.org/ipashchenko/jetshow.git``

``cd jetshow``

``mkdir build-dir; cd build_dir``

``cmake ..``

``make``

## Using

Configuration is managed via ``config.json`` file which is explained below.
After building and optional editing of this file we ran executable:

``./jetshow``

Result files are appeared in the same directory.



## Parameters of ``config.json``
 
 * ``observation`` describes observer-related parameters such as
 frequency, angle of jet axis to observer LOS, redshift
 of the source.
 
 * ``image`` describes resulting image-specific parameters.
 
 * ``jet`` describes parameters of the jet.
 
    * ``geometry`` with possible ``type``:
        * ``cone`` with parameters ``angle`` and ``scale_pc``.
        * (**not using for now**)``cylinder`` with parameters ``r`` and scale.
        
        The last one determines geometrical path in case of infinite
        intersections.
        
    * ``bfield`` - magnetic field with possible types:
        * ``spiral_conical`` with parameters ``b_1`` and ``pitch_angle``
        * ``radial_conical`` with parameters ``b_1`` and ``n_b``.
        
        ``b1`` - value of magnetic field [G] on 1 pc from SBH, ``pitch_angle`` -
        ratio of toroidal to longitudinal components (in observer's frame),
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
        * ``const_central`` with parameter ``gamma`` (lorentz factor of bulk
        motion). This corresponds to radial velocity field.
        
    * ``nfield`` - density field with possible types:
        * ``bk`` - power law dependence on distance from SBH with parameters
        ``n_1`` and ``n_n``.
        
        ``n1`` - density [cm^(-3)] on 1 pc from SBH, `n_n`` - minus exponent of
        distance dependence.
  
 * ``integration`` describes parameters of numerical solving of ODEs.
 
    * ``step_type`` - Step type with possible values:
        * ``constant`` with parameters ``tau_max``, ``n``, ``tau_n_min``,
        ``n_tau_max``
        * ``adaptive`` with parameters ``dl_max_pc``, ``n``, ``abs`` and ``rel``
        errors.
    * ``tau_max`` - maximum optical depth to reach. All part of jet further than
    region with this optical depth will be discarded in transport.
    * ``dl_max_pc`` - maximum length of adaptive step [pc].
    * ``log10_tau_min`` - log10 of minimum optical depth to calculate flux.
    * ``n`` - preliminary number of steps. The real number will be different for
    adaptive steps.
    * ``tau_n_min`` - border value of the optical depth: if depth is less then
     this value - ``n`` steps will be used with ``constant`` step type when
     integrating Stokes parameters. If depth is higher - value higher then ``n``
     will be used.
    * ``n_tau_max`` - number of steps to take while integrating Stokes
    parameters when optical depth is higher or equal to ``tau_max``.
    * ``..._abs/rel_error...`` - absolute or relative error for adaptive type.
    
 * ``output`` describes output files.
 
    * ``file_tau`` - name of the file with image of optical depth.
    * ``file_length`` - name of the file with image of geometrical depth [pc].
    * ``file_i`` - name of the file with image of Stokes *I* intensity
    [Jy/pixel].
    
    
## Warnings

* Currently, models of random B-field could be used only in single threaded run.
To do it - comment out line 27 at ``Observation.cpp`` and re-compile.
 
* Using of ``adaptive`` integration type slows down the running time seriously
in case of ``RandomPointBField``.