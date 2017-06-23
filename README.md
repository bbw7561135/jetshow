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
        
    * ``vfield`` - velocity field with possible types:
        * ``const_central`` with parameter ``gamma`` (lorentz factor of bulk
        motion). This corresponds to radial velocity field.
        
    * ``nfield`` - density field with possible types:
        * ``bk`` - power law dependence on distance from SBH with parameters
        ``n_1`` and ``n_n``.
        
        ``n1`` - density [cm^(-3)] on 1 pc from SBH, `n_n`` - minus exponent of
        distance dependence.
  
 * ``integration`` describes parameters of numerical solving of ODEs.
 
    * ``tau_max`` - maximum optical depth to reach. All part of jet further than
    region with this optical depth will be discarded in transport.
    * ``dl_max_pc`` - maximum length of adaptive step [pc].
    * ``log10_tau_min`` - log10 of minimum optical depth to calculate flux.
    * ``n`` - preliminary number of steps. The real number will be different.
    
 * ``output`` describes output files.
 
    * ``file_tau`` - name of the file with image of optical depth.
    * ``file_length`` - name of the file with image of geometrical depth.
    * ``file_i`` - name of the file with image of Stokes ``I`` intensity (in
    units of [Jy/pixel]).
    