File descriptions:\
NOTE: All elements are geometrically nonlinear

    CBARX.*         Extended bar element with additional parameters, such as state, mechanism, timeToNextEvent, etc.
    
    CBARX2.*        Same as CBARX.*, but assumes only X-displacement is active, and only computes the corresponding component of residual vector and stiffness matrix.
    
    CBEAMX.*        Extended beam element with additional parameters, such as state, mechanism, timeToNextEvent, etc.
    
    MECH_EL01.*     Element mechanism. Simple detachment and reattachment to same nodes
    
    MECH_EL02.*     Element mechanism. Simple detachment and reattachment to different nodes based on element restlength
    
    MECH_EL03.*     Element mechanism. Dynein mechanism
    
    MECH_EL04.*     Element mechanism. Myosin mechanism
    
    MECH_EL05.*     OUTDATED. Myosin mechanism
    
    MECH_EL06.*     Element mechanism. Bell model mechanism, similar to MECH_EL02
    
    MECH_MT01.*     Microtubule mechanism. Polymerization and depolymerization
    
    MECH_MT02.*     Microtubule mechanism. Extends MECH_MT01.*, and adds interaction with growth cone node.
    
    Mechanism.*     Parent object for mechanisms.
    
    MT.*            Microtubule
    
    NodeX.*         Extended node with additional parameter, such as elPlus and elMinus
    
    SOL2X.*         Extended Newton-Raphson solver that facilates mechanisms
    
    
