File descriptions:
NOTE: All elements are geometrically nonlinear

    Amplitude.*     Stores time versus amplitude arrays. Used to create time dependent SPC or LOAD

    CBAR.*          Bar/truss element
    
    CBAR1.*         OUTDATED Equivalent to CBAR, but different solution method. And does not allow for variable restlength
    
    CBAR7.*         OUTDATED
    
    CBEAM.*         Beam element. Based on https://link.springer.com/article/10.1007/s00466-003-0422-7
    
    CHEXA.*         Hexahedral element with growth
    
    CHEXA7.*        Hexahedral element with growth and additional dof for cell density
    
    CONTACT#.*      Different node-surface contact implementations. CONTACT4.* performs best
    
    CQUAD.*         Quadrilateral element with growth
    
    CQUAD5.*        Quadrilateral axisymmetric element around Z=0 (coordinates in r,Z plane) with growth.
    
    CQUAD7.*        Quadrilateral element with growth and additional dof for cell density
    
    DataContainer.* Contains all data that is updated during simulations, such as dof, state variables, time, dt
    
    Element.*       Parent object for all elements
    
    ElementHelper.* Contains functions shared by many elements, such as quadrature calculations.
    
    LOAD.*          Point load
    
    MAT_GENT.*      Gent material model
    
    MAT_LE.*        Saint-Venant Kirchhoff material model
    
    MAT_NEOH.*      Neo hookean material model
    
    MAT_OGDEN.*     Ogden material model
    
    MAT_VISC.*      Viscous material model (nonlinear damper)
    
    Material.*      Parent object for all material models
    
    ModelContainer.*    Parent object for model description
    
    MPC.*           Multiple point constraint
    
    Node.*          Node object
    
    OutputHelper.*  Contains functions for writing output, such as paraview files
    
    PBAR.*          Bar property
    
    PBAR1.*         Bar property with internal force
    
    PBEAM.*         Beam property
    
    PDAMAGE.*       Damage property (damage softening, very basic, no length scale)
    
    PMAXWELL.*      Maxwell property
    
    PRESSURE.*      Pressure load (follows deformation)
    
    Property.*      Parent object for all properties
    
    PSOLID10.*      Isotropic prescribed growth property
    
    PSOLID11.*      Anisotropic prescribed fiber growth property
    
    PSOLID12.*      Anisotropic prescribed area growth property
    
    PSOLID20.*      Isotropic stretch driven growth property
    
    PSOLID21.*      Anisotropic stretch driven fiber growth property
    
    PSOLID71.*      Orthotropic growth and cell migration/diffusion property (used in https://www.sciencedirect.com/science/article/pii/S0022509617310347)
    
    SOL2.*          Newton-Raphson solver
    
    Solver.*        Parent object for all solvers
    
    SPC.*           Single point constraint
