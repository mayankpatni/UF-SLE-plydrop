# UF-SLE-plydrop
# Softwares required to run the code
# 1. Visual studio 2017 Enterprise 
# 2. Intel parallel studio XE 2018 or 2019 version
# Visual studio 2017 is compatible with Intel parallel studio 2018 version, and won’t work with the previous intel versions. 

# K_mat.f90 - element stiffness matrix computation and global assembly
# load_vector.f90 - Load application and Boundary conditions
# mod_gauss.f90 - Integration points 
# post_process.f90 - post-processing data (displacement, stress, strains)
# solver.f90 - Pardiso and eigenvalue solver
# struc_3D.f90 - Variable kinematics 3D element based on Unified formulation - definition, geometry and element connectivity
# integ_2D_cs.f90 - cross-section shape functions (Serendipity Lagrange expansion)
# integ_1D_beam.f90 - 1D Lagrange shape functions (along the beam axis)
# UF_SLE_lin_tapered.f90 - main file (calling all the subroutines)
# read_input_SLE.f90 - input data processing (geometry - beam and cross-sections, element connectivity, materials)
# integ_3D.f90 - 3D integration ((for non-prismatic and curved elements)
# MOD_math.f90 - mathematical functions definitions
# variable_declare_SLE.f90 - input/output variables declaration

# The following functions/subroutines are developed/extended by the researcher during his time spent within the Wind Blade Research Hub (WBRH) (a collaboration between the University of Bristol and ORE Catapult).

1. subroutine brickshapefun (in module struc_3D.f90)
2. subroutine brick3D (in module struc_3D.f90)
3. subroutine cs_elem_3D (in module read_input_SLE.f90)
4. subroutine sle_elem_connect (in module read_input_SLE.f90)
5. subroutine gauss_cs_tri (in module mod_gauss.f90)
6. subroutine expfun_cs (in module integ_3D.f90) (for tri element)
