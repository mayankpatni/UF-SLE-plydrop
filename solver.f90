module solver
contains

!----------------------------------------------------------------
!   FUNCTION solve_disp_sparse FOR calculating displacement vector
!----------------------------------------------------------------
 SUBROUTINE solve_disp_sparse 
    use var_inputdata
    use var_analysis
    
      IMPLICIT NONE 
      include 'mkl_pardiso.fi'
!-----Internal variables------------------------------------------------
      REAL*8,allocatable            ::RESULTS(:)                  !SOLUTION
      INTEGER                       ::I
      
       !SOLVER PARAMTERS
      TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE:: PT(:)   ! INTERNAL SOLVER MEMORY POINTER
      
      INTEGER                 :: MAXFCT               ! MAXIMUM NUMBER OF FACTORS IN MEMORY (1)
      INTEGER                 :: MNUM                 ! NUMBER OF MATRIX ( 1<=MAXFCT) (1) 
      INTEGER                 :: MTYPE                ! MATRIX TYPE(11- REAL UNSYMMETRIC, 2 -REAL SYMM)
      INTEGER                 :: PHASE                ! CONTROLS THE EXECUTION OF SOLVER 
      INTEGER                 :: NRHS                 ! NUMBER OF RHS
      INTEGER                 :: ERROR                ! ERROR FLAG
      INTEGER                 :: MSGLVL               ! INFORMATION FROM SOLVER (0/1)
      INTEGER                 :: IDUM(1)              !
      INTEGER,ALLOCATABLE     :: IPARM(:)             ! ARRAY WITH VARIOUS PARAMETERS FOR MKL PARADISO (SIZE:64)
      REAL*8                  :: DDUM(1)
      
      !SETTING UP PARADISO CONTROL PARAMETERS
      
      ALLOCATE(IPARM(64))
      IPARM =     0
      IPARM(1) =  1 ! no solver default
      IPARM(2) =  2 ! fill-in reordering from METIS
      IPARM(4) =  0 ! no iterative-direct algorithm
      IPARM(5) =  0 ! no user fill-in reducing permutation
      IPARM(6) =  0 ! =0 solution on the first n components of x
      IPARM(8) =  2 ! numbers of iterative refinement steps
      IPARM(10) = 8 ! perturb the pivot elements with 1E-13 or 1E-8 (13 or 8)
      IPARM(11) = 1 ! use nonsymmetric permutation and scaling MPS
      IPARM(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
      IPARM(14) = 0 ! Output: number of perturbed pivots
      IPARM(18) =-1 ! Output: number of nonzeros in the factor LU
      IPARM(19) =-1 ! Output: Mflops for LU factorization
      IPARM(20) = 0 ! Output: Numbers of CG Iterations 
      IPARM(60) = 0 ! OOC 
      
      ERROR       = 0 
      MSGLVL      = 0 !PRINT STATISTICAL INFORMATION
      MAXFCT      = 1
      NRHS        = 1
      MNUM        = 1
      MTYPE       = -2
      
      !IF (M_TYPE == 'SP') MTYPE= 2     ! symmetric positive definite
      !IF (M_TYPE == 'SI') MTYPE=-2     ! symm indefinite
      !IF (M_TYPE == 'NS') MTYPE=11     ! non symmetric
      
      ALLOCATE(PT(64))
      DO I =1,64
          PT(I)%DUMMY = 0
      ENDDO
      
      
      allocate(RESULTS(tdof))
      RESULTS = 0.0D0
      del_uvect=0.0
      
! important link for incore and outofcore memory settings
!scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-431916D5-B76D-48A1-ABB5-1A0613FDC0FA.htm#GUID-431916D5-B76D-48A1-ABB5-1A0613FDC0FA  
    
      !SOLVER - PHASE 1 - REORDERING AND SYMBOLIC FACTORIZATION - MEMEORY ALLOCATION IS DONE
      PHASE = 11
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,DDUM,DDUM,ERROR)
      IF (error.ne.0) THEN
        WRITE(*,*) 'The following ERROR was detected (PHASE 1): ', error
        goto 1000
      ENDIF
      
      !SOLVER - PHASE 2 - FACTORIZATION
      PHASE = 22
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,DDUM,DDUM,ERROR)
      IF (error.ne.0) THEN
        WRITE(*,*) 'The following ERROR was detected (PHASE 2): ', error, max(iparm(15),iparm(16)+iparm(63)), max(iparm(15), iparm(16)+iparm(17))
        !do i=1,tdof
        !    write(11,*)k_val(diag_idx(i)),i,diag_idx(i)
        !enddo
        goto 1000
      ENDIF
      
      !SOLVER - PHASE 3 - SOLUTION
      IPARM(8) = 2 ! MAXIMUM NUMBER IF INTERATIVE REFINEMENT STEPS
      PHASE = 33
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,pvect_res,RESULTS,ERROR)
      del_uvect=RESULTS
                
      IF (error.ne.0) THEN
        WRITE(*,*) 'The following ERROR was detected (PHASE 3): ', error
        goto 1000
      ENDIF
      
1000  CONTINUE    
      !SOLVER - PHASE 4 - TERMINATION AND RELEASE
      PHASE = -1
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,DDUM,DDUM,ERROR)
      
      deallocate(RESULTS)
      deallocate(k_val,row_ptr,col_idx)                            !!!!!comment for modified NR
      RETURN
 END SUBROUTINE
 
 
!----------------------------------------------------------------
!   FUNCTION solve_eigen_value FOR Eigen value and eigen vector
!----------------------------------------------------------------
 SUBROUTINE solve_eigen_value 
    use var_inputdata
    use var_analysis
    
      IMPLICIT NONE 
      include 'mkl_pardiso.fi'
      
      !-----Input PARAMTERS------------------------------------------------
      Character*1             :: uplo                       ! Must be 'U' or 'L' or 'F' for upper, lower triangular, full matrices A and B respectively
      INTEGER                 :: m0                         ! specifies the initial guess for subspace dimension 0 < m0 ≤ n. Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]
      real*8,ALLOCATABLE      :: x_eigvec(:,:)              ! if fpm (5) =1, array x_eigvec (n, m) contains a basis of guess subspace where n is the order of the input matrix.
      INTEGER,ALLOCATABLE     :: fpm(:)                     ! to pass various parameters to Extended Eigensolver routines
      REAL*8                  :: emin,emax                  ! The lower and upper bounds of the interval to be searched for eigenvalues
      
      !-----Output PARAMTERS------------------------------------------------
      REAL*8                  :: epsout                     ! The lower and upper bounds of the interval to be searched for eigenvalues
      INTEGER                 :: loop                       ! On output, contains the number of refinement loop executed. Ignored on input.
      real*8,allocatable      :: eigenval(:)                ! Array of length m0. On output, the first m entries of eigenval are eigenvalues found in the interval.
      integer                 :: m                          ! total number of eigenvalues found in the interval [emin, emax]: 0 ≤ m ≤ m0.
      integer                 :: info                       ! If info=0, the execution is successful
      real*8,allocatable      :: res(:)                     ! Array of length m0. On exit, the first m components contain the relative residual vector
      
      m0=10
      
      ALLOCATE(fpm(128),eigenval(m0),res(m0),x_eigvec(tdof,m0))
      
      call feastinit (fpm)                                  ! Initialize Extended Eigensolver input parameters with default values
      
      !call dfeast_scsrgv(uplo, tdof, a, ia, ja, b, ib, jb, fpm, epsout, loop, emin, emax, m0, eigenval, x_eigvec, m, res, info)
      
      
      RETURN
 END SUBROUTINE

end module

