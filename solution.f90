module solution
    contains
    
!---------------------------------------------------------------------------
!   FUNCTION lin_static_3D FOR nonlinear static analysis with 3D integral
!---------------------------------------------------------------------------
    subroutine lin_static_3D
    use var_inputdata
    use var_analysis
    use k_mat
    use load_vector
    use solver
    use post_process
    
    implicit none
    real            ::t1,t2
    integer         ::i,t1_sys,t2_sys,t0_rate
    
    !!!!!!!!!!!!!!!!!!!!!!!! k matrix assembly !!!!!!!!!!!!!!!!!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
                
    ALLOCATE(MTR_R(tdof))
    DO i=1,tdof
        ALLOCATE(MTR_R(i)%VAL(INT(tdof*2.5/100)))
        ALLOCATE(MTR_R(i)%COL(INT(tdof*2.5/100))) 
        MTR_R(i)%VAL=0.0D0
        MTR_R(i)%NNZ=1
        MTR_R(i)%COL=0
        MTR_R(i)%COL(1)=i
    ENDDO
        
    call k_mat3D
    call kmat_BC
    call kmat_sparse
    DEALLOCATE(MTR_R)    
    
    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,*)'K matrix assembly',real(t2_sys-t1_sys)/t0_rate,'s',t2-t1,'s'
    
    !!!!!!!!!!!!!!!!!!!!!!!! Solve (3D int) !!!!!!!!!!!!!!!!!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    
    pvect_res(:)=fvect(:)
    call solve_disp_sparse
    call uvect_BC(del_uvect)
    uvect(:)=uvect(:)+del_uvect(:)
    
    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,*)'Solution',real(t2_sys-t1_sys)/t0_rate,'s',t2-t1,'s'
    
    write(13,*)'uvect'
    do i=tdof-3*maxexp+1,tdof
        write(13,*)uvect(i)
    enddo
    
    return
    end subroutine
    
end module