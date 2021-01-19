!     Last edited:  MP  17 Oct 2019 
!===================================================================================================================================
!      CUF Implementation using Serendipity Lagrange Expansion 
!      Laminated composites (LW), linear, 3D integration, 3D structures, ply drops
!-----------------------------------------------------------------------------------------------------------------------------------
!      UF-SLE-lin-tapered PROGRAMME
!===================================================================================================================================
    PROGRAM UF_SLE_lin_tapered
    use var_inputdata
    use var_analysis
    use read_input_data
    use struc_3D
	use integ_3D
    use k_mat
    use mod_gauss
    use mod_math
    use load_vector
    use solution
    use solver
    use post_process
    
    implicit none
    integer::i,ti_sys,tf_sys,t_rate,t1_sys,t2_sys,t0_rate,el,j,el1,c,c1
    integer, allocatable::elemcon_replace1(:,:),tri_conn(:,:)
    real:: t1,t2,tf,ti
    real*8::alp,bet
	
    open(11,file="output/SLE-output_comp.txt")
    open(12,file="output/para_SLE_comp.vtk")
    open(13,file="output/uvect_comp.txt")
    open(14,file="output/post_p_out.dat")
    
     write(11,60)
     write(*,70)'Step','Elapsed Time'
60   format(10X,'**********************************************'/,10X,&
                '*               CUF FEM                      *'/,10X,&
                '**********************************************'/)
70  format(A17,15X,A13)   
80  format(A22,11X,F15.2,A1,F15.2,A1)
    
    CALL SYSTEM_CLOCK(COUNT_RATE=t_rate); CALL SYSTEM_CLOCK(COUNT=ti_sys); call cpu_time(ti)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Read input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
     
    call read_input                                                         ! sle_elem_replace subroutine change - include tri element
    call PRINTI(elemcon,ane,mexp+1)
 
    call brick3D

    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,80)'Read input',real(t2_sys-t1_sys)/t_rate,'sec',t2-t1,'sec'
    
   call PRINTI(elemcon,ane,mexp+1)
   call PRINTI(cs_node_bc,num_bcinp,nexp+1)
   
    write(11,*) ' Node no    Y-cord'
   do i=1,tnodes
        write(11,10)i,y(i)
   enddo
10 format(2x,i4,5x,f21.10)

   !allocate(elemcon_replace1(ane,mexp+1))
   !do i=1,tnodes
   !    write(11,*)'beam node=',i
   !    elemcon_replace1(:,:)=elemcon_replace(:,:,i)
   !    call PRINTI(elemcon_replace1,ane,mexp+1)
   !enddo

    call phy_to_norCS_1    
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	allocate (f(mexp),fa(mexp),fb(mexp))
!	allocate (tri_conn(phynormesh*phynormesh,8))
!	tri_conn=0
!	tri_conn(:,2)=2
!	tri_conn(:,3)=2
!	tri_conn(:,4)=1
!	tri_conn(:,5)=2
!	el=1
!	write(14,*)'el=1'
!	write(13,*)'el=1'
!	do i=1,tpts(2)
!		alp=normcs_tri(i,1)
!		bet=normcs_tri(i,2)
!		call expfun_cs(el,alp,bet)
!		write(14,*)phycs(i,2),phycs(i,3)   !alp,bet		
!		write(13,*)f(1),f(2),f(3)
!	enddo
!	
!	write(14,*)'el=2'
!	write(13,*)'el=2'
!	el=2
!	do i=tpts(2)+1,tpts(2)*2
!		alp=normcs_tri(i-tpts(2),1)
!		bet=normcs_tri(i-tpts(2),2)
!		call expfun_cs(el,alp,bet)
!		write(14,*)phycs(i,2),phycs(i,3)   !alp,bet		
!		write(13,*)f(1),f(2),f(3)
!	enddo
!	
!	write(14,*)'el=3'
!	write(13,*)'el=3'
!	el=3
!	do i=tpts(2)*2+1,tpts(2)*3
!		alp=normcs_tri(i-tpts(2)*2,1)
!		bet=normcs_tri(i-tpts(2)*2,2)
!		call expfun_cs(el,alp,bet)
!		write(14,*)phycs(i,2),phycs(i,3)   !alp,bet		
!		write(13,*)f(1),f(2),f(3)
!	enddo
!	
!	el=0 
!    el1=0
!	c=0
!	c1=0
!    do i=1,phynormesh
!		c1=c1+phynormesh+2-i
!        do j=1,phynormesh-el
!            el1=el1+1
!			tri_conn(el1,1)=el1		!el1
!			tri_conn(el1,6)=j+c			!el1+i-1
!			tri_conn(el1,7)=j+c+1			!el1+i
!			tri_conn(el1,8)=j+c1			!el1+6
!        enddo
!        el=el+1
!		c=tri_conn(el1,7)
!	enddo
!	
!	el=0
!	c=0
!	c1=0
!	do i=1,phynormesh-1
!		c1=c1+phynormesh+2-i
!        do j=1,phynormesh-1-el
!            el1=el1+1
!			tri_conn(el1,1)=el1
!			tri_conn(el1,6)=j+c+1
!			tri_conn(el1,7)=j+c1+1
!			tri_conn(el1,8)=j+c1
!        enddo
!        el=el+1
!		c=c+phynormesh+2-i
!	enddo
!	
!	DO I = 1, phynormesh*phynormesh
!        WRITE(11,1001)  tri_conn(I,:)
!	1001 format(24(I4,1x))
!	END DO
!	
!	!call PRINTI(tri_conn,phynormesh*phynormesh,8)
!	
!	deallocate(f,fa,fb)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (fullorpost.eq.2) then
        allocate (uvect(tdof))
        uvect=0.0d0
        open(8,file="input/uvect.dat")
        do i=1,tdof
            read(8,*)uvect(i)
        enddo
        close(8)
    else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
       CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
       allocate (fvect(tdof))
       fvect=0.0d0

       call loadvector
   
       write(13,*)'fvect'
       do i=3,tdof,3
           write(13,*)fvect(i)
       enddo
   
       CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
       write(*,80)'Load vector',real(t2_sys-t1_sys)/t0_rate,'sec',t2-t1,'sec'
       
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Solution !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Linear Static (analysis_type=1) !!!!!!!!!!!!!!!!!!!!!
        allocate (uvect(tdof),del_uvect(tdof),pvect_res(tdof))
        allocate (trcmat(6,6),sf(nne),dsf(nne),f(mexp),fa(mexp),fb(mexp))
		allocate(kmat(band_width,tdof))
        uvect=0.0d0
        del_uvect=0.0d0
        pvect_res=0.0d0

        call lin_static_3D 
        
        deallocate (trcmat,sf,dsf,f,fa,fb,del_uvect,pvect_res,fvect)
    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Post-process !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate (trcmat(6,6),sf(nne),dsf(nne),f(mexp),fa(mexp),fb(mexp))
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!! paraview !!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    !
    !call post_processing
    !call PARA_WRITE_8
    !
    !CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    !write(*,80)'Post-process paraview',real(t2_sys-t1_sys)/t0_rate,'sec',t2-t1,'sec'
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!! at desired points !!!!!!!!!!!!!!!!!!!!!!!
   CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    
   call post_processing2
   !call post_processing2_section
   !call post_processing2_abn
    
   CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
   write(*,80)'Post-process points',real(t2_sys-t1_sys)/t0_rate,'sec',t2-t1,'sec'

   
   deallocate (trcmat,sf,dsf,f,fa,fb,uvect)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !call cond_number  
   
    close(11)
    close(12)
    close(13)
    close(14)
    
    CALL SYSTEM_CLOCK(COUNT=tf_sys); call cpu_time(tf)
    write(*,80)'Full analysis',real(tf_sys-ti_sys)/t_rate,'sec',tf-ti,'sec'

	stop

    END PROGRAM UF_SLE_lin_tapered