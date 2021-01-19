module load_vector
contains

!---------------------------------------------------------
!      FUNCTION loadvector FOR getting the load vector
!---------------------------------------------------------
	subroutine loadvector  
    use var_inputdata
    use var_analysis
    use integ_3D
    use struc_3D
    use mod_math
    use mod_gauss
    
    implicit none
    integer::i9,j9,k9,div,i11,i10,j10,ele_beam,ele_cs,count,ndnm,loadtyp,ndt
    integer::count1,cs,i,j,ndnm1,beam_nd,count2,ij,ij1,ij2,cs_node,ele_bm,bne,agpt,bgpt
    real*8::xp,zp,yp,eta,alpha,beta,pi,jac,wt,xa,xb,za,zb,wth,wtk
    real*8::loadx,loady,loadz,aa
    real*8,allocatable,dimension(:,:)   ::ft_int_beam,ft_int
    real*8,allocatable,dimension(:,:)   ::sf_int_beam
    real*8,allocatable,dimension(:)     ::xcord,cs_area
    
    allocate (sf(nne),dsf(nne),f(mexp),fa(mexp),fb(mexp))

    j=0
    pi=4.0*ATAN(1.0)
        
    do i11=1,np
        xp=loadinp(1,i11)    
        yp=loadinp(2,i11)   
        zp=loadinp(3,i11)   
        loadx=loadinp(4,i11) 
        loady=loadinp(5,i11)  
        loadz=loadinp(6,i11)                     
        loadtyp=loadinp(7,i11)
        
        call phy_to_norBeam(yp,eta,ele_beam)
        call expfun_beam(eta)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Point Load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (loadtyp.eq.1) then
            call phy_to_norCS_2_load(xp,zp,alpha,beta,ele_cs)
            call expfun_cs(ele_cs,alpha,beta)
            
            if (ele_cs.eq.0) then
                write(11,*)xp,yp,zp
                write(11,*)'check for the load location'
            else
                do i9=1,nne 
                ndt=(ele_beam-1)*(nne-1)+i9
                do i10=1,mexp
                    ndnm=elemcon_replace(ele_cs,i10+1,ndt)
                    if (ndnm.eq.0) cycle
                    
                    count1=sum(maxexp_nodes(1:ndt-1))
                    count=(3*count1)+(3*ndnm)
                
                    fvect(count-2)=fvect(count-2)+loadx*f(i10)*sf(i9)
                    fvect(count-1)=fvect(count-1)+loady*f(i10)*sf(i9)
                    fvect(count)=fvect(count)+loadz*f(i10)*sf(i9)
                enddo
                enddo
            endif
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! cross-section Surface Load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (loadtyp.eq.2) then
            allocate(ft_int(mexp,ane),cs_area(ane))
            ft_int=0.0d0
            cs_area=0.0d0
            
            do ele_cs=1,ane
            if (cs_ele_type(ele_cs,2).eq.3) then
                agpt=agp_tri
                bgpt=1
            elseif (cs_ele_type(ele_cs,2).eq.4) then
                agpt=agp_quad
                bgpt=agp_quad
            endif
            !do j10=1,ane_load_num(i11)                 ! for load applied at particular cs elements
            !ele_cs=ane_load(i11,j10)
                if (elem_beam_coll(ele_cs,2).gt.0.and.elem_beam_coll(ele_cs,2).le.ele_beam) cycle
                bne=(ele_beam-1)*ane+ele_cs             
                do j9=1,agpt
                    do k9=1,bgpt
                        if (cs_ele_type(ele_cs,2).eq.3) then
                            alpha=agpt_cs_t(j9)
                            beta=bgpt_cs_t(j9)
                            wth=wpt_cs_t(j9)
                            wtk=1.0d0  
                        elseif (cs_ele_type(ele_cs,2).eq.4) then
                            alpha=gpt_cs(j9)
                            beta=gpt_cs(k9)
                            wth=wpt_cs(j9)
                            wtk=wpt_cs(k9)
                        endif
                        call brickshapefun(bne,alpha,eta,beta)
                        call expfun_cs(ele_cs,alpha,beta)
                        jac=jac3D(1,1)*jac3D(3,3)-jac3D(3,1)*jac3D(1,3)
                        do i9=1,mexp
                            ft_int(i9,ele_cs)=ft_int(i9,ele_cs)+(jac*wth*wtk*f(i9))
							!if (i9.eq.9.and.ele_cs.eq.1) then
							!	write(11,*)ft_int(i9,ele_cs),alpha,beta
							!endif
							
                        enddo 
                        cs_area(ele_cs)=cs_area(ele_cs)+abs(jac)*wth*wtk
                    enddo
                enddo
            enddo
                       
            do ele_cs=1,ane
                if (elem_beam_coll(ele_cs,2).gt.0.and.elem_beam_coll(ele_cs,2).le.ele_beam) cycle
                do i9=1,nne
                    ndt=(ele_beam-1)*(nne-1)+i9
                    do i10=1,mexp
                        ndnm=elemcon_replace(ele_cs,i10+1,ndt)
                        if (ndnm.eq.0) cycle
                        count1=sum(maxexp_nodes(1:ndt-1))
                        count=(3*count1)+(3*ndnm)
                
                        fvect(count-2)=fvect(count-2)+loadx*ft_int(i10,ele_cs)*sf(i9)
                        fvect(count-1)=fvect(count-1)+loady*ft_int(i10,ele_cs)*sf(i9)
                        fvect(count)=fvect(count)+loadz*ft_int(i10,ele_cs)*sf(i9)
                    enddo
                enddo
             enddo
            
            deallocate(ft_int)
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! beam Surface Load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! will not work for tri elements in cs !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! need to change ft_int_beam calculation part !!!!!!!!!!!!!!
        elseif (loadtyp.eq.3) then
            
            allocate(ft_int_beam(mexp,ane),sf_int_beam(nne,ne))
            count2=0
            ij=0
            ft_int_beam=0.0d0
            sf_int_beam=0.0d0
            
            do cs_node=1,anode
                if(cscord(cs_node,3).eq.zp) then
                    count2=count2+1
                endif
            enddo
            allocate(xcord(count2))
            
            do cs_node=1,anode
                if(cscord(cs_node,3).eq.zp) then
                    ij=ij+1
                    xcord(ij)=cscord(cs_node,2)
                endif
            enddo
            
            call Sort(xcord,count2)
            
            do ij1=1,count2-1
                xp=0.5d0*(xcord(ij1)+xcord(ij1+1))
                jac=(xcord(ij1+1)-xcord(ij1))/2.0d0
                call phy_to_norCS_2_load(xp,zp,alpha,beta,ele_cs)
                do ij2=1,agp_quad
    	            alpha=gpt_cs(ij2)
     	            wt=wpt_cs(ij2)
                    call expfun_cs(ele_cs,alpha,beta)
                    do i9=1,mexp
                        ft_int_beam(i9,ele_cs)=ft_int_beam(i9,ele_cs)+(wt*f(i9)*jac)
                    enddo 
                enddo
            enddo
            
            do ij1=1,ne
                jac=(y((nne-1)*ij1+1)-y((nne-1)*(ij1-1)+1))/2.0d0
                do ij2=1,ngp
    	            eta=gpt_bm(ij2)
     	            wt=wpt_bm(ij2)
                    call expfun_beam(eta)
                    do i9=1,nne
                        sf_int_beam(i9,ij1)=sf_int_beam(i9,ij1)+(wt*sf(i9)*jac)
                    enddo 
                enddo
            enddo
            
            do ele_bm=1,ne
                do ele_cs=1,ane
                    do i9=1,nne
                        ndt=(ele_bm-1)*(nne-1)+i9
                        loadz=loadinp(6,i11)    !*sin(pi*y(ndt)/(ylengthf-ylengthi))              !!!!!! for Sinusiodal load !!!!!!!!!!!!!!!!!!
                        
                        do i10=1,mexp
                            ndnm=elemcon(ele_cs,i10+1)
                            if (ndnm.eq.0) cycle
                        
                            count1=sum(maxexp_nodes(1:ndt-1))
                            count=(3*count1)+(3*ndnm)
                
                            fvect(count-2)=fvect(count-2)+loadx*ft_int_beam(i10,ele_cs)*sf_int_beam(i9,ele_bm)
                            fvect(count-1)=fvect(count-1)+loady*ft_int_beam(i10,ele_cs)*sf_int_beam(i9,ele_bm)
                            fvect(count)=fvect(count)+loadz*ft_int_beam(i10,ele_cs)*sf_int_beam(i9,ele_bm)
                        enddo
                    enddo
                enddo
            enddo
            
         deallocate(xcord,ft_int_beam,sf_int_beam)
        endif
    enddo
    
    deallocate(f,fa,fb,sf,dsf)

     return
    end subroutine
 
    
    
!-----------------------------------------------------------
!      FUNCTION kmat_BC FOR getting Kmat after applying BCs
!-----------------------------------------------------------
	subroutine kmat_BC  
    use var_inputdata
    use var_analysis
    
    implicit none
    integer::row_nz,i,j,beam_nd,row_nnz,ii,ij,j1,cs_nd
    
    !do ii=1,tdof
    !    row_nnz=MTR_R(ii)%NNZ
    !    do ij=1,row_nnz
    !        j=MTR_R(ii)%COL(ij) 
    !        MTR_R(ii)%col_real(j)=ij
    !    enddo
    !enddo
    
    do i=1,num_bcinp
        beam_nd=node_bc(i,1)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Fixed BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (node_bc(i,2).eq.0.and.node_bc(i,3).eq.0) then           
           if (beam_nd.eq.1) then
               do j=1,3*sum(maxexp_nodes(1:beam_nd))
                    row_nz=MTR_R(j)%NNZ
                    MTR_R(j)%VAL(2:row_nz)=0.0d0
                    MTR_R(j)%VAL(1)=1.0d0
              enddo 
           elseif (beam_nd.eq.tnodes) then  
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                  row_nz=MTR_R(j)%NNZ
                  MTR_R(j)%VAL(2:row_nz)=0.0d0
                  MTR_R(j)%VAL(1)=1.0d0
                  
                  !do ii=1,j-1
                  !    ij=MTR_R(ii)%col_real(j)
                  !    if (ij.ne.0) then
                  !       MTR_R(ii)%VAL(ij)=0.0d0
                  !    endif
                  !enddo
                  
                  do ii=1,j-1
                    row_nnz=MTR_R(ii)%NNZ
                    do ij=2,row_nnz
                        if (MTR_R(ii)%COL(ij).eq.j) then
                            MTR_R(ii)%VAL(ij)=0.0d0
                            exit
                        endif
                    enddo
                  enddo
               enddo
           else
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                    row_nz=MTR_R(j)%NNZ
                    MTR_R(j)%VAL(2:row_nz)=0.0d0
                    MTR_R(j)%VAL(1)=1.0d0
                    
                    !do ii=1,j-1
                    !  ij=MTR_R(ii)%col_real(j)
                    !  if (ij.ne.0) then
                    !     MTR_R(ii)%VAL(ij)=0.0d0
                    !  endif
                    !enddo
                  
                    do ii=1,j-1
                        row_nnz=MTR_R(ii)%NNZ
                        do ij=2,row_nnz
                            if (MTR_R(ii)%COL(ij).eq.j) then
                                MTR_R(ii)%VAL(ij)=0.0d0
                                exit
                            endif
                        enddo
                    enddo 
               enddo
           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simply-supported BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (node_bc(i,2).eq.100.and.node_bc(i,3).eq.100) then
            if (beam_nd.eq.1) then
               do j=1,3*sum(maxexp_nodes(1:beam_nd))
                   if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                                   !!!!!! x and z constrained !!!!!!
                       row_nz=MTR_R(j)%NNZ
                       MTR_R(j)%VAL(2:row_nz)=0.0d0
                       MTR_R(j)%VAL(1)=1.0d0
                       
                       !do ii=1,j-1
                       !    ij=MTR_R(ii)%col_real(j)
                       !    if (ij.ne.0) then
                       !        MTR_R(ii)%VAL(ij)=0.0d0
                       !    endif
                       !enddo
                       
                       do ii=1,j-1
                        row_nnz=MTR_R(ii)%NNZ
                        do ij=2,row_nnz
                            if (MTR_R(ii)%COL(ij).eq.j) then
                                MTR_R(ii)%VAL(ij)=0.0d0
                                exit
                            endif
                        enddo
                      enddo
                   endif
               enddo 
           elseif (beam_nd.eq.tnodes) then  
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                    if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                                  !!!!!! x and z constrained !!!!!!
                      row_nz=MTR_R(j)%NNZ
                      MTR_R(j)%VAL(2:row_nz)=0.0d0
                      MTR_R(j)%VAL(1)=1.0d0
                      
                      !do ii=1,j-1
                      !     ij=MTR_R(ii)%col_real(j)
                      !     if (ij.ne.0) then
                      !         MTR_R(ii)%VAL(ij)=0.0d0
                      !     endif
                      !enddo
                  
                      do ii=1,j-1
                        row_nnz=MTR_R(ii)%NNZ
                        do ij=2,row_nnz
                            if (MTR_R(ii)%COL(ij).eq.j) then
                                MTR_R(ii)%VAL(ij)=0.0d0
                                exit
                            endif
                        enddo
                      enddo
                    endif                    
                enddo
           else
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                   if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                                         !!!!!! x and z constrained !!!!!! 
                        row_nz=MTR_R(j)%NNZ
                        MTR_R(j)%VAL(2:row_nz)=0.0d0
                        MTR_R(j)%VAL(1)=1.0d0
                        
                        !do ii=1,j-1
                        !   ij=MTR_R(ii)%col_real(j)
                        !   if (ij.ne.0) then
                        !       MTR_R(ii)%VAL(ij)=0.0d0
                        !   endif
                        !enddo
                  
                        do ii=1,j-1
                            row_nnz=MTR_R(ii)%NNZ
                            do ij=2,row_nnz
                                if (MTR_R(ii)%COL(ij).eq.j) then
                                    MTR_R(ii)%VAL(ij)=0.0d0
                                    exit
                                endif
                            enddo
                        enddo 
                   endif
                enddo
           endif
           

           ! if (beam_nd.eq.1) then
           !    do j=1,3*sum(maxexp_nodes(1:beam_nd))
           !        if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                           !!!!!! x and z constrained !!!!!!
           !            row_nz=MTR_R(j)%NNZ
           !            MTR_R(j)%VAL(2:row_nz)=0.0d0
           !            MTR_R(j)%VAL(1)=1.0d0
           !            
           !            do ii=1,j-1
           !             row_nnz=MTR_R(ii)%NNZ
           !             do ij=2,row_nnz
           !                 if (MTR_R(ii)%COL(ij).eq.j) then
           !                     MTR_R(ii)%VAL(ij)=0.0d0
           !                 endif
           !             enddo
           !           enddo
           !        endif
           !    enddo 
           !elseif (beam_nd.eq.tnodes) then  
           !    do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
           !         if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                          !!!!!! x and z constrained !!!!!!
           !           row_nz=MTR_R(j)%NNZ
           !           MTR_R(j)%VAL(2:row_nz)=0.0d0
           !           MTR_R(j)%VAL(1)=1.0d0
           !       
           !           do ii=1,j-1
           !             row_nnz=MTR_R(ii)%NNZ
           !             do ij=2,row_nnz
           !                 if (MTR_R(ii)%COL(ij).eq.j) then
           !                     MTR_R(ii)%VAL(ij)=0.0d0
           !                 endif
           !             enddo
           !           enddo
           !         endif                    
           !     enddo
           !else
           !    do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
           !        if (mod(j+1,3).eq.0) then                                           !!!!!! y constrained !!!!!! 
           !             row_nz=MTR_R(j)%NNZ
           !             MTR_R(j)%VAL(2:row_nz)=0.0d0
           !             MTR_R(j)%VAL(1)=1.0d0
           !       
           !             do ii=1,j-1
           !                 row_nnz=MTR_R(ii)%NNZ
           !                 do ij=2,row_nnz
           !                     if (MTR_R(ii)%COL(ij).eq.j) then
           !                         MTR_R(ii)%VAL(ij)=0.0d0
           !                     endif
           !                 enddo
           !             enddo 
           !        endif
           !    enddo
           !endif
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Constrained nodes (or edges) BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else                                            
            if (beam_nd.eq.1) then
                do j1=1,nexp+1
                    cs_nd=cs_node_bc(i,j1)
                    if (cs_nd.eq.0) then
                        goto 1010
                    endif
                    
                    do j=3*(cs_nd-1)+1,3*cs_nd
                        if (i.le.2) then
                            if (mod(j+1,3).eq.0) then                                   !!!!!!y constrained !!!!!!
                                row_nz=MTR_R(j)%NNZ
                                MTR_R(j)%VAL(2:row_nz)=0.0d0
                                MTR_R(j)%VAL(1)=1.0d0
                        
                                do ii=1,j-1
                                    row_nnz=MTR_R(ii)%NNZ
                                    do ij=2,row_nnz
                                        if (MTR_R(ii)%COL(ij).eq.j) then
                                            MTR_R(ii)%VAL(ij)=0.0d0
                                        endif
                                    enddo
                                enddo
                            endif                                                     !!!!!!x constrained !!!!!! uncomment to remove
                        elseif (i.eq.3) then
                            if (mod(j,3).eq.0) then                                   !!!!!!z constrained !!!!!!
                                row_nz=MTR_R(j)%NNZ
                                MTR_R(j)%VAL(2:row_nz)=0.0d0
                                MTR_R(j)%VAL(1)=1.0d0
                        
                                do ii=1,j-1
                                    row_nnz=MTR_R(ii)%NNZ
                                    do ij=2,row_nnz
                                        if (MTR_R(ii)%COL(ij).eq.j) then
                                            MTR_R(ii)%VAL(ij)=0.0d0
                                        endif
                                    enddo
                                enddo
                            endif
                        else
                            if (mod(j-1,3).eq.0) then                                   !!!!!!x constrained !!!!!!
                                row_nz=MTR_R(j)%NNZ
                                MTR_R(j)%VAL(2:row_nz)=0.0d0
                                MTR_R(j)%VAL(1)=1.0d0
                        
                                do ii=1,j-1
                                    row_nnz=MTR_R(ii)%NNZ
                                    do ij=2,row_nnz
                                        if (MTR_R(ii)%COL(ij).eq.j) then
                                            MTR_R(ii)%VAL(ij)=0.0d0
                                        endif
                                    enddo
                                enddo
                            endif
                        endif
                    enddo
1010             enddo
            elseif (beam_nd.eq.tnodes) then
                do j1=1,nexp+1
                    cs_nd=cs_node_bc(i,j1)
                    if (cs_nd.eq.0) then
                        goto 1012
                    endif
                    
                    do j=3*sum(maxexp_nodes(1:beam_nd-1))+3*(cs_nd-1)+1,3*sum(maxexp_nodes(1:beam_nd-1))+3*cs_nd
                        if (i.le.2) then
                            if (mod(j+1,3).eq.0) then                                   !!!!!!y constrained !!!!!!
                                row_nz=MTR_R(j)%NNZ
                                MTR_R(j)%VAL(2:row_nz)=0.0d0
                                MTR_R(j)%VAL(1)=1.0d0
                        
                                do ii=1,j-1
                                    row_nnz=MTR_R(ii)%NNZ
                                    do ij=2,row_nnz
                                        if (MTR_R(ii)%COL(ij).eq.j) then
                                            MTR_R(ii)%VAL(ij)=0.0d0
                                        endif
                                    enddo
                                enddo
                            endif                                                     !!!!!!x constrained !!!!!! uncomment to remove
                        elseif (i.eq.3) then
                            if (mod(j,3).eq.0) then                                   !!!!!!z constrained !!!!!!
                                row_nz=MTR_R(j)%NNZ
                                MTR_R(j)%VAL(2:row_nz)=0.0d0
                                MTR_R(j)%VAL(1)=1.0d0
                        
                                do ii=1,j-1
                                    row_nnz=MTR_R(ii)%NNZ
                                    do ij=2,row_nnz
                                        if (MTR_R(ii)%COL(ij).eq.j) then
                                            MTR_R(ii)%VAL(ij)=0.0d0
                                        endif
                                    enddo
                                enddo
                            endif
                        else
                            if (mod(j-1,3).eq.0) then                                   !!!!!!x constrained !!!!!!
                                row_nz=MTR_R(j)%NNZ
                                MTR_R(j)%VAL(2:row_nz)=0.0d0
                                MTR_R(j)%VAL(1)=1.0d0
                        
                                do ii=1,j-1
                                    row_nnz=MTR_R(ii)%NNZ
                                    do ij=2,row_nnz
                                        if (MTR_R(ii)%COL(ij).eq.j) then
                                            MTR_R(ii)%VAL(ij)=0.0d0
                                        endif
                                    enddo
                                enddo
                            endif
                        endif
                        !if (mod(j-1,3).eq.0) then                                   !!!!!!x constrained !!!!!!
                        !    row_nz=MTR_R(j)%NNZ
                        !    MTR_R(j)%VAL(2:row_nz)=0.0d0
                        !    MTR_R(j)%VAL(1)=1.0d0
                        !
                        !    do ii=1,j-1
                        !        row_nnz=MTR_R(ii)%NNZ
                        !        do ij=2,row_nnz
                        !            if (MTR_R(ii)%COL(ij).eq.j) then
                        !                MTR_R(ii)%VAL(ij)=0.0d0
                        !            endif
                        !        enddo
                        !    enddo
                        !endif                                                     !!!!!!x constrained !!!!!! uncomment to remove
                    enddo
1012            enddo
            else
                do j1=1,nexp+1
                    cs_nd=cs_node_bc(i,j1)
                    if (cs_nd.eq.0) then
                        goto 1011
                    endif
                    
                    do j=3*sum(maxexp_nodes(1:beam_nd-1))+3*(cs_nd-1)+1,3*sum(maxexp_nodes(1:beam_nd-1))+3*cs_nd
                        if (mod(j-1,3).eq.0) then                                   !!!!!!x constrained !!!!!!
                            row_nz=MTR_R(j)%NNZ
                            MTR_R(j)%VAL(2:row_nz)=0.0d0
                            MTR_R(j)%VAL(1)=1.0d0
                        
                            do ii=1,j-1
                                row_nnz=MTR_R(ii)%NNZ
                                do ij=2,row_nnz
                                    if (MTR_R(ii)%COL(ij).eq.j) then
                                        MTR_R(ii)%VAL(ij)=0.0d0
                                    endif
                                enddo
                            enddo
                        endif                                                     !!!!!!x constrained !!!!!! uncomment to remove
                    enddo
1011            enddo
            endif
            
        endif
    enddo
    return
    end subroutine
    
!------------------------------------------------------------------
!      FUNCTION uvect_BC FOR enforcing uvect to zero to model BCs
!------------------------------------------------------------------
	subroutine uvect_BC(vector)  
    use var_inputdata
    use var_analysis
    
    implicit none
    integer::row_nz,i,j,beam_nd,row_nnz,ii,ij,j1,cs_nd
    real*8,allocatable,dimension(:)::vector
    
    do i=1,num_bcinp
        beam_nd=node_bc(i,1)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Fixed BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (node_bc(i,2).eq.0.and.node_bc(i,3).eq.0) then           
           if (beam_nd.eq.1) then
               do j=1,3*sum(maxexp_nodes(1:beam_nd))
                    vector(j)=0.0d0
              enddo 
           else
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                    vector(j)=0.0d0
               enddo
           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simply-supported BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (node_bc(i,2).eq.100.and.node_bc(i,3).eq.100) then
            if (beam_nd.eq.1) then
               do j=1,3*sum(maxexp_nodes(1:beam_nd))
                   if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                           !!!!!! x and z constrained !!!!!!
                       vector(j)=0.0d0
                   endif
               enddo 
           elseif (beam_nd.eq.tnodes) then  
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                    if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                          !!!!!! x and z constrained !!!!!!
                        vector(j)=0.0d0
                    endif                    
                enddo
           else
               do j=3*sum(maxexp_nodes(1:beam_nd-1))+1,3*sum(maxexp_nodes(1:beam_nd))
                   if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                           !!!!!! x and z constrained !!!!!! 
                        vector(j)=0.0d0
                   endif
               enddo
           endif
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Constrained nodes (or edges) BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else                                            
            if (beam_nd.eq.1) then
                do j1=1,nexp+1
                    cs_nd=cs_node_bc(i,j1)
                    if (cs_nd.eq.0) then
                        goto 1010
                    endif
                    
                    do j=3*(cs_nd-1)+1,3*cs_nd
                        if (i.le.2) then
                            if (mod(j+1,3).eq.0) then                                   !!!!!!y constrained !!!!!!
                               vector(j)=0.0d0
                            endif                                                     !!!!!!x constrained !!!!!! uncomment to remove
                        elseif (i.eq.3) then
                            if (mod(j,3).eq.0) then                                   !!!!!!z constrained !!!!!!
                                vector(j)=0.0d0
                            endif
                        else
                            if (mod(j-1,3).eq.0) then                          !!!!!! x constrained !!!!!!
                            vector(j)=0.0d0
                            endif
                        endif   
                    enddo
1010             enddo
            elseif (beam_nd.eq.tnodes) then
                do j1=1,nexp+1
                    cs_nd=cs_node_bc(i,j1)
                    if (cs_nd.eq.0) then
                        goto 1012
                    endif
                    
                    do j=3*sum(maxexp_nodes(1:beam_nd-1))+3*(cs_nd-1)+1,3*sum(maxexp_nodes(1:beam_nd-1))+3*cs_nd
                        if (i.le.2) then
                            if (mod(j+1,3).eq.0) then                                   !!!!!!y constrained !!!!!!
                               vector(j)=0.0d0
                            endif                                                     !!!!!!x constrained !!!!!! uncomment to remove
                        elseif (i.eq.3) then
                            if (mod(j,3).eq.0) then                                   !!!!!!z constrained !!!!!!
                                vector(j)=0.0d0
                            endif
                        else
                            if (mod(j-1,3).eq.0) then                          !!!!!! x constrained !!!!!!
                            vector(j)=0.0d0
                            endif
                        endif
                    enddo
1012            enddo
            else
                do j1=1,nexp+1
                    cs_nd=cs_node_bc(i,j1)
                    if (cs_nd.eq.0) then
                        goto 1011
                    endif
                    
                    do j=3*sum(maxexp_nodes(1:beam_nd-1))+3*(cs_nd-1)+1,3*sum(maxexp_nodes(1:beam_nd-1))+3*cs_nd
                        if (mod(j-1,3).eq.0) then                          !!!!!! x constrained !!!!!!
                        vector(j)=0.0d0
                        endif   
                    enddo
1011            enddo
            endif
            
        endif
    enddo
              
     return
    end subroutine
    
    
end module