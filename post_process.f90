module post_process
contains

!-------------------------------------------------------------------------
!      FUNCTION post_processing2 FOR post processing at desired points
!-------------------------------------------------------------------------
	subroutine post_processing2
    use var_inputdata
    use var_analysis
    use mod_math
    
    implicit none
    integer::k3,k4,ele_cs,alpbet_cnt,k2,i,npost_p,nne_post,tnodes_post,tgauss_beam
    character*1::mark
    real*8::xm1,zm1,alpha,beta
    real*8,ALLOCATABLE, DIMENSION (:)::nodeval,alpha_p,beta_p
    real*8,allocatable,dimension(:,:,:)::disp_p,strn_p_G,strs_p_G,strs_p_Le,strs_p_Ls,strn_p_G_gauss,strs_p_G_gauss,strs_p_Le_gauss,strs_p_Ls_gauss,strn_p_G_1,Git_temp_sec_1
    real*8,allocatable,dimension(:,:)::post_p,u_G,eps_G,sig_G,sig_Le,sig_Ls,xyz_post,eps_G_gauss,sig_G_gauss,sig_Le_gauss,sig_Ls_gauss,xyz_post_gauss,eps_G_1
    integer, ALLOCATABLE, DIMENSION (:)::ele_cs_p
    allocate (alpha_p(4),beta_p(4),ele_cs_p(4))
    
    
    open(5,file="input/post_p.dat")
    read(5,*)mark
    read(5,*)npost_p,be_div_post_p
    allocate (post_p(npost_p,2))
    read(5,*)mark
    do i=1,npost_p
        read(5,*)post_p(i,1),post_p(i,2)
    enddo   
    
    nne_post=be_div_post_p*(nne-1)+1
    tnodes_post=ne*(nne_post-1)+1
    tgauss_beam=ne*ngp
    
    allocate (strs_p_Le(tnodes_post,9,npost_p),strs_p_Ls(tnodes_post,9,npost_p),strn_p_G(tnodes_post,9,npost_p),strs_p_G(tnodes_post,9,npost_p))
    allocate (disp_p(tnodes_post,6,npost_p),strn_p_G_1(tnodes_post,2,npost_p))
    allocate (strs_p_Le_gauss(tgauss_beam,9,npost_p),strs_p_Ls_gauss(tgauss_beam,9,npost_p),strn_p_G_gauss(tgauss_beam,9,npost_p),strs_p_G_gauss(tgauss_beam,9,npost_p))
	
    
    disp_p=0.0d0
    strs_p_Le=0.0d0
    strs_p_Ls=0.0d0
    strn_p_G=0.0d0
    strs_p_G=0.0d0
	strn_p_G_1=0.0d0
    strs_p_Le_gauss=0.0d0
    strs_p_Ls_gauss=0.0d0
    strn_p_G_gauss=0.0d0
    strs_p_G_gauss=0.0d0
    

        allocate (eps_G(tnodes_post,6),u_G(tnodes_post,3),xyz_post(tnodes_post,3),eps_G_1(tnodes_post,6))
        allocate (sig_G(tnodes_post,6),sig_Le(tnodes_post,6),sig_Ls(tnodes_post,6),nodeval(tnodes_post))
        
        allocate (eps_G_gauss(tgauss_beam,6),xyz_post_gauss(tgauss_beam,3))
        allocate (sig_G_gauss(tgauss_beam,6),sig_Le_gauss(tgauss_beam,6),sig_Ls_gauss(tgauss_beam,6))
		allocate (Git_temp_sec_1(tnodes_post,3,mexp))
		
        do k3=1,npost_p
            xm1=post_p(k3,1)
            zm1=post_p(k3,2)

            call phy_to_norCS_2_post(xm1,zm1,alpha_p,beta_p,ele_cs_p,alpbet_cnt)            
           
            do k4=1,alpbet_cnt
                alpha=alpha_p(k4)
                beta=beta_p(k4)
                ele_cs=ele_cs_p(k4)
				!write(11,*)alpha,beta,ele_cs
                
                call post_disp_strs_strn(xm1,zm1,alpha,beta,ele_cs,nne_post,xyz_post,u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval,eps_G_1,Git_temp_sec_1) 
                call post_strs_strn_gauss(xm1,zm1,alpha,beta,ele_cs,xyz_post_gauss,eps_G_gauss,sig_G_gauss,sig_Le_gauss,sig_Ls_gauss)   
                
                if(k3.gt.1) then
                    if(post_p(k3-1,2).eq.post_p(k3,2)) then
                        write(14,*)sig_G(((tnodes_post-1)*50/100)+1,2),alpha,beta
                        write(14,*)sig_G(((tnodes_post-1)*50/100)+1,4),ele_cs
                    endif 
                endif
            
                
                do k2=1,tnodes_post
                    disp_p(k2,1,k3)=xyz_post(k2,1)
                    disp_p(k2,2,k3)=xyz_post(k2,2) 
                    disp_p(k2,3,k3)=xyz_post(k2,3)
                    
                    do i=1,3
                        disp_p(k2,i+3,k3)=(disp_p(k2,i+3,k3)*(k4-1)+u_G(k2,i))/(k4*1.0d0)
                    enddo
                    
                    strn_p_G(k2,1,k3)=xyz_post(k2,1)
                    strn_p_G(k2,2,k3)=xyz_post(k2,2)
                    strn_p_G(k2,3,k3)=xyz_post(k2,3)
                    
                    strs_p_G(k2,1,k3)=xyz_post(k2,1)
                    strs_p_G(k2,2,k3)=xyz_post(k2,2)
                    strs_p_G(k2,3,k3)=xyz_post(k2,3)
                    
                    strs_p_Le(k2,1,k3)=xyz_post(k2,1)
                    strs_p_Le(k2,2,k3)=xyz_post(k2,2)
                    strs_p_Le(k2,3,k3)=xyz_post(k2,3)
                    
                    strs_p_Ls(k2,1,k3)=xyz_post(k2,1)
                    strs_p_Ls(k2,2,k3)=xyz_post(k2,2) 
                    strs_p_Ls(k2,3,k3)=xyz_post(k2,3)
                    
                    do i=1,6
                        strn_p_G(k2,i+3,k3)=(strn_p_G(k2,i+3,k3)*(k4-1)+eps_G(k2,i))/(k4*1.0d0)                      
                        strs_p_G(k2,i+3,k3)=(strs_p_G(k2,i+3,k3)*(k4-1)+sig_G(k2,i))/(k4*1.0d0)
                        strs_p_Le(k2,i+3,k3)=(strs_p_Le(k2,i+3,k3)*(k4-1)+sig_Le(k2,i))/(k4*1.0d0)
                        strs_p_Ls(k2,i+3,k3)=(strs_p_Ls(k2,i+3,k3)*(k4-1)+sig_Ls(k2,i))/(k4*1.0d0)
					enddo
					strn_p_G_1(k2,1,k3)=(strn_p_G(k2,1,k3)*(k4-1)+eps_G_1(k2,1))/(k4*1.0d0)    
					strn_p_G_1(k2,2,k3)=(strn_p_G(k2,2,k3)*(k4-1)+eps_G_1(k2,2))/(k4*1.0d0)
                enddo 
                
                do k2=1,tgauss_beam
                    strn_p_G_gauss(k2,1,k3)=xyz_post_gauss(k2,1)
                    strn_p_G_gauss(k2,2,k3)=xyz_post_gauss(k2,2)
                    strn_p_G_gauss(k2,3,k3)=xyz_post_gauss(k2,3)
                    
                    strs_p_G_gauss(k2,1,k3)=xyz_post_gauss(k2,1)
                    strs_p_G_gauss(k2,2,k3)=xyz_post_gauss(k2,2)
                    strs_p_G_gauss(k2,3,k3)=xyz_post_gauss(k2,3)
                    
                    strs_p_Le_gauss(k2,1,k3)=xyz_post_gauss(k2,1)
                    strs_p_Le_gauss(k2,2,k3)=xyz_post_gauss(k2,2)
                    strs_p_Le_gauss(k2,3,k3)=xyz_post_gauss(k2,3)
                    
                    strs_p_Ls_gauss(k2,1,k3)=xyz_post_gauss(k2,1)
                    strs_p_Ls_gauss(k2,2,k3)=xyz_post_gauss(k2,2) 
                    strs_p_Ls_gauss(k2,3,k3)=xyz_post_gauss(k2,3)
                    
                    do i=1,6
                        strn_p_G_gauss(k2,i+3,k3)=(strn_p_G_gauss(k2,i+3,k3)*(k4-1)+eps_G_gauss(k2,i))/(k4*1.0d0)                      
                        strs_p_G_gauss(k2,i+3,k3)=(strs_p_G_gauss(k2,i+3,k3)*(k4-1)+sig_G_gauss(k2,i))/(k4*1.0d0)
                        strs_p_Le_gauss(k2,i+3,k3)=(strs_p_Le_gauss(k2,i+3,k3)*(k4-1)+sig_Le_gauss(k2,i))/(k4*1.0d0)
                        strs_p_Ls_gauss(k2,i+3,k3)=(strs_p_Ls_gauss(k2,i+3,k3)*(k4-1)+sig_Ls_gauss(k2,i))/(k4*1.0d0)
                    enddo
                enddo 
            enddo                       
        enddo
        
        deallocate (u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval,xyz_post,eps_G_gauss,sig_G_gauss,sig_Le_gauss,sig_Ls_gauss,xyz_post_gauss)
    
    !write(14,*)'uy'
    !write(11,*)abs(disp_p(tnodes_post,4,1)),abs(disp_p(tnodes_post,5,1)),abs(disp_p(tnodes_post,6,1)),app_load   
    
    !write(14,*)'uy at load=',app_load
    !do k3=1,tnodes_post
    !    write(14,*)disp_p(k3,5,1)
    !enddo
    
        write(14,*)'ux, uy and uz along y at bottom face'
        do k3=1,tnodes_post
            write(14,*)disp_p(k3,4,6),disp_p(k3,5,6),disp_p(k3,6,6)
        enddo
        
        !write(14,*)'y _ gauss, sig yy along y at bottom face _ gauss'
        !do k3=1,tgauss_beam
        !    write(14,*)strs_p_G_gauss(k3,2,1),strs_p_G_gauss(k3,5,81)
        !enddo
    
        
        write(14,*)'y, sig yy along y at bottom face'
        do k3=1,tnodes_post
            write(14,*)strs_p_G(k3,2,1),strs_p_G(k3,5,1)
		enddo
  
		write(14,*)'sig xx, sig yy, sig zz along y at bottom face'
        do k3=1,tnodes_post
            write(14,*)strs_p_G(k3,4,1),strs_p_G(k3,5,1),strs_p_G(k3,6,1)
		enddo
		
		write(14,*)'tau yz along y at mid face'
        do k3=1,tnodes_post
            write(14,*)strs_p_G(k3,7,6)
		enddo

		write(14,*)'tau yz along z at mid span'
        do k3=1,npost_p
            write(14,*)strs_p_G(((tnodes_post-1)*50/100)+1,5,k3),strs_p_G(((tnodes_post-1)*50/100)+1,7,k3),strs_p_G(((tnodes_post-1)*50/100)+1,6,k3)
		enddo
        
		write(14,*)'ux, uy and uz along z at mid span'
        do k3=1,npost_p
            write(14,*)disp_p(((tnodes_post-1)*50/100)+1,4,k3),disp_p(((tnodes_post-1)*50/100)+1,5,k3),disp_p(((tnodes_post-1)*50/100)+1,6,k3)
		enddo
		
		write(14,*)'strain yz components and Eyz along z at mid span'
        do k3=1,npost_p
            write(14,*)strn_p_G_1(((tnodes_post-1)*50/100)+1,1,k3),strn_p_G_1(((tnodes_post-1)*50/100)+1,2,k3),strn_p_G(((tnodes_post-1)*50/100)+1,7,k3)
		enddo
		
		write(14,*)'dUy/dz'
        do k3=2,npost_p-1
            write(14,*)(disp_p(((tnodes_post-1)*50/100)+1,5,k3+1)-disp_p(((tnodes_post-1)*50/100)+1,5,k3-1))/(disp_p(((tnodes_post-1)*50/100)+1,3,k3+1)-disp_p(((tnodes_post-1)*50/100)+1,3,k3-1))
		enddo
		
		write(14,*)'dUz/dy'
        do k3=2,npost_p-1
            write(14,*)(disp_p(((tnodes_post-1)*50/100)+2,6,k3)-disp_p(((tnodes_post-1)*50/100),6,k3))/(disp_p(((tnodes_post-1)*50/100)+2,2,k3)-disp_p(((tnodes_post-1)*50/100),2,k3))
		enddo
		
		
        !write(14,*)'sig yz along y at interface top'
        !do k3=1,tnodes_post
        !    write(14,*)strs_p_G(k3,7,1),strs_p_Le(k3,7,1)
        !enddo
        !
        !write(14,*)'sig zz along y at interface top'
        !do k3=1,tnodes_post
        !    write(14,*)strs_p_G(k3,6,1),strs_p_Le(k3,6,1)
        !enddo
        
        !write(14,*)'sig yy and yz Global at (50%)',app_load
        !do k3=1,npost_p
        !    write(14,*)strs_p_G(((tnodes_post-1)*50/100)+1,3,k3),strs_p_G(((tnodes_post-1)*50/100)+1,5,k3),strs_p_G(((tnodes_post-1)*50/100)+1,7,k3)
        !enddo
        !
        ! write(14,*)'sig yy and yz Local at (50%)',app_load
        !do k3=1,npost_p
        !    write(14,*)strs_p_Le(((tnodes_post-1)*50/100)+1,5,k3),strs_p_Le(((tnodes_post-1)*50/100)+1,7,k3)
        !enddo

        
        !write(14,*)'sig yy and yz and zz at (90%)',app_load
        !do k3=1,npost_p
        !    write(14,*)strs_p_G(((tnodes_post-1)*90/100)+1,5,k3),strs_p_G(((tnodes_post-1)*90/100)+1,7,k3),SQRT(((strs_p_G(((tnodes_post-1)*90/100)+1,4,k3)-strs_p_G(((tnodes_post-1)*90/100)+1,5,k3))**2+(strs_p_G(((tnodes_post-1)*90/100)+1,5,k3)-strs_p_G(((tnodes_post-1)*90/100)+1,6,k3))**2+(strs_p_G(((tnodes_post-1)*90/100)+1,6,k3)-strs_p_G(((tnodes_post-1)*90/100)+1,4,k3))**2+6*(strs_p_G(((tnodes_post-1)*90/100)+1,7,k3)**2+strs_p_G(((tnodes_post-1)*90/100)+1,8,k3)**2+strs_p_G(((tnodes_post-1)*90/100)+1,9,k3)**2))/2)
        !enddo
    !
    ! write(14,*)'sigzz along y'
    !do k3=1,tnodes_post
    !    write(14,*)strs_p_G(k3,6,1)
    !enddo
    !
    ! write(14,*)'sigxx along y'
    !do k3=1,tnodes_post
    !    write(14,*)strs_p_G(k3,4,1)
    !enddo
    
    !write(14,*)'manual residual: dSyz_dz at z=0.999'
    !do k3=7,tnodes-6,6
    !    write(14,*)(strs_p_G(k3,7,npost_p-2)-strs_p_G(k3,7,npost_p-1))/(post_p(npost_p-2,2)-post_p(npost_p-1,2))
    !enddo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    DEALLOCATE (alpha_p,beta_p,ele_cs_p,disp_p,strs_p_G,strn_p_G,strs_p_Le,strs_p_Ls,post_p)
    close(5)
       
     return
    end subroutine
    
!----------------------------------------------------------------------------
!      FUNCTION post_disp_strs_strn FOR calculating disp, strains and stress
!----------------------------------------------------------------------------
	subroutine post_disp_strs_strn(xm1,zm1,alpha,beta,ele_cs,nne_post,xyz_post,u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval,eps_G_1,Git_temp_sec_1)        
    use var_inputdata
    use var_analysis
    use integ_3D
    use read_input_data
    use mod_math
    use struc_3D
    
    implicit none
    integer::i15,j15,c3,l15,o15,cn,ele_cs,ndnm,nd_var,i,nne_post,i14,j14,bne
    real*8::alpha,beta,eta,unode(3),uvw_xyz(9),epsnode(6),signode(6),signode_Le(6),signode_Ls(6),etainc,sljacdet,epsnode_mat(3,3),signode_mat(3,3),signode_Le_mat(3,3),signode_Ls_mat(3,3)
    real*8::eta_phy,xm1,zm1,tmin,tmax,epsnode_1(2)
    real*8::Git,Git_temp(3),Term_temp(3)
    real*8,ALLOCATABLE, DIMENSION (:)::nodeval
    real*8,allocatable,dimension(:,:)::eps_G,sig_G,sig_Le,sig_Ls,u_G,xyz_post,eps_G_1
	real*8,allocatable,dimension(:,:,:)::Git_temp_sec_1
        
    etainc=2.0d0/(nne_post-1)
    nodeval=0.0d0
    u_G=0.0d0
    eps_G=0.0d0
    sig_G=0.0d0
    sig_Le=0.0d0
    sig_Ls=0.0d0
	Git_temp_sec_1=0.0d0

    xyz_post=0.0d0

    tmin=y(1)                               !yloc (or t) at first node of beam division
    tmax=y(ne*(nne-1)+1)                    !yloc (or t) at last node of beam division
    call expfun_cs(ele_cs,alpha,beta)
    
    do l15=1,ne
        if (elem_beam_coll(ele_cs,2).gt.0.and.elem_beam_coll(ele_cs,2).le.l15) cycle
        bne=(l15-1)*ane+ele_cs
        do o15=1,nne_post
            cn=o15+(l15-1)*(nne_post-1)
            eta=-1.0d0+((o15-1)*etainc)   
            eta_phy=((y((nne-1)*l15+1)-y((nne-1)*(l15-1)+1))/2.0d0)*(eta+1.0d0)+y((nne-1)*(l15-1)+1)
            call parametric(xm1,eta_phy,tmin,tmax,zm1,xyz_post(cn,1),xyz_post(cn,2),xyz_post(cn,3))
            call expfun_beam(eta)
            call brickshapefun(bne,alpha,eta,beta)
			
			if (cn.eq.31.and.l15.eq.10) then
				if (xm1.eq.0.0d0.and.zm1.eq.0.0d0) then
					write(11,*)bne
					call PRINTF(jac3D_inv,3,3)
				endif
			endif
			
            
            call mat_stiff(ele_cs,eta_phy)
            unode=0.0d0
            epsnode=0.0d0
			epsnode_1=0.0d0
            uvw_xyz=0.0d0
            
            do i15=1,nne
                nd_var=(l15-1)*(nne-1)+i15
                do j15=1,mexp
                    ndnm=elemcon_replace(ele_cs,j15+1,nd_var)
                    if (ndnm.eq.0) cycle
                        
                    c3=(3*sum(maxexp_nodes(1:nd_var-1)))+(3*ndnm) 
                    
                    Git=f(j15)*sf(i15)
                    unode(1)=unode(1)+uvect(c3-2)*Git
                    unode(2)=unode(2)+uvect(c3-1)*Git
                    unode(3)=unode(3)+uvect(c3)*Git

                    
                    Term_temp(1)=fa(j15)*sf(i15) 
                    Term_temp(2)=f(j15)*dsf(i15)
                    Term_temp(3)=fb(j15)*sf(i15)
                    Git_temp=matmul(jac3D_inv,Term_temp)
					
					if (cn.eq.nd_var) then
						Git_temp_sec_1(cn,1,j15)=Git_temp(1)
						Git_temp_sec_1(cn,2,j15)=Git_temp(2)
						Git_temp_sec_1(cn,3,j15)=Git_temp(3)
					endif
					
                    do i14=1,3
                        do j14=1,3
                            uvw_xyz((i14-1)*3+j14)=uvw_xyz((i14-1)*3+j14)+Git_temp(i14)*uvect(c3-3+j14)
                        enddo
					enddo
					!if (cn.eq.31.and.l15.eq.10.and.i15.eq.4) then
					!	if (xm1.eq.0.0d0.and.zm1.eq.0.0d0) then
					!		write(14,*)Git_temp(3),Term_temp(1),Term_temp(3)  !uvw_xyz(8)
					!	endif
					!endif
				enddo
				
			enddo            
			
            nodeval(cn)=nodeval(cn)+1
            !!!!!!!!!!!!!!!!!! global displacements !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i=1,3
                u_G(cn,i)=(((nodeval(cn)-1)*u_G(cn,i))+unode(i))/nodeval(cn)
			enddo
	
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!! strain vector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            epsnode(1)=uvw_xyz(1)                                                                                       !epsxx               
            epsnode(2)=uvw_xyz(5)                                                                                       !epsyy
            epsnode(3)=uvw_xyz(9)                                                                                       !epszz
            epsnode(4)=uvw_xyz(6)+uvw_xyz(8)                                                                            !epsyz
            epsnode(5)=uvw_xyz(7)+uvw_xyz(3)                                                                            !epsxz
            epsnode(6)=uvw_xyz(4)+uvw_xyz(2)                                                                            !epsxy
			epsnode_1(1)=uvw_xyz(6)
			epsnode_1(2)=uvw_xyz(8)

            !!!!!!!!!!!!!!! compute material stiffness and local cys inverse !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call matinv3(loc_e,loc_e_inv)
            call matinv3(loc_struc,loc_struc_inv)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!! Global stress vector and tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode=matmul(trcmat,epsnode)
            signode_mat(1,1)=signode(1);    signode_mat(1,2)=signode(6);    signode_mat(1,3)=signode(5)
            signode_mat(2,1)=signode(6);    signode_mat(2,2)=signode(2);    signode_mat(2,3)=signode(4)
            signode_mat(3,1)=signode(5);    signode_mat(3,2)=signode(4);    signode_mat(3,3)=signode(3)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!! Local stress tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!! Local csys wrt to material csys!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode_Le_mat=matmul(transpose(loc_e_inv),matmul(signode_mat,loc_e_inv))
            signode_Le(1)=signode_Le_mat(1,1);  signode_Le(2)=signode_Le_mat(2,2);  signode_Le(3)=signode_Le_mat(3,3)
            signode_Le(4)=signode_Le_mat(2,3);  signode_Le(5)=signode_Le_mat(1,3);  signode_Le(6)=signode_Le_mat(1,2)
            !!!!!!!!!!!!!!!! Local csys wrt to structure csys !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode_Ls_mat=matmul(transpose(loc_struc_inv),matmul(signode_mat,loc_struc_inv))
            signode_Ls(1)=signode_Ls_mat(1,1);  signode_Ls(2)=signode_Ls_mat(2,2);    signode_Ls(3)=signode_Ls_mat(3,3)
            signode_Ls(4)=signode_Ls_mat(2,3);  signode_Ls(5)=signode_Ls_mat(1,3);    signode_Ls(6)=signode_Ls_mat(1,2)
          
            do i=1,6
                eps_G(cn,i)=(((nodeval(cn)-1)*eps_G(cn,i))+epsnode(i))/nodeval(cn)
                sig_G(cn,i)=(((nodeval(cn)-1)*sig_G(cn,i))+signode(i))/nodeval(cn)
                sig_Le(cn,i)=(((nodeval(cn)-1)*sig_Le(cn,i))+signode_Le(i))/nodeval(cn)
                sig_Ls(cn,i)=(((nodeval(cn)-1)*sig_Ls(cn,i))+signode_Ls(i))/nodeval(cn)
			enddo
			eps_G_1(cn,1)=(((nodeval(cn)-1)*eps_G_1(cn,1))+epsnode_1(1))/nodeval(cn)
			eps_G_1(cn,2)=(((nodeval(cn)-1)*eps_G_1(cn,2))+epsnode_1(2))/nodeval(cn)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		enddo
	enddo  
           
	
!!!!!!!!!!!!!*******************************************************************!!!!!!!!!!!!
     !deallocate(jacob) 
     
    return
    end subroutine
    
!----------------------------------------------------------------------------
!      FUNCTION post_strs_strn_gauss FOR calculating disp, strains and stress
!----------------------------------------------------------------------------
	subroutine post_strs_strn_gauss(xm1,zm1,alpha,beta,ele_cs,xyz_post_gauss,eps_G_gauss,sig_G_gauss,sig_Le_gauss,sig_Ls_gauss)        
    use var_inputdata
    use var_analysis
    use integ_3D
    use read_input_data
    use mod_math
    use struc_3D
    use mod_gauss
    
    implicit none
    integer::i15,j15,c3,l15,o15,cn,ele_cs,ndnm,nd_var,i,nne_post,i14,j14,bne
    real*8::alpha,beta,eta,uvw_xyz(9),epsnode(6),signode(6),signode_Le(6),signode_Ls(6),sljacdet,epsnode_mat(3,3),signode_mat(3,3),signode_Le_mat(3,3),signode_Ls_mat(3,3)
    real*8::eta_phy,xm1,zm1,tmin,tmax
    real*8::Git_temp(3),Term_temp(3)
    real*8,allocatable,dimension(:,:)::eps_G_gauss,sig_G_gauss,sig_Le_gauss,sig_Ls_gauss,xyz_post_gauss
    
    eps_G_gauss=0.0d0
    sig_G_gauss=0.0d0
    sig_Le_gauss=0.0d0
    sig_Ls_gauss=0.0d0

    xyz_post_gauss=0.0d0

    tmin=y(1)                               !yloc (or t) at first node of beam division
    tmax=y(ne*(nne-1)+1)                    !yloc (or t) at last node of beam division
    call expfun_cs(ele_cs,alpha,beta)
    
    do l15=1,ne
        if (elem_beam_coll(ele_cs,2).gt.0.and.elem_beam_coll(ele_cs,2).le.l15) cycle
        bne=(l15-1)*ane+ele_cs
        do o15=1,ngp
            eta=gpt_bm(o15)
            cn=o15+(l15-1)*ngp
            eta_phy=((y((nne-1)*l15+1)-y((nne-1)*(l15-1)+1))/2.0d0)*(eta+1.0d0)+y((nne-1)*(l15-1)+1)
            call parametric(xm1,eta_phy,tmin,tmax,zm1,xyz_post_gauss(cn,1),xyz_post_gauss(cn,2),xyz_post_gauss(cn,3))
            call expfun_beam(eta)
            call brickshapefun(bne,alpha,eta,beta)
            
            call mat_stiff(ele_cs,eta_phy)
         
            epsnode=0.0d0   
            uvw_xyz=0.0d0
            
            do i15=1,nne
                nd_var=(l15-1)*(nne-1)+i15
                do j15=1,mexp
                    ndnm=elemcon_replace(ele_cs,j15+1,nd_var)
                    if (ndnm.eq.0) cycle
                        
                    c3=(3*sum(maxexp_nodes(1:nd_var-1)))+(3*ndnm) 
                    
                    Term_temp(1)=fa(j15)*sf(i15) 
                    Term_temp(2)=f(j15)*dsf(i15)
                    Term_temp(3)=fb(j15)*sf(i15)
                    Git_temp=matmul(jac3D_inv,Term_temp)

                    do i14=1,3
                        do j14=1,3
                            uvw_xyz((i14-1)*3+j14)=uvw_xyz((i14-1)*3+j14)+Git_temp(i14)*uvect(c3-3+j14)
                        enddo
                    enddo
                enddo
            enddo        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!! strain vector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            epsnode(1)=uvw_xyz(1)                                                                                       !epsxx               
            epsnode(2)=uvw_xyz(5)                                                                                       !epsyy
            epsnode(3)=uvw_xyz(9)                                                                                       !epszz
            epsnode(4)=uvw_xyz(8)+uvw_xyz(6)                                                                            !epsyz
            epsnode(5)=uvw_xyz(7)+uvw_xyz(3)                                                                            !epsxz
            epsnode(6)=uvw_xyz(4)+uvw_xyz(2)                                                                            !epsxy

            !!!!!!!!!!!!!!! compute material stiffness and local cys inverse !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call matinv3(loc_e,loc_e_inv)
            call matinv3(loc_struc,loc_struc_inv)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!! Global stress vector and tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode=matmul(trcmat,epsnode)
            signode_mat(1,1)=signode(1);    signode_mat(1,2)=signode(6);    signode_mat(1,3)=signode(5)
            signode_mat(2,1)=signode(6);    signode_mat(2,2)=signode(2);    signode_mat(2,3)=signode(4)
            signode_mat(3,1)=signode(5);    signode_mat(3,2)=signode(4);    signode_mat(3,3)=signode(3)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!! Local stress tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!! Local csys wrt to material csys!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode_Le_mat=matmul(transpose(loc_e_inv),matmul(signode_mat,loc_e_inv))
            signode_Le(1)=signode_Le_mat(1,1);  signode_Le(2)=signode_Le_mat(2,2);  signode_Le(3)=signode_Le_mat(3,3)
            signode_Le(4)=signode_Le_mat(2,3);  signode_Le(5)=signode_Le_mat(1,3);  signode_Le(6)=signode_Le_mat(1,2)
            !!!!!!!!!!!!!!!! Local csys wrt to structure csys !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode_Ls_mat=matmul(transpose(loc_struc_inv),matmul(signode_mat,loc_struc_inv))
            signode_Ls(1)=signode_Ls_mat(1,1);  signode_Ls(2)=signode_Ls_mat(2,2);    signode_Ls(3)=signode_Ls_mat(3,3)
            signode_Ls(4)=signode_Ls_mat(2,3);  signode_Ls(5)=signode_Ls_mat(1,3);    signode_Ls(6)=signode_Ls_mat(1,2)
          
            do i=1,6
                eps_G_gauss(cn,i)=epsnode(i)
                sig_G_gauss(cn,i)=signode(i)
                sig_Le_gauss(cn,i)=signode_Le(i)
                sig_Ls_gauss(cn,i)=signode_Ls(i)
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
    enddo  
                    
!!!!!!!!!!!!!*******************************************************************!!!!!!!!!!!!
     
    return
    end subroutine
    
!---------------------------------------------------------------------------------------------------------------
!      FUNCTION post_processing2_abn FOR post processing at all gauss points in the structure (alpha,beta,eta)
!---------------------------------------------------------------------------------------------------------------
	subroutine post_processing2_abn
    use var_inputdata
    use var_analysis
    use integ_3D
    use read_input_data
    use mod_math
    use mod_gauss
    use struc_3D
    
    implicit none
    integer::total_gpt,i,ig,j,ja,jb,i14,j14,i15,j15,gptnum,bne,nd_var,ndnm,vec,c3,agpt,bgpt,ele_tri
    real*8::xm,zm,tmin,tmax,alpha,beta,eta,eta_phy,uvw_xyz(9),epsnode(6),signode(6),signode_Le(6),signode_Ls(6),Term_temp(3),Git_temp(3)
    real*8::signode_mat(3,3),signode_Le_mat(3,3),signode_Ls_mat(3,3)
    real*8,allocatable,dimension(:,:)::xyz_post,strn_G_gp,strs_G_gp,strs_Le_gp,strs_Ls_gp
    
    total_gpt=(count(cs_ele_type(:,2).eq.3)*agp_tri+count(cs_ele_type(:,2).eq.4)*agp_quad*agp_quad)*ngp*ne
    tmin=y(1)                               !yloc (or t) at first node of beam division
    tmax=y(ne*(nne-1)+1)                    !yloc (or t) at last node of beam division
    
    allocate (xyz_post(total_gpt,3),strs_Le_gp(total_gpt,6),strs_Ls_gp(total_gpt,6),strn_G_gp(total_gpt,6),strs_G_gp(total_gpt,6))
    
    xyz_post=0.0d0
    strn_G_gp=0.0d0
    strs_G_gp=0.0d0
    strs_Le_gp=0.0d0
    strs_Ls_gp=0.0d0
    xm=0.0d0
    zm=0.0d0
    gptnum=0

    do i=1,ne
        do ig=1,ngp
            eta=gpt_bm(ig)
            call expfun_beam(eta)
            eta_phy=((y((nne-1)*i+1)-y((nne-1)*(i-1)+1))/2.0d0)*(eta+1.0d0)+y((nne-1)*(i-1)+1)
            do j=1,ane
                if (cs_ele_type(j,2).eq.3) then
                    agpt=agp_tri
                    bgpt=1
                elseif (cs_ele_type(j,2).eq.4) then
                    agpt=agp_quad
                    bgpt=agp_quad
                endif
                if (elem_beam_coll(j,2).gt.0.and.elem_beam_coll(j,2).le.i) cycle
                bne=(i-1)*ane+j  
                do ja=1,agpt
                    do jb=1,bgpt
                        gptnum=gptnum+1
                        if (cs_ele_type(j,2).eq.3) then
                            alpha=agpt_cs_t(ja)
                            beta=bgpt_cs_t(ja)
                        elseif (cs_ele_type(j,2).eq.4) then
                            alpha=gpt_cs(ja)
                            beta=gpt_cs(jb)
                        endif
                        call expfun_cs(j,alpha,beta)
                        xm=f(1)*cscord(elemcon(j,2),2)+f(2)*cscord(elemcon(j,3),2)+f(3)*cscord(elemcon(j,4),2)+f(4)*cscord(elemcon(j,5),2)
                        zm=f(1)*cscord(elemcon(j,2),3)+f(2)*cscord(elemcon(j,3),3)+f(3)*cscord(elemcon(j,4),3)+f(4)*cscord(elemcon(j,5),3)
                        call parametric(xm,eta_phy,tmin,tmax,zm,xyz_post(gptnum,1),xyz_post(gptnum,2),xyz_post(gptnum,3))
                        call brickshapefun(bne,alpha,eta,beta)
                        call mat_stiff(j,eta_phy)
 
                        uvw_xyz=0.0d0
                        epsnode=0.0d0  
                        signode=0.0d0
                        signode_Le=0.0d0
                        signode_Ls=0.0d0
                        do i15=1,nne
                            nd_var=(i-1)*(nne-1)+i15
                            do j15=1,mexp
                                ndnm=elemcon_replace(j,j15+1,nd_var)
                                if (ndnm.eq.0) cycle
                                c3=(3*sum(maxexp_nodes(1:nd_var-1)))+(3*ndnm) 
  
                                Term_temp(1)=fa(j15)*sf(i15) 
                                Term_temp(2)=f(j15)*dsf(i15)
                                Term_temp(3)=fb(j15)*sf(i15)
                                Git_temp=matmul(jac3D_inv,Term_temp)

                                do i14=1,3
                                    do j14=1,3
                                        uvw_xyz((i14-1)*3+j14)=uvw_xyz((i14-1)*3+j14)+Git_temp(i14)*uvect(c3-3+j14)
                                    enddo
                                enddo
                            enddo
                        enddo            
            
                        !!!!!!!!!!!!!!!!!!!!! strain vector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        epsnode(1)=uvw_xyz(1)                                                                                       !epsxx               
                        epsnode(2)=uvw_xyz(5)                                                                                       !epsyy
                        epsnode(3)=uvw_xyz(9)                                                                                       !epszz
                        epsnode(4)=uvw_xyz(8)+uvw_xyz(6)                                                                            !epsyz
                        epsnode(5)=uvw_xyz(7)+uvw_xyz(3)                                                                            !epsxz
                        epsnode(6)=uvw_xyz(4)+uvw_xyz(2)                                                                            !epsxy

                        !!!!!!!!!!!!!!! compute material stiffness and local cys inverse !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        call matinv3(loc_e,loc_e_inv)
                        call matinv3(loc_struc,loc_struc_inv)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!!!!!!!!!!!!!!!!!!! Global stress vector and tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        signode=matmul(trcmat,epsnode)
                        signode_mat(1,1)=signode(1);    signode_mat(1,2)=signode(6);    signode_mat(1,3)=signode(5)
                        signode_mat(2,1)=signode(6);    signode_mat(2,2)=signode(2);    signode_mat(2,3)=signode(4)
                        signode_mat(3,1)=signode(5);    signode_mat(3,2)=signode(4);    signode_mat(3,3)=signode(3)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!!!!!!!!!!!!!!!!!!!!!!!!!! Local stress tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!!!!!!!!!!!!!!!!!!!!!!!!!! Local csys wrt to material csys!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        signode_Le_mat=matmul(transpose(loc_e_inv),matmul(signode_mat,loc_e_inv))
                        signode_Le(1)=signode_Le_mat(1,1);  signode_Le(2)=signode_Le_mat(2,2);  signode_Le(3)=signode_Le_mat(3,3)
                        signode_Le(4)=signode_Le_mat(2,3);  signode_Le(5)=signode_Le_mat(1,3);  signode_Le(6)=signode_Le_mat(1,2)
                        !!!!!!!!!!!!!!!! Local csys wrt to structure csys !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        signode_Ls_mat=matmul(transpose(loc_struc_inv),matmul(signode_mat,loc_struc_inv))
                        signode_Ls(1)=signode_Ls_mat(1,1);  signode_Ls(2)=signode_Ls_mat(2,2);    signode_Ls(3)=signode_Ls_mat(3,3)
                        signode_Ls(4)=signode_Ls_mat(2,3);  signode_Ls(5)=signode_Ls_mat(1,3);    signode_Ls(6)=signode_Ls_mat(1,2)
          
                        do vec=1,6
                            strn_G_gp(gptnum,vec)=epsnode(vec)
                            strs_G_gp(gptnum,vec)=signode(vec)
                            strs_Le_gp(gptnum,vec)=signode_Le(vec)
                            strs_Ls_gp(gptnum,vec)=signode_Ls(vec)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    
        
        write(14,*)'x y z (gauss points)'
        do i=1,total_gpt
            !if (xyz_post(i,1).eq.0.0d0.and.abs(xyz_post(i,2)-244.8036775).le.1e-6) then
                write(14,*)xyz_post(i,1),xyz_post(i,2),xyz_post(i,3)
            !endif
        enddo
        
        write(14,*)'Syy Syz Szz local structural (gauss points)'
        do i=1,total_gpt
             !if (xyz_post(i,1).eq.0.0d0.and.abs(xyz_post(i,2)-244.8036775).le.1e-6) then
                write(14,*)strs_G_gp(i,2),strs_G_gp(i,4),strs_G_gp(i,3)
             !endif
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    DEALLOCATE (xyz_post,strs_Le_gp,strs_Ls_gp,strn_G_gp,strs_G_gp,gpt_cs,wpt_cs,gpt_bm,wpt_bm)
       
     return
    end subroutine
    
!-----------------------------------------------------------------------------------------------------
!      FUNCTION post_processing2_section FOR post processing across a section element by element wise
!-----------------------------------------------------------------------------------------------------
	subroutine post_processing2_section
    use var_inputdata
    use var_analysis
    use mod_math
    
    implicit none
    integer::k3,ele_cs,k2,i,nne_post,tnodes_post,alpbet_cnt,j,pts
    real*8::xm1,zm1,alpha,beta
    real*8,ALLOCATABLE, DIMENSION (:)::nodeval,alpha_p,beta_p
    real*8,allocatable,dimension(:,:,:)::disp_p_sec,strn_p_G_sec,strs_p_G_sec,strs_p_Le_sec,strs_p_Ls_sec,Git_temp_sec_1
    real*8,allocatable,dimension(:,:)::u_G,eps_G,sig_G,sig_Le,sig_Ls,xyz_post,eps_G_1
	real*8,allocatable,dimension(:,:,:,:)::Git_temp_sec
	integer, ALLOCATABLE, DIMENSION (:)::ele_cs_p
    allocate (alpha_p(4),beta_p(4),ele_cs_p(4))

	be_div_post_p=1
    nne_post=be_div_post_p*(nne-1)+1
    tnodes_post=ne*(nne_post-1)+1
	pts=tpts(2)
    
    allocate (strs_p_Le_sec(tnodes_post,9,ane*pts),strs_p_Ls_sec(tnodes_post,9,ane*pts),strn_p_G_sec(tnodes_post,9,ane*pts),strs_p_G_sec(tnodes_post,9,ane*pts))
    allocate (disp_p_sec(tnodes_post,6,ane*pts))
	allocate (Git_temp_sec(tnodes_post,3,ane*pts,mexp))
    disp_p_sec=0.0d0
	strs_p_Le_sec=0.0d0
	strs_p_Ls_sec=0.0d0
	strn_p_G_sec=0.0d0
	strs_p_G_sec=0.0d0
	Git_temp_sec=0.0d0	
    do ele_cs=1,ane

        allocate (eps_G(tnodes_post,6),u_G(tnodes_post,3),xyz_post(tnodes_post,3),eps_G_1(tnodes_post,6))
        allocate (sig_G(tnodes_post,6),sig_Le(tnodes_post,6),sig_Ls(tnodes_post,6),nodeval(tnodes_post))
		allocate (Git_temp_sec_1(tnodes_post,3,mexp))
        
        do k3=1,pts
			alpha=normcs_tri(k3,1)				!normcs(k3,1) 			!normcs_tri(k3,1)   
			beta=normcs_tri(k3,2)				!normcs(k3,2) 			!normcs_tri(k3,2) 
            xm1=phycs(k3+pts*(ele_cs-1),2)
            zm1=phycs(k3+pts*(ele_cs-1),3)
			
			call post_disp_strs_strn(xm1,zm1,alpha,beta,ele_cs,nne_post,xyz_post,u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval,eps_G_1,Git_temp_sec_1) 

            do k2=1,tnodes_post
                do i=1,3
					disp_p_sec(k2,i,k3+pts*(ele_cs-1))=xyz_post(k2,i)
					strn_p_G_sec(k2,i,k3+pts*(ele_cs-1))=xyz_post(k2,i)
					strs_p_G_sec(k2,i,k3+pts*(ele_cs-1))=xyz_post(k2,i)
					strs_p_Le_sec(k2,i,k3+pts*(ele_cs-1))=xyz_post(k2,i)
					strs_p_Ls_sec(k2,i,k3+pts*(ele_cs-1))=xyz_post(k2,i)
						
                    disp_p_sec(k2,i+3,k3+pts*(ele_cs-1))=u_G(k2,i)
					
					do j=1,mexp
						Git_temp_sec(k2,i,k3+pts*(ele_cs-1),j)=Git_temp_sec_1(k2,i,j)
					enddo
					
				enddo

					
                do i=1,6
                    strn_p_G_sec(k2,i+3,k3+pts*(ele_cs-1))=eps_G(k2,i)                     
                    strs_p_G_sec(k2,i+3,k3+pts*(ele_cs-1))=sig_G(k2,i)
                    strs_p_Le_sec(k2,i+3,k3+pts*(ele_cs-1))=sig_Le(k2,i)
                    strs_p_Ls_sec(k2,i+3,k3+pts*(ele_cs-1))=sig_Ls(k2,i)
                enddo
            enddo 
		enddo   
		
		deallocate (u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval,xyz_post,eps_G_1,Git_temp_sec_1)
	enddo
	
	do k3=1,ane*pts
        write(14,*)k3,disp_p_sec(((tnodes_post-1)*50/100)+1,1,k3),disp_p_sec(((tnodes_post-1)*50/100)+1,3,k3)
	enddo
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gitx-1"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,1,k3,1)
	enddo
	write(14,"(a)")'$EndNodeData'	
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gity-1"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,2,k3,1)
	enddo
	write(14,"(a)")'$EndNodeData'
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gitz-1"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,3,k3,1)
	enddo
	write(14,"(a)")'$EndNodeData'
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gitx-2"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,1,k3,2)
	enddo
	write(14,"(a)")'$EndNodeData'	
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gity-2"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,2,k3,2)
	enddo
	write(14,"(a)")'$EndNodeData'
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gitz-2"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,3,k3,2)
	enddo
	write(14,"(a)")'$EndNodeData'
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gitx-3"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,1,k3,3)
	enddo
	write(14,"(a)")'$EndNodeData'	
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gity-3"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,2,k3,3)
	enddo
	write(14,"(a)")'$EndNodeData'
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Gitz-3"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,Git_temp_sec(((tnodes_post-1)*50/100)+1,3,k3,3)
	enddo
	write(14,"(a)")'$EndNodeData'
	
    write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Ux"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,disp_p_sec(((tnodes_post-1)*50/100)+1,4,k3)
	enddo
	write(14,"(a)")'$EndNodeData'	
		
	
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Uy"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,disp_p_sec(((tnodes_post-1)*50/100)+1,5,k3)
	enddo
	write(14,"(a)")'$EndNodeData'	
		
	write(14,"(a)")'$NodeData'
	write(14,"(a)")'1'
	write(14,"(a)")'"Uz"'
	write(14,"(a)")'1'
	write(14,"(a)")'0.0'
	write(14,"(a)")'3'
	write(14,"(a)")'0'
	write(14,"(a)")'1'
	write(14,"(I5)")ane*pts
	do k3=1,ane*pts
        write(14,*)k3,disp_p_sec(((tnodes_post-1)*50/100)+1,6,k3)
	enddo
	write(14,"(a)")'$EndNodeData'	
	
 !   write(14,*)'$NodeData'
	!write(14,*)'1'
	!write(14,*)'"Sxx"'
	!write(14,*)'1'
	!write(14,*)'0.0'
	!write(14,*)'3'
	!write(14,*)'0'
	!write(14,*)'1'
	!write(14,*)ane*pts
 !   do k3=1,ane*pts
 !       write(14,*)k3,strs_p_G_sec(((tnodes_post-1)*50/100)+1,4,k3)
	!enddo
	!write(14,*)'$EndNodeData'
	!	
	!
 !   write(14,*)'$NodeData'
	!write(14,*)'1'
	!write(14,*)'"Syy"'
	!write(14,*)'1'
	!write(14,*)'0.0'
	!write(14,*)'3'
	!write(14,*)'0'
	!write(14,*)'1'
	!write(14,*)ane*pts
 !   do k3=1,ane*pts
 !       write(14,*)k3,strs_p_G_sec(((tnodes_post-1)*50/100)+1,5,k3)
	!enddo
	!write(14,*)'$EndNodeData'
	!	
	!
 !   write(14,*)'$NodeData'
	!write(14,*)'1'
	!write(14,*)'"Szz"'
	!write(14,*)'1'
	!write(14,*)'0.0'
	!write(14,*)'3'
	!write(14,*)'0'
	!write(14,*)'1'
	!write(14,*)ane*pts
 !   do k3=1,ane*pts
 !       write(14,*)k3,strs_p_G_sec(((tnodes_post-1)*50/100)+1,6,k3)
	!enddo
	!write(14,*)'$EndNodeData'
	!	
 !
	!write(14,*)'$NodeData'
	!write(14,*)'1'
	!write(14,*)'"Eyz"'
	!write(14,*)'1'
	!write(14,*)'0.0'
	!write(14,*)'3'
	!write(14,*)'0'
	!write(14,*)'1'
	!write(14,*)ane*pts
 !   do k3=1,ane*pts
 !       write(14,*)k3,strn_p_G_sec(((tnodes_post-1)*50/100)+1,7,k3)
	!enddo
	!write(14,*)'$EndNodeData'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    DEALLOCATE (disp_p_sec,strs_p_G_sec,strn_p_G_sec,strs_p_Le_sec,strs_p_Ls_sec,Git_temp_sec)
       
     return
	end subroutine
	
!-----------------------------------------------------------------------------------
!      FUNCTION Rvector_mat required for transforming Cauchy's stress in deformed Csys
!-----------------------------------------------------------------------------------
    subroutine Rvector_mat(defgrad,Rvector)
    use mod_math
    
    implicit none
    real*8::defgrad(3,3),Rvector(3,3),FTF(3,3),eigval(3),eigval_mat(3,3),sqrtFTF(3,3),sqrtFTF_inv(3,3),eigvec(3,3),eigvec_inv(3,3)
    real*8,allocatable::work(:)
    integer::lwork,info

    lwork=8                         ! lwork=max(1,3*k-1), k is the dimension of FTF
    allocate(work(lwork))
    FTF=0.0d0
    eigval_mat=0.0d0
    eigvec=0.0d0
    eigvec_inv=0.0d0
    sqrtFTF=0.0d0
    sqrtFTF_inv=0.0d0
    Rvector=0.0d0
    
    FTF=matmul(transpose(defgrad),defgrad)

    call dsyev('V','U',3,FTF,3,eigval,work,lwork,info)
    if(info.ne.0) then
        write(*,*)'check Eigen value solver for Cauchy stress computation'
    endif
    
    eigval_mat(1,1)=sqrt(eigval(1)); eigval_mat(2,2)=sqrt(eigval(2)); eigval_mat(3,3)=sqrt(eigval(3))
    eigvec(:,:)=FTF(:,:)
    call matinv3(eigvec,eigvec_inv)
    
    sqrtFTF=matmul(eigvec,matmul(eigval_mat,eigvec_inv))
    
    call matinv3(sqrtFTF,sqrtFTF_inv)
    Rvector=matmul(defgrad,sqrtFTF_inv)
    
    deallocate(work)
    
    return
    end subroutine
 
    
end module
