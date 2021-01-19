module read_input_data
    contains
      
subroutine read_input   
use var_inputdata
use var_analysis
use mod_math
use mod_gauss
    
    implicit none
!!!!!!!!!!!!!! internal variables !!!!!!!!!!!!!!!!!!!!!
    real*8,allocatable,dimension(:,:)   ::elas_constant

    integer                             ::i,j,k,ij,ik,jk,ijk,ndnm,elbm,count
    character*1                         ::symbol,mark
    character*25                        ::filename_cs,filename_pp
     
    integer                             ::anodes,cstype,maxnode,bmnod_coll,ane_coll
    real*8,ALLOCATABLE,DIMENSION(:)     ::anod,xloc,yloc,zloc
    real*8,ALLOCATABLE,DIMENSION(:,:)   ::elem_rot
    integer,ALLOCATABLE,DIMENSION(:,:)  ::elem1con,elecon,elem2con,elemcon_test
    
    integer                             ::nset,totnldpt
    real*8,ALLOCATABLE,DIMENSION(:)     ::tloadx,tloady,tloadz
    real*8,allocatable,dimension(:,:,:) ::load
    integer,ALLOCATABLE,DIMENSION(:)    ::set,nldpt,loadtype
    
    integer,ALLOCATABLE,DIMENSION(:,:)  ::post_elecon1
    integer                             ::post_elem_max
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    symbol='!'
    
    open(1,file="input/SLE-input_comp.dat")
    open(2,file="input/cs_mesh.msh")
    open(3,file="input/post_points.dat")
    open(4,file="input/load_inp.dat")
    
999 format(A1)
     
!!!!!!!!!!!!!!!!!!!!!!!!!! number of cross-sections !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
99 read(1,999)mark
    if(mark.eq.symbol) go to 99
    read(1,*)analysis_type,fullorpost
     
!!!!!!!!!!!!!!!!!!!!!!!!!! Material !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
101 read(1,999)mark
    if(mark.eq.symbol) go to 101
    read(1,*)nmat,nlayer
    allocate(matin(nlayer),angle(nlayer,2),elas_constant(nmat,9),mat_strength(nmat,9))
     mat_strength=0.0d0
     
102  read(1,999)mark
     if(mark.eq.symbol) go to 102
     do k=1,nlayer
 	    read(1,*)matin(k),angle(k,1),angle(k,2)
     enddo
     
103  read(1,999)mark
     if(mark.eq.symbol) go to 103
     do i=1,nmat
        read(1,*)(elas_constant(i,j),j=1,9),(mat_strength(i,j),j=1,9)
     enddo 
     
     call inp1(elas_constant)
     deallocate(elas_constant)
     
!!!!!!!!!!!!!!!!!!!!!!!!!! Beam - LE inputs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
104  read(1,999)mark
     if(mark.eq.symbol) go to 104
     read(1,*)nexp
     
105  read(1,999)mark
     if(mark.eq.symbol) go to 105
     read(1,*)ne,nne,beam_mesh
     
    SELECT CASE(nne)
        CASE (2)
        ngp=2
        CASE (3)
        ngp=3
        CASE (4)
        ngp=4
        CASE (5)
        ngp=5
    END SELECT
        
     allocate(gpt_bm(ngp),wpt_bm(ngp))
     call gauss_bm(ngp)
     
106  read(1,999)mark
     if(mark.eq.symbol) go to 106
     read(1,*)ylengthi,ylengthf
     
     call inp2
     
107  read(1,999)mark
     if(mark.eq.symbol) go to 107
     read(1,*)phynormesh,be_div_post
     
!!!!!!!!! number of terms in the cross-section expansion and gauss points !!!!!!!!!!!!!     
    if (nexp.le.3) then
        mexp=4*nexp
    else
        mexp=4*nexp+0.5*(nexp-2)*(nexp-3)
    endif
    
    agp_quad=nexp+1
    agp_tri=25
    allocate(gpt_cs(agp_quad),wpt_cs(agp_quad),agpt_cs_t(agp_tri),bgpt_cs_t(agp_tri),wpt_cs_t(agp_tri),agpt_cs_t_test(agp_tri),bgpt_cs_t_test(agp_tri),wpt_cs_t_test(agp_tri))
    call gauss_cs(agp_quad)
    call gauss_cs_tri(agp_tri)
!!!!!!!!!!!!!!!!!!!!!!!!!! Load and BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
108  read(4,999)mark
     if(mark.eq.symbol) go to 108
     read(4,*)nset,totnldpt
     
     allocate (load(9,totnldpt,nset),set(nset),nldpt(nset),tloadx(nset),tloady(nset),tloadz(nset),loadtype(nset))
     allocate(ane_load_num(nset),ane_load(nset,20))
     ane_load=0
      
109  read(4,999)mark
     if(mark.eq.symbol) go to 109
     do i=1,nset
          read(4,*)set(i),nldpt(i),tloadx(i),tloady(i),tloadz(i),loadtype(i),ane_load_num(i),(ane_load(i,j),j=1,ane_load_num(i))
          do k=1,nldpt(i)
            read(4,*)(load(j,k,i),j=1,9)
          enddo
     enddo
     
110  read(4,999)mark
     if(mark.eq.symbol) go to 110
     read(4,*)num_bcinp 
     allocate(node_bc(num_bcinp,3))
     
111  read(4,999)mark
     if(mark.eq.symbol) go to 111
     do i=1,num_bcinp
        read(4,*)node_bc(i,1),node_bc(i,2),node_bc(i,3)   !,node_bc(i,4),node_bc(i,5),node_bc(i,6)
     enddo
     
     allocate(cs_node_bc(num_bcinp,nexp+1))
     cs_node_bc=0
     cs_node_bc(:,1)=node_bc(:,2)
     cs_node_bc(:,2)=node_bc(:,3)
     
!!!!!!!! load input processing !!!!!!!!!    
    call inp_load(nset,nldpt,tloadx,tloady,tloadz,load,loadtype)
     
    deallocate (load,set,nldpt,tloadx,tloady,tloadz,loadtype)  
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!! area mesh input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    read(2,999)mark
    read(2,*)anode,ane
  
    allocate(cscord(anode,3),elemcon(ane,mexp+1),layerno(ane,2),cs_ele_type(ane,2),elemcon_replace(ane,mexp+1,tnodes),elemcon_test(ane,4+1),elem_beam_coll(ane,3))
    allocate(anod(anode),xloc(anode),yloc(anode),zloc(anode))
    allocate (elecon(ane,4+5),elem1con(ane,4+1))
    cscord=0.0d0
    elemcon=0
    elecon=0
    layerno=0
    elemcon_replace=0
    elemcon_test=0
    elem_beam_coll=0
        
    read(2,999)mark
    do i=1,anode
        read(2,*)anod(i),xloc(i),yloc(i),zloc(i)
    enddo
    cscord(1:anode,1)=anod(:)
    cscord(1:anode,2)=xloc(:)
    cscord(1:anode,3)=zloc(:)
    
    read(2,*)mark
    do i=1,ane
        read(2,*)(elecon(i,j),j=1,4+5)
    enddo

    elem1con(:,1)=elecon(:,1)
    elem1con(:,2:4+1)=elecon(:,1+5:4+5)
    layerno(1:ane,1)=elecon(:,1)
    layerno(1:ane,2)=elecon(:,4)
    cs_ele_type(1:ane,1)=elecon(:,1)
    cs_ele_type(1:ane,2)=elecon(:,5)
    call sle_elem_connect(elem1con)
    
    elemcon_test(:,1:5)=elemcon(:,1:5)
    do i=1,tnodes
        elemcon_replace(:,:,i)=elemcon(:,:)
    enddo
    do i=1,ane
        elem_beam_coll(i,1)=i
    enddo
    
    deallocate(elecon,elem1con,xloc,yloc,zloc,anod)
    
    maxnode=maxval(elemcon(1:ane,1:9))                      ! needed for quadratic (L8 and L6 cs connectivity)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! generate L9 mesh in cross-section for tangent and normals !!!!!!!!!!!!!!!!!!!!
    anodes3D=maxnode                                    ! number of nodes in order SL2 
    allocate(cscord3D(anodes3D,3),elemcon3D(ane,9))
    cscord3D=0.0d0
    elemcon3D=0
    call cs_elem_3D
    
    !!!!!!!!!!!!!!!!!!!!!!!! tdof !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     maxexp=maxexp_cs
     ALLOCATE(maxexp_nodes(tnodes))
     maxexp_nodes=0
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    read(2,*)mark
    read(2,*)numcs_collap
    read(2,*)mark
    do i=1,numcs_collap
        read(2,*)ane_coll,bmnod_coll
        elbm=1+(bmnod_coll-1)/(nne-1)
        elem_beam_coll(ane_coll,2)=elbm
        elem_beam_coll(ane_coll,3)=bmnod_coll

        do ij=1,ane
            if (ij.lt.ane_coll) then
                elemcon_test(ij,1+1)=elemcon_replace(ij,1+1,bmnod_coll)
                elemcon_test(ij,2+1)=elemcon_replace(ij,2+1,bmnod_coll)
                elemcon_test(ij,3+1)=elemcon_replace(ij,3+1,bmnod_coll)
                elemcon_test(ij,4+1)=elemcon_replace(ij,4+1,bmnod_coll)
            elseif (ij.eq.ane_coll) then
                elemcon_test(ij,4+1)=elemcon_replace(ij,1+1,bmnod_coll)
                elemcon_test(ij,3+1)=elemcon_replace(ij,2+1,bmnod_coll)
            else
                count=0
                do ijk=1,ij-1
                    if (elem_beam_coll(ijk,2).gt.0) then
                        count=count+1
                    endif
                enddo
                if (elem_beam_coll(ij,2).gt.0.and.elem_beam_coll(ij,3).lt.bmnod_coll) then
                    elemcon_test(ij,1+1)=elemcon(ij-count,1+1)
                    elemcon_test(ij,2+1)=elemcon(ij-count,2+1)
                    elemcon_test(ij,3+1)=elemcon(ij-count,2+1)
                    elemcon_test(ij,4+1)=elemcon(ij-count,1+1)
                else
                    elemcon_test(ij,1+1)=elemcon(ij-count,1+1)
                    elemcon_test(ij,2+1)=elemcon(ij-count,2+1)
                    elemcon_test(ij,3+1)=elemcon(ij-count,3+1)
                    elemcon_test(ij,4+1)=elemcon(ij-count,4+1)
                endif
            endif
        enddo
        call sle_elem_connect_replace(elemcon_test,bmnod_coll,k)
        do ik=bmnod_coll+1,tnodes
            elemcon_replace(:,:,ik)=elemcon_replace(:,:,bmnod_coll)
        enddo
    enddo
     
    !allocate(elem_collap(tnodes,ane+2,4+1))
    !elem_collap=0
    !do i=1,tnodes
    !    elem_collap(i,1,1)=i
    !enddo
    
    !read(2,*)mark
    !read(2,*)numcs_collap                                           ! total number of cs elements collapsed
    
    !read(2,*)mark
    !do i=1,numcs_collap
    !    read(2,*)bmnod_coll,elem_collap(bmnod_coll,2,1)             ! beam node at which element collapsed, number of cs elements collaped at that beam node
    !    elbm=1+(bmnod_coll-1)/(nne-1)                               ! beam element at which the cs element is collapsed
    !    do k=1,elem_collap(bmnod_coll,2,1)
    !        read(2,*)elem_collap(bmnod_coll,2+k,1),(elem_collap(bmnod_coll,2+k,j+1),j=1,4)   ! element cs number which is collapsed, cs nodes which are collapsed
    !        elem_beam_coll(elem_collap(bmnod_coll,2+k,1),2)=elbm
    !        do ij=elem_collap(bmnod_coll,2+k,1),ane
    !            if (ij.eq.elem_collap(bmnod_coll,2+k,1)) then
    !                do jk=1,4
    !                    ndnm=elemcon_test(ij,jk+1)
    !                    if (ndnm.eq.elem_collap(bmnod_coll,2+k,3)) then
    !                        elemcon_test(ij,jk+1)=elem_collap(bmnod_coll,2+k,2)
    !                    elseif (ndnm.eq.elem_collap(bmnod_coll,2+k,5)) then
    !                        elemcon_test(ij,jk+1)=elem_collap(bmnod_coll,2+k,4) 
    !                    endif
    !                enddo
    !            elseif (elem_beam_coll(ij,2).gt.0.and.elem_beam_coll(ij,2).lt.elbm) then
    !                if (ij-1.eq.elem_collap(bmnod_coll,2+k,1)) then
    !                    elemcon_test(ij,2:5)=elemcon_test(ij-1,2:5)
    !                endif
    !            else
    !                count=0
    !                do ijk=1,ij
    !                    if (elem_beam_coll(ijk,2).gt.0) then
    !                        count=count+1
    !                    endif
    !                enddo
    !                elemcon_test(ij,2:5)=elemcon(ij-count,2:5)
    !            endif
    !        enddo
    !        call sle_elem_connect_replace(elemcon_test,bmnod_coll,k)
    !    enddo
    !    elemcon_test(:,1:5)=elemcon_replace(:,1:5,bmnod_coll)
    !    do ik=bmnod_coll+1,tnodes
    !        elemcon_replace(:,:,ik)=elemcon_replace(:,:,bmnod_coll)
    !    enddo
    !enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     do i=1,ne
         do j=1,nne
            ijk=(i-1)*(nne-1)+j
            maxexp_nodes(ijk)=maxval(elemcon_replace(1:ane,2:mexp+1,ijk))
            !if (elem_collap(ijk,2,1).eq.0) then             ! no element collapsed at the beam node 'ijk'
            !    maxexp_nodes(ijk)=maxexp
            !else
            !    maxexp_nodes(ijk)=maxexp-mexp*elem_collap(ijk,2,1)+(nexp+1)*elem_collap(ijk,2,1)
            !endif            
         enddo
     enddo
     
     tdof=3*sum(maxexp_nodes(1:tnodes))
	 band_width=3*maxexp*nne
	 tdof1=3*sum(maxexp_nodes(2:tnodes))
     write(11,*)'tdof',tdof
     
!!!!!!!!!!!!!! post-processing input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
    read(3,999)mark
    read(3,*)npts,post_elem

    allocate(postpoint(npts,4),post_elecon(post_elem,5))
    allocate(post_elecon1(post_elem,9))
    
    postpoint=0.0
    post_elecon=0
     
112     read(3,999)mark
    if(mark.eq.symbol) go to 112
    do i=1,npts
        read(3,*)postpoint(i,1),postpoint(i,2),postpoint(i,3),postpoint(i,4)
    enddo
     
113     read(3,999)mark
    if(mark.eq.symbol) go to 113
    do i=1,post_elem
        read(3,*)(post_elecon1(i,j),j=1,9)
    enddo    
    
    post_elecon(:,1)=post_elecon1(:,1)
    post_elecon(:,2:5)=post_elecon1(:,6:9)
    
    deallocate(post_elecon1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    close(1)
    close(2)
    close(3)
    close(4)
    
end subroutine


!----------------------------------------------------------------
!      subroutine inp2 FOR beam NODE COORDINATES
!----------------------------------------------------------------
subroutine inp2
    use var_inputdata
    implicit none
    Integer::ii,jj,el,tnds,kk
    real*8::tl,dbn,lpe,pi,y1,y2
    pi=4.0*ATAN(1.0)
        
    tnodes=((nne-1)*ne)+1
    allocate(y(tnodes))
        
    y=0.0d0 
    y(1)=ylengthi
        
    SELECT CASE (beam_mesh)
    CASE ('CHEB','Cheb','cheb')
        tl=ylengthf-ylengthi
        do ii=2,tnodes-1
            y(ii)=(tl/2.0)-(0.5d0*tl*cos(pi*(ii-0.5d0)/tnodes))
        enddo
        y(tnodes)=ylengthf
            
    CASE ('UNIF','Unif','unif')
        tl=ylengthf-ylengthi
        lpe=tl/ne            
        dbn=lpe/(nne-1)    
        do ii=2,tnodes-1
            y(ii)=y(ii-1)+dbn
        enddo
        y(tnodes)=ylengthf
        
    CASE ('MANU','Manu','manu')
        open(5,file="input/beam_yinp.dat")
        do ii=1,tnodes
            read(5,*)y(ii)
        enddo
        close(5)
        
    END SELECT

    return
end subroutine
    
    
!------------------------------------------------------------------------------
!      Definition of subroutine inp1 Material Properties
!------------------------------------------------------------------------------
    subroutine inp1(elas_constant)
    use var_inputdata
    use mod_math
    Implicit none
    integer::i
    real*8::E1,E2,E3,G23,G13,G12,nu12,nu13,nu23,delta,nu21,nu31,nu32
    real*8,allocatable,dimension(:,:)::elas_constant
    allocate(cmat(6,6,nmat))
    cmat=0.0
    
    do i=1,nmat
        E1=elas_constant(i,1)  
        E2=elas_constant(i,2) 
        E3=elas_constant(i,3) 
        G23=elas_constant(i,4) 
        G13=elas_constant(i,5) 
        G12=elas_constant(i,6) 
        nu12=elas_constant(i,7) 
        nu13=elas_constant(i,8) 
        nu23=elas_constant(i,9) 

        nu21=nu12*E2/E1
	    nu31=nu13*E3/E1
	    nu32=nu23*E3/E2
    
	    delta=(1.0-(nu12*nu21)-(nu23*nu32)-(nu31*nu13)-(2.*nu21*nu13*nu32))
        cmat(1,1,i) = (1.0 - nu23*nu32)*E1/delta
        cmat(2,2,i) = (1.0 - nu13*nu31)*E2/delta
        cmat(3,3,i) = (1.0 - nu12*nu21)*E3/delta
        cmat(1,2,i) = (nu12 + nu32*nu13)*E2/delta
        cmat(2,1,i) = (nu21 + nu23*nu31)*E1/delta
        cmat(1,3,i) = (nu13 + nu12*nu23)*E3/delta
        cmat(3,1,i) = (nu31 + nu21*nu32)*E1/delta
        cmat(2,3,i) = (nu23 + nu21*nu13)*E3/delta
        cmat(3,2,i) = (nu32 + nu12*nu31)*E2/delta
    
        if (E1.eq.E2.and.E2.eq.E3) then
            cmat(4,4,i) = (E1/(2.0d0*(1+nu12)))
            cmat(5,5,i) = (E1/(2.0d0*(1+nu12)))
            cmat(6,6,i) = (E1/(2.0d0*(1+nu12)))
        else  
            cmat(4,4,i) = G23
            cmat(5,5,i) = G13
            cmat(6,6,i) = G12
        endif
    enddo
    
    
    return
    end subroutine

    
!------------------------------------------------------------------------------
!      Definition of subroutine mat_stiff Material Properties
!------------------------------------------------------------------------------
    subroutine mat_stiff(ele,yloc)
    use var_inputdata
    use var_analysis
    use mod_math
    Implicit none
    integer::mat,lay,ele,i30,j30
    real*8::yloc,ang,x1,x2,x3,x4,x5,x6,x7,x8,x9
    real*8,allocatable,dimension(:)::V,N,P,T,Tcos,Tsin
    real*8,allocatable,dimension(:,:)::cmat1,Qrot,lmn,lmn1
    allocate(cmat1(6,6),Qrot(6,6),lmn(5,5),lmn1(3,3))
    allocate(V(3),N(3),T(3),P(3),Tcos(3),Tsin(3))
    
    lay=layerno(ele,2)
    
    !x1=32.0d0; x2=34.0d0; x3=36.0d0   
    !x4=38.0d0; x5=40.0d0; x6=42.0d0 
    !x7=44.0d0; x8=46.0d0; x9=48.0d0
    !if (ele.eq.5.and.yloc-x2.gt.1e-8.and.yloc-x3.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.6.and.yloc-x6.gt.1e-8.and.yloc-x7.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.7.and.yloc-x8.gt.1e-8.and.yloc-x9.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.8.and.yloc-x4.gt.1e-8.and.yloc-x5.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.11.and.yloc-x3.gt.1e-8.and.yloc-x4.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.12.and.yloc-x7.gt.1e-8.and.yloc-x8.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.13.and.yloc-x5.gt.1e-8.and.yloc-x6.le.1e-8) then
    !    lay=17
    !endif
    !if (ele.eq.14.and.yloc-x1.gt.1e-8.and.yloc-x2.le.1e-8) then
    !    lay=17
    !endif
    
    mat=matin(lay)
    ang=(2*(angle(lay,1)-angle(lay,2))*abs(yloc-0.5d0*(y(tnodes)-y(1)))/(y(tnodes)-y(1)))+angle(lay,2)
    cmat1(:,:)=cmat(:,:,mat)
    
    !!!!! fibre direction for straight beams
    V(1)=0.0d0; V(2)=1.0d0; V(3)=0.0d0    
    
    !!!!! fibre direction for curved beams (fibres should follow the structure path)
    !V(1)=jac3D(2,1)/(sqrt(jac3D(2,1)**2+jac3D(2,2)**2+jac3D(2,3)**2))
    !V(2)=jac3D(2,2)/(sqrt(jac3D(2,1)**2+jac3D(2,2)**2+jac3D(2,3)**2))
    !V(3)=jac3D(2,3)/(sqrt(jac3D(2,1)**2+jac3D(2,2)**2+jac3D(2,3)**2))
    
    ! normal direction_temporary (check for an orthotropic case or for laminates, if tapered... in that case better to have a 2D Jacobain for cross-section)
    N(1)=jac3D(3,1)/(sqrt(jac3D(3,1)**2+jac3D(3,2)**2+jac3D(3,3)**2))
    N(2)=jac3D(3,2)/(sqrt(jac3D(3,1)**2+jac3D(3,2)**2+jac3D(3,3)**2))
    N(3)=jac3D(3,3)/(sqrt(jac3D(3,1)**2+jac3D(3,2)**2+jac3D(3,3)**2))
    
    ! tangent direction
    P(1)=V(2)*N(3)-N(2)*V(3); P(2)=V(3)*N(1)-N(3)*V(1); P(3)=V(1)*N(2)-N(1)*V(2)
    
    ! Normal direction
    N(1)=P(2)*V(3)-V(2)*P(3); N(2)=P(3)*V(1)-V(3)*P(1); N(3)=P(1)*V(2)-V(1)*P(2)
    
    ! axial direction
    T(1)=V(1); T(2)=V(2); T(3)=V(3)
    
    ! considering fibre orientation in the xy plane
    Tsin(1)=T(1)*sind(-ang); Tsin(2)=T(2)*sind(-ang); Tsin(3)=T(3)*sind(-ang)
    Tcos(1)=T(1)*cosd(-ang); Tcos(2)=T(2)*cosd(-ang); Tcos(3)=T(3)*cosd(-ang)
    
!!!!!!!!!!!!! ang is the fibre angle with respect to y axis (beam axis) (+ve input means anti-cs) !!!!!!!!!!!!!!!!!!!!!!!!
    !!! loc_e is the coordinate system such that y-axis is always along the fibre direction
    !!! It can be called as a material csys, Z is normal, y is along fibre, x is tangent. 
    !!! loc_e will be diff for laminates with diff fibre orientation
    !!! loc_struc is another local csys which is only depends on the structure orientation, which is used in post-processing
    !!!!!!!!!!!!!!!!!! local x, y and z coordinate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    loc_e(1,1)=(Tcos(2)*N(3)-Tcos(3)*N(2))-Tsin(1)
    loc_e(1,2)=(Tcos(3)*N(1)-Tcos(1)*N(3))-Tsin(2)
    loc_e(1,3)=(Tcos(1)*N(2)-Tcos(2)*N(1))-Tsin(3)
    
    loc_e(2,1)=Tsin(2)*N(3)-Tsin(3)*N(2)+Tcos(1)
    loc_e(2,2)=Tsin(3)*N(1)-Tsin(1)*N(3)+Tcos(2)
    loc_e(2,3)=Tsin(1)*N(2)-Tsin(2)*N(1)+Tcos(3) 
    
    loc_e(3,1)=N(1)
    loc_e(3,2)=N(2)
    loc_e(3,3)=N(3)
    
    loc_struc(1,1)=jac3D(1,1)/(sqrt(jac3D(1,1)**2+jac3D(1,2)**2+jac3D(1,3)**2))
    loc_struc(1,2)=jac3D(1,2)/(sqrt(jac3D(1,1)**2+jac3D(1,2)**2+jac3D(1,3)**2))
    loc_struc(1,3)=jac3D(1,3)/(sqrt(jac3D(1,1)**2+jac3D(1,2)**2+jac3D(1,3)**2))
    
    loc_struc(2,1)=jac3D(2,1)/(sqrt(jac3D(2,1)**2+jac3D(2,2)**2+jac3D(2,3)**2))
    loc_struc(2,2)=jac3D(2,2)/(sqrt(jac3D(2,1)**2+jac3D(2,2)**2+jac3D(2,3)**2))
    loc_struc(2,3)=jac3D(2,3)/(sqrt(jac3D(2,1)**2+jac3D(2,2)**2+jac3D(2,3)**2))
    
    loc_struc(3,1)=jac3D(3,1)/(sqrt(jac3D(3,1)**2+jac3D(3,2)**2+jac3D(3,3)**2))
    loc_struc(3,2)=jac3D(3,2)/(sqrt(jac3D(3,1)**2+jac3D(3,2)**2+jac3D(3,3)**2))
    loc_struc(3,3)=jac3D(3,3)/(sqrt(jac3D(3,1)**2+jac3D(3,2)**2+jac3D(3,3)**2))
    
    !!!!!!!!!!!!!!!!!! global x, y and z coordinate !!!!!!!!!!!!!!!!!!!!!
    glb_e(1,1)=1.0d0; glb_e(1,2)=0.0d0; glb_e(1,3)=0.0d0
    glb_e(2,1)=0.0d0; glb_e(2,2)=1.0d0; glb_e(2,3)=0.0d0
    glb_e(3,1)=0.0d0; glb_e(3,2)=0.0d0; glb_e(3,3)=1.0d0
    
    lmn1=matmul(loc_e,transpose(glb_e))
    lmn(1:3,1:3)=lmn1(:,:)
    lmn(1:3,4)=lmn1(1:3,1)
    lmn(1:3,5)=lmn1(1:3,2)
    lmn(4:5,1:5)=lmn(1:2,1:5)
    
    do i30=1,3
        do j30=1,3
            Qrot(i30,j30)=lmn(i30,j30)*lmn(i30,j30)
        enddo
    enddo
   
    do i30=1,3
        do j30=4,6
            Qrot(i30,j30)=lmn(i30,j30-2)*lmn(i30,j30-1)
        enddo
    enddo
   
    do i30=4,6
        do j30=1,3
            Qrot(i30,j30)=2*lmn(i30-2,j30)*lmn(i30-1,j30)
        enddo
    enddo
    
    do i30=4,6
        do j30=4,6
            Qrot(i30,j30)=lmn(i30-2,j30-2)*lmn(i30-1,j30-1)+lmn(i30-1,j30-2)*lmn(i30-2,j30-1)
        enddo
    enddo
        
    trcmat=matmul(transpose(Qrot),matmul(cmat1,Qrot))
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Plane strain condition!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !trcmat(1,2:6)=0.0d0
    !trcmat(5,1:4)=0.0d0
    !trcmat(5,6)=0.0d0
    !trcmat(6,1:5)=0.0d0
    !trcmat(2:6,1)=0.0d0
    !trcmat(1:4,5)=0.0d0
    !trcmat(6,5)=0.0d0
    !trcmat(1:5,6)=0.0d0
    !
    !trcmat(1:2,6)=0.0d0
    !trcmat(6,1:2)=0.0d0
    !trcmat(4,5)=0.0d0
    !trcmat(5,4)=0.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    
   ! write(11,*)ang,ele,mat,lay
   ! call PRINTF(trcmat,6,6)

    
    
    deallocate(cmat1,Qrot,lmn,lmn1,V,N,T,P,Tcos,Tsin)
    return
    end subroutine
    
!-----------------------------------------------------------------------------------------
!      FUNCTION sle_elem_connect FOR element connectivity for SLE
!-----------------------------------------------------------------------------------------
	subroutine sle_elem_connect(elem1con)
    use var_inputdata
    use mod_math
        implicit none 
        integer::l,kk,ii,jj,max,m1,m2,i1,j1,i,nodnm,mval,mval1(1),ij
        integer,allocatable,dimension(:,:)::connect,nodarray,elem1con,elem3con
        integer,allocatable,dimension(:,:)::elem2con
        real*8,allocatable::x(:),z(:),xz(:)
        real*8::cross
        
        allocate(connect(ane,mexp),nodarray(anode,anode),x(4),z(4),elem3con(ane,4),xz(4),elem2con(ane,4+1))
        connect=0
        elem2con=0
        elem3con=0
!!!!!!!!!!!!!!!!!!!!!!!!! anti-clockwise connectivity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        do jj=1,ane
            do ii=1,cs_ele_type(jj,2)           ! cs_ele_type (3 or 4; tri or quad) 
                nodnm=elem1con(jj,ii+1)
                x(ii)=cscord(nodnm,2)
                z(ii)=cscord(nodnm,3)
                xz(ii)=x(ii)+z(ii)
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!! bottom left corner to be the first node !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !mval1=minloc(xz)
            !mval=mval1(1)
            !do kk=1,cs_ele_type(jj,2)
            !    if (mval.le.cs_ele_type(jj,2)) then
            !        elem3con(jj,kk)=elem1con(jj,mval+1)
            !        mval=mval+1
            !    else
            !        elem3con(jj,kk)=elem1con(jj,mval+1-cs_ele_type(jj,2))
            !        mval=mval+1
            !    endif
            !enddo
            !elem1con(jj,2:5)=elem3con(jj,1:4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            i=1
1234		cross=((x(i+1)-x(i))*(z(i+2)-z(i+1)))-((x(i+2)-x(i+1))*(z(i+1)-z(i)))
            if (cross.eq.0.0) then
                i=i+1
                goto 1234
            
                ! clockwise
            else if (cross.lt.0.0.and.cs_ele_type(jj,2).eq.4) then
                    elem2con(jj,1)=elem1con(jj,1)
                    elem2con(jj,2)=elem1con(jj,2)
                    elem2con(jj,3)=elem1con(jj,5)
                    elem2con(jj,4)=elem1con(jj,4)
                    elem2con(jj,5)=elem1con(jj,3)
            else if (cross.lt.0.0.and.cs_ele_type(jj,2).eq.3) then
                    elem2con(jj,1)=elem1con(jj,1)
                    elem2con(jj,2)=elem1con(jj,2)
                    elem2con(jj,3)=elem1con(jj,4)
                    elem2con(jj,4)=elem1con(jj,3)
                    elem2con(jj,5)=elem1con(jj,5)
                ! anti-clockwise
            else if (cross.gt.0.0) then                             
                elem2con(jj,:)=elem1con(jj,:)
            endif
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do l=1,nexp
           nodarray=0
           if (l.eq.1) then
                do kk=1,ane
                    do ii=1,cs_ele_type(kk,2)
                        connect(kk,ii)=elem2con(kk,ii+1)
                    enddo
                enddo
                max=maxval(connect)
                
           elseif (l.eq.2.or.l.eq.3) then
                do kk=1,ane
                    if (cs_ele_type(kk,2).eq.4) then            ! quad element
                        m1=4*(l-1)
                        m2=4*l
                    elseif (cs_ele_type(kk,2).eq.3) then        ! tri element
                        !m1=3*(l-1)
                        !m2=0.5*(l+1)*(l+2)
                        m1=4*(l-1)
                        m2=4*l
                    endif
                    
                    do ii=1,cs_ele_type(kk,2)
                       i1=elem2con(kk,ii+1)
                       if (ii.le.cs_ele_type(kk,2)-1) then
                        j1=elem2con(kk,ii+2)
                       else
                        j1=elem2con(kk,2)
                       endif 
               
                        if (nodarray(i1,j1).eq.0.and.i1.ne.j1) then
                            connect(kk,m1+ii)=max+1
                            nodarray(i1,j1)=connect(kk,m1+ii)
                            nodarray(j1,i1)=nodarray(i1,j1)
                            max=max+1
                        else
                            connect(kk,m1+ii)=nodarray(i1,j1)
                        endif
                        
                        do ij=1,num_bcinp
                            if (cs_node_bc(ij,1).eq.i1.and.cs_node_bc(ij,2).eq.j1.or.cs_node_bc(ij,1).eq.j1.and.cs_node_bc(ij,2).eq.i1) then
                                cs_node_bc(ij,l+1)=nodarray(i1,j1)
                                goto 999
                            endif
999                     enddo
                    enddo
                    if (cs_ele_type(kk,2).eq.3.and.l.eq.3) then        ! tri element
                        do jj=m1+4,m2
                            connect(kk,jj)=max+1
                            max=max+1
                        enddo
                    endif
                enddo
           else
               do kk=1,ane
                   if (cs_ele_type(kk,2).eq.4) then            ! quad element
                        m1=4*(l-1)+0.5*(l-3)*(l-4)
                        m2=4*l+0.5*(l-2)*(l-3)
                   elseif (cs_ele_type(kk,2).eq.3) then        ! tri element
                        !m1=0.5*l*(l+1)
                        !m2=0.5*(l+1)*(l+2)
                        m1=4*(l-1)+0.5*(l-3)*(l-4)
                        m2=4*l+0.5*(l-2)*(l-3)
                    endif
                    do ii=1,cs_ele_type(kk,2)
                       i1=elem2con(kk,ii+1)
                       if (ii.le.cs_ele_type(kk,2)-1) then
                        j1=elem2con(kk,ii+2)
                       else
                        j1=elem2con(kk,2)
                       endif 
               
                        if (nodarray(i1,j1).eq.0.and.i1.ne.j1) then
                            connect(kk,m1+ii)=max+1
                            nodarray(i1,j1)=connect(kk,m1+ii)
                            nodarray(j1,i1)=nodarray(i1,j1)
                            max=max+1
                        else
                            connect(kk,m1+ii)=nodarray(i1,j1)
                        endif
                        
                        do ij=1,num_bcinp
                            if (cs_node_bc(ij,1).eq.i1.and.cs_node_bc(ij,2).eq.j1.or.cs_node_bc(ij,1).eq.j1.and.cs_node_bc(ij,2).eq.i1) then
                                cs_node_bc(ij,l+1)=nodarray(i1,j1)
                                goto 998
                            endif
998                     enddo
                    enddo
                    
                    do jj=m1+cs_ele_type(kk,2)+1,m2
                        connect(kk,jj)=max+1
                        max=max+1
                    enddo
               enddo
           endif
       enddo
       maxexp_cs=max
       elemcon(1:ane,1)=elem2con(:,1)
       elemcon(1:ane,2:mexp+1)=connect(:,:)
           
       deallocate(connect,nodarray,x,z,elem3con,xz,elem2con)
      return 
    end subroutine

    
!-----------------------------------------------------------------------------------------
!      FUNCTION sle_elem_connect FOR element connectivity for SLE (replacement for ply drops)
!-----------------------------------------------------------------------------------------
	subroutine sle_elem_connect_replace(elemcon_test,bmnod_coll,k)
    use var_inputdata
    use mod_math
        implicit none 
        integer::l,kk,ii,jj,max,m1,m2,i1,j1,i,nodnm,mval,mval1(1),ij,bmnod_coll,k
        integer,allocatable,dimension(:,:)::connect,nodarray,elemcon_test
        integer,allocatable,dimension(:,:)::elem2con
        
        allocate(connect(ane,mexp),nodarray(anode,anode),elem2con(ane,4+1))
        connect=0
!!!!!!!!!!!!!!!!!!!!!!!!! anti-clockwise connectivity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        do jj=1,ane
            elem2con(jj,:)=elemcon_test(jj,:)
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do l=1,nexp
           nodarray=0
           if (l.eq.1) then
                do kk=1,ane
                    do ii=1,4
                        connect(kk,ii)=elem2con(kk,ii+1)
                    enddo
                enddo
                max=maxval(connect)
                
           elseif (l.eq.2.or.l.eq.3) then
               m1=4*(l-1)
               m2=4*l
        
                do kk=1,ane
                    do ii=1,4
                        i1=elem2con(kk,ii+1)
                       if (ii.le.3) then
                        j1=elem2con(kk,ii+2)
                       else
                        j1=elem2con(kk,2)
                       endif 
               
                        if (nodarray(i1,j1).eq.0.and.i1.ne.j1) then
                            connect(kk,m1+ii)=max+1
                            nodarray(i1,j1)=connect(kk,m1+ii)
                            nodarray(j1,i1)=nodarray(i1,j1)
                            max=max+1
                        else
                            connect(kk,m1+ii)=nodarray(i1,j1)
                        endif
                    enddo
                enddo
           else
               m1=4*(l-1)+0.5*(l-3)*(l-4)
               m2=4*l+0.5*(l-2)*(l-3)
               
               do kk=1,ane
                    do ii=1,4
                        i1=elem2con(kk,ii+1)
                       if (ii.le.3) then
                        j1=elem2con(kk,ii+2)
                       else
                        j1=elem2con(kk,2)
                       endif 
               
                        if (nodarray(i1,j1).eq.0.and.i1.ne.j1) then
                            connect(kk,m1+ii)=max+1
                            nodarray(i1,j1)=connect(kk,m1+ii)
                            nodarray(j1,i1)=nodarray(i1,j1)
                            max=max+1
                        else
                            connect(kk,m1+ii)=nodarray(i1,j1)
                        endif
                    enddo
                    
                    do jj=m1+5,m2
                        if (elem_beam_coll(kk,2).ne.0.and.elem_beam_coll(kk,3).le.bmnod_coll) then
                            connect(kk,jj)=0
                        else
                            connect(kk,jj)=max+1
                            max=max+1
                        endif
                    enddo
               enddo
           endif
       enddo
       elemcon_replace(1:ane,1,bmnod_coll)=elem2con(:,1)
       elemcon_replace(1:ane,2:mexp+1,bmnod_coll)=connect(:,:)
           
       deallocate(connect,nodarray,elem2con)
      return 
    end subroutine
    
!-----------------------------------------------------------------------------------------
!      FUNCTION cs_elem_3D FOR 3D brick element cross-section coordinates
!-----------------------------------------------------------------------------------------
	subroutine cs_elem_3D
    use var_inputdata
        implicit none 
        integer::i,m
        
        elemcon3D(1:ane,1:9)=elemcon(1:ane,1:9)
        cscord3D(1:anode,1)=cscord(1:anode,2)
        cscord3D(1:anode,3)=cscord(1:anode,3)
        
        do i=1,ane
            if (cs_ele_type(i,2).eq.4) then            ! quad element
                cscord3D(elemcon3D(i,6),1:3)=0.5d0*(cscord3D(elemcon3D(i,2),1:3)+cscord3D(elemcon3D(i,3),1:3))  
                cscord3D(elemcon3D(i,7),1:3)=0.5d0*(cscord3D(elemcon3D(i,3),1:3)+cscord3D(elemcon3D(i,4),1:3))
                cscord3D(elemcon3D(i,8),1:3)=0.5d0*(cscord3D(elemcon3D(i,4),1:3)+cscord3D(elemcon3D(i,5),1:3))
                cscord3D(elemcon3D(i,9),1:3)=0.5d0*(cscord3D(elemcon3D(i,5),1:3)+cscord3D(elemcon3D(i,2),1:3))
            elseif (cs_ele_type(i,2).eq.3) then        ! tri element
                cscord3D(elemcon3D(i,6),1:3)=0.5d0*(cscord3D(elemcon3D(i,2),1:3)+cscord3D(elemcon3D(i,3),1:3))  
                cscord3D(elemcon3D(i,7),1:3)=0.5d0*(cscord3D(elemcon3D(i,3),1:3)+cscord3D(elemcon3D(i,4),1:3))
                cscord3D(elemcon3D(i,8),1:3)=0.5d0*(cscord3D(elemcon3D(i,4),1:3)+cscord3D(elemcon3D(i,2),1:3))
            endif
        enddo

      return 
    end subroutine

!-----------------------------------------------------------------------------------------------------
!      FUNCTION componentwise_cs_pos FOR element connectivity for SLE at different CS interface
!-----------------------------------------------------------------------------------------------------
!	subroutine componentwise_cs_pos
!    use var_inputdata
!        implicit none 
!        integer::i,j,k,ij,cs1,cs2,el1,el2,el3,el4,pos_el1,pos_el2,pos_el3,pos_el4,beam_node,counter,vir_ele
!        real*8 ::xcs,zcs
!        
!        allocate(cs_pos(maxexp,2,numcs))
!            cs_pos=0
!            
!        do i=1,be_div-1
!            counter=1
!            beam_node=be_cs_type(i,2)*(nne-1)+1
!            cs1=be_cs_type(i,3)
!            cs2=be_cs_type(i+1,3)
!            
!            do j=1,anode(cs2)
!                cs_pos(j,1,cs2)=j
!                xcs=cscord(j,2,cs2)
!                zcs=cscord(j,3,cs2)
!                do ij=1,anode(cs1)
!                    if (xcs.eq.cscord(ij,2,cs1).and.zcs.eq.cscord(ij,3,cs1)) then
!                        cs_pos(j,2,cs2)=ij
!                        goto 2001
!                    elseif (ij.eq.anode(cs1)) then
!                        cs_pos(j,2,cs2)=maxexp_nodes(beam_node)+counter
!                        counter=counter+1
!                    endif
!                enddo
!2001       enddo
!            
!            do j=1,ane(cs2)
!                el1=elemcon(j,2,cs2)               
!                el2=elemcon(j,3,cs2)               
!                el3=elemcon(j,4,cs2)
!                el4=elemcon(j,5,cs2)
!                
!                pos_el1=cs_pos(el1,2,cs2)
!                pos_el2=cs_pos(el2,2,cs2)
!                pos_el3=cs_pos(el3,2,cs2)
!                pos_el4=cs_pos(el4,2,cs2)
!                
!                 do ij=1,ane(cs1)
!                    if (pos_el1.eq.elemcon(ij,2,cs1).and.pos_el2.eq.elemcon(ij,3,cs1).and.pos_el3.eq.elemcon(ij,4,cs1).and.pos_el4.eq.elemcon(ij,5,cs1)) then
!                        do k=6,mexp+1
!                            vir_ele=elemcon(j,k,cs2)
!                            cs_pos(vir_ele,1,cs2)=vir_ele
!                            cs_pos(vir_ele,2,cs2)=elemcon(ij,k,cs1)
!                        enddo
!                        goto 2002
!                    endif
!                 enddo
!2002        enddo
!            
!            do j=1,maxexp_nodes(beam_node)
!                if (cs_pos(j,2,cs2).eq.0) then
!                    cs_pos(j,1,cs2)=j
!                    cs_pos(j,2,cs2)=maxexp_nodes(beam_node)+counter
!                    counter=counter+1
!                endif
!            enddo
!            
!            write(11,*)'cs_pos cs1',cs1
!            do j=1,maxexp
!                write(11,*)cs_pos(j,1,cs1),cs_pos(j,2,cs1)
!            enddo
!        
!            write(11,*)'cs_pos cs2',cs2
!            do j=1,maxexp
!                write(11,*)cs_pos(j,1,cs2),cs_pos(j,2,cs2)
!            enddo
!        
!        enddo
!        
!
!      return 
!    end subroutine

!-----------------------------------------------------------------------------------------
!      FUNCTION inp_load FOR creating the input file for load
!-----------------------------------------------------------------------------------------    
    subroutine inp_load(nset,nldpt,tloadx,tloady,tloadz,load,loadtype)
    use var_inputdata
        implicit none
        integer::ii,jj,ii1,ii2,ii3,nset,set,tsetpt(nset),ndx,ndy,ndz,cnt
        real*8::xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz
        real*8,allocatable,dimension(:)::tloadx,tloady,tloadz
        real*8,allocatable,dimension(:,:,:)::load
        integer,allocatable,dimension(:)::nldpt,loadtype
        
        np=0
        tsetpt=0
        do jj=1,nset
            if (loadtype(jj).eq.2) then
                load(3,:,jj)=1
                load(9,:,jj)=1
            endif
            
            do ii=1,nldpt(jj)                
                np=np+load(3,ii,jj)*load(6,ii,jj)*load(9,ii,jj)
                tsetpt(jj)=tsetpt(jj)+load(3,ii,jj)*load(6,ii,jj)*load(9,ii,jj)
            enddo
        enddo
        
        allocate(loadinp(7,np))
            
        cnt=1
        do jj=1,nset
            do ii=1,nldpt(jj)
                xmin=load(1,ii,jj)
                xmax=load(2,ii,jj)
                ymin=load(4,ii,jj)
                ymax=load(5,ii,jj)
                zmin=load(7,ii,jj)
                zmax=load(8,ii,jj)
                
                ndx=load(3,ii,jj)
                ndy=load(6,ii,jj)
                ndz=load(9,ii,jj)
                
                if (ndx.eq.1) then
                    dx=(xmax-xmin)/ndx
                else
                    dx=(xmax-xmin)/(ndx-1)
                endif
                
                if (ndy.eq.1) then
                    dy=(ymax-ymin)/ndy
                else
                    dy=(ymax-ymin)/(ndy-1)
                endif
                
                if (ndz.eq.1) then
                    dz=(zmax-zmin)/ndz
                else
                    dz=(zmax-zmin)/(ndz-1)
                endif
            
                do ii1=1,ndx
                    do ii2=1,ndy
                        do ii3=1,ndz
                            loadinp(1,cnt)=xmin+dx*(ii1-1)
                            loadinp(2,cnt)=ymin+dy*(ii2-1)
                            loadinp(3,cnt)=zmin+dz*(ii3-1)
                            loadinp(4,cnt)=tloadx(jj)/tsetpt(jj)
                            loadinp(5,cnt)=tloady(jj)/tsetpt(jj)
                            loadinp(6,cnt)=tloadz(jj)/tsetpt(jj)
                            loadinp(7,cnt)=loadtype(jj)
                            cnt=cnt+1
                        enddo
                    enddo
                enddo
            enddo
        enddo
        
        return
    end subroutine
    
!-----------------------------------------------------------------------------------------------------
!      FUNCTION componentwise_cs_pos_post FOR element connectivity for postpoints (paraview)
!-----------------------------------------------------------------------------------------------------
!	subroutine componentwise_cs_pos_post
!    use var_inputdata
!        implicit none 
!        integer::i,j,ij,cs1,cs2,counter,cs0
!        real*8 ::xcs,zcs
!        real*8, ALLOCATABLE, DIMENSION (:,:,:)  ::cs_pos_post_real
!        
!        allocate(cs_pos_post_real(npts_max,4,be_div-1),cs_pos_post(npts_max,2,be_div-1))
!        cs_pos_post=0
!        
!        !cs0=be_cs_type(1,3)
!        !do i=1,npts(cs0)
!        !    cs_pos_post(i,1,cs0)=i
!        !    cs_pos_post(i,2,cs0)=i
!        !    cs_pos_post(i,3,cs0)=postpoint(i,2,cs0)
!        !    cs_pos_post(i,4,cs0)=postpoint(i,4,cs0)
!        !enddo
!        
!        
!        do i=2,be_div
!            counter=1
!            cs1=be_cs_type(i-1,3)
!            cs2=be_cs_type(i,3)
!            if (npts(cs1).gt.npts(cs2)) then
!                do j=1,npts(cs1)
!                    xcs=postpoint(j,2,cs1)
!                    zcs=postpoint(j,4,cs1)
!                    do ij=1,npts(cs2)
!                        if (xcs.eq.postpoint(ij,2,cs2).and.zcs.eq.postpoint(ij,4,cs2)) then
!                            cs_pos_post_real(ij,1,i-1)=j
!                            cs_pos_post_real(ij,2,i-1)=ij
!                            cs_pos_post_real(ij,3,i-1)=xcs
!                            cs_pos_post_real(ij,4,i-1)=zcs
!                            cs_pos_post(ij,1,i-1)=j
!                            cs_pos_post(ij,2,i-1)=ij
!                            goto 2003
!                        elseif (ij.eq.npts(cs2)) then
!                            cs_pos_post_real(npts(cs2)+counter,1,i-1)=j
!                            cs_pos_post_real(npts(cs2)+counter,2,i-1)=npts(cs2)+counter
!                            cs_pos_post_real(npts(cs2)+counter,3,i-1)=xcs
!                            cs_pos_post_real(npts(cs2)+counter,4,i-1)=zcs                         
!                            cs_pos_post(npts(cs2)+counter,1,i-1)=j
!                            cs_pos_post(npts(cs2)+counter,2,i-1)=npts(cs2)+counter
!                            counter=counter+1
!                        endif
!                    enddo
!2003            enddo
!            elseif (npts(cs1).lt.npts(cs2)) then
!                do j=1,npts(cs2)
!                    xcs=postpoint(j,2,cs2)
!                    zcs=postpoint(j,4,cs2)
!                    do ij=1,npts(cs1)
!                        if (xcs.eq.postpoint(ij,2,cs1).and.zcs.eq.postpoint(ij,4,cs1)) then
!                            cs_pos_post_real(ij,1,i-1)=j
!                            cs_pos_post_real(ij,2,i-1)=ij
!                            cs_pos_post_real(ij,3,i-1)=xcs
!                            cs_pos_post_real(ij,4,i-1)=zcs
!                            cs_pos_post(ij,1,i-1)=j
!                            cs_pos_post(ij,2,i-1)=ij
!                            goto 2004
!                        elseif (ij.eq.npts(cs1)) then
!                            cs_pos_post_real(npts(cs1)+counter,1,i-1)=j
!                            cs_pos_post_real(npts(cs1)+counter,2,i-1)=npts(cs1)+counter
!                            cs_pos_post_real(npts(cs1)+counter,3,i-1)=xcs
!                            cs_pos_post_real(npts(cs1)+counter,4,i-1)=zcs
!                            cs_pos_post(npts(cs1)+counter,1,i-1)=j
!                            cs_pos_post(npts(cs1)+counter,2,i-1)=npts(cs1)+counter
!                            counter=counter+1
!                        endif
!                    enddo
!2004            enddo
!            elseif (npts(cs1).eq.npts(cs2)) then
!                do j=1,npts(cs2)
!                    cs_pos_post_real(j,1,i-1)=j
!                    cs_pos_post_real(j,2,i-1)=j
!                    cs_pos_post_real(j,2,i-1)=xcs
!                    cs_pos_post_real(j,4,i-1)=zcs
!                    cs_pos_post(j,1,i-1)=j
!                    cs_pos_post(j,2,i-1)=j
!                enddo
!            endif  
!        enddo
!             
!        
!        write(11,*)'cs_pos_post cs1'
!        do i=1,npts_max
!            write(11,*)int(cs_pos_post_real(i,1,1)),int(cs_pos_post_real(i,2,1)),cs_pos_post_real(i,3,1),cs_pos_post_real(i,4,1)
!        enddo
!        
!        !write(11,*)'cs_pos_post cs2'
!        !do i=1,npts_max
!        !    write(11,*)int(cs_pos_post_real(i,1,2)),int(cs_pos_post_real(i,2,2)),cs_pos_post_real(i,3,2),cs_pos_post_real(i,4,2)
!        !enddo
!        
!        deallocate(cs_pos_post_real)
!
!      return 
!    end subroutine
    
end module
    