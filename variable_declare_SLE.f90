module var_inputdata
    
    implicit none
    integer                                 ::analysis_type, failure_model, fullorpost
    
    integer                                 ::nlayer,nmat
    integer, ALLOCATABLE, DIMENSION (:)     ::matin
    real*8                                  ::glb_e(3,3),loc_e(3,3),loc_struc(3,3),loc_e_inv(3,3),loc_struc_inv(3,3)
    real*8, ALLOCATABLE, DIMENSION (:,:)    ::angle,trcmat,mat_strength
    real*8, ALLOCATABLE, DIMENSION (:,:,:)  ::cmat
    
    integer                                 ::ne,nne,ngp,tnodes,tdof,tdof1,band_width
    real*8                                  ::ylengthi,ylengthf
    real*8, ALLOCATABLE, DIMENSION (:)      ::y
    CHARACTER*4                             ::beam_mesh
    
    integer                                 ::nexp,mexp,maxexp,agp_quad,agp_tri
    integer                                 ::anodes3D,maxexp_cs,ane,anode,bnode3D,bne3D,numcs_collap
    integer, ALLOCATABLE, DIMENSION (:)     ::maxexp_nodes,belemcon3D_type
    integer, ALLOCATABLE, DIMENSION (:,:)   ::elemcon,layerno,cs_ele_type,elemcon3D,belemcon3D,elem_beam_coll
    integer,allocatable,dimension (:,:,:)   ::elem_collap,elemcon_replace
    real*8, ALLOCATABLE, DIMENSION (:,:)    ::cscord,cscord3D,bcord3D
    
    integer                                 ::phynormesh
    real*8, ALLOCATABLE, DIMENSION (:,:)    ::normcs,normcs_tri
    real*8, ALLOCATABLE, DIMENSION (:,:)    ::phycs
    integer                                 ::tpts(2)
    
    integer                                 ::np,num_bcinp
    real*8, ALLOCATABLE, DIMENSION (:,:)    ::loadinp
    integer, ALLOCATABLE, DIMENSION (:,:)   ::node_bc,cs_node_bc,ane_load
    integer, ALLOCATABLE, DIMENSION (:)     ::ane_load_num
    
    integer                                 ::npts_max,be_div_post,be_div_post_p,npts,post_elem
    integer,ALLOCATABLE, DIMENSION (:,:)    ::tnod_npts
    real*8, ALLOCATABLE, DIMENSION (:,:)    ::postpoint
    integer, ALLOCATABLE, DIMENSION (:,:)   ::post_elecon
    
        
    !type my_type
    !    integer :: col1,col2
    !    real*8:: col3,col4
    !endtype
    !
    !type(mytype), dimension(:,:,:), allocatable:: cs_pos_post

end module


module var_analysis

implicit none
    integer                                         ::totnodes,tot_nnz
    integer,allocatable, dimension (:)              ::col_idx,row_ptr
    real*8, ALLOCATABLE, DIMENSION (:)              ::gpt_cs,wpt_cs,agpt_cs_t,bgpt_cs_t,wpt_cs_t,gpt_bm,wpt_bm,agpt_cs_t_test,bgpt_cs_t_test,wpt_cs_t_test
    real*8, ALLOCATABLE, DIMENSION (:)              ::jacob,jacobinv,sf,dsf,ft,fs,ftx,fsx,ftz,fsz,ftxx,ftxz,ftzz,ftxxx,ftxxz,dddsf,f,fa,fb,fq,fqx,fqz,fp,fpx,fpz
    real*8, allocatable, dimension(:,:)             ::ft_test,sf_test
    real*8, ALLOCATABLE, DIMENSION (:,:,:,:,:)      ::ninj,ninjd,nidnj,nidnjd
    real*8, ALLOCATABLE, DIMENSION (:,:)            ::ftzfs,ftfsz,ftzfsz,ftxfsz,ftzfsx,ftfs,ftxfs,ftfsx,ftxfsx,kmat
    
    real*8, allocatable, dimension (:)              ::ninjnlnm_test,ninjnlnm_testS,ninj_test,ninjnl_test,ninjnl_test2,ninjnl_testS
    real*8, ALLOCATABLE, DIMENSION (:)              ::ftfsfpfq_test,ftfs_test,ftfsfp_test

    
    real*8, ALLOCATABLE, DIMENSION (:,:,:)          ::ktsij
    real*8                                          ::jac3D(3,3),jac3D_inv(3,3),det_jac_3D
    real                                            ::app_load
    integer,allocatable, dimension (:)              ::diag_idx
    integer                                         ::cpos(33)
    
    
    real*8, ALLOCATABLE, DIMENSION (:)              ::uvect,k_val,del_uvect,pvect_res,fvect,uconverged
    real*8,ALLOCATABLE, DIMENSION (:,:)             ::strn_G,strs_G,strs_Le,strs_Ls,disp,failure_val,strs_C,strs_Cd
    
    TYPE MTRR
     real*8, DIMENSION(:), ALLOCATABLE              ::VAL
     real*8, DIMENSION(:), ALLOCATABLE              ::TEMP
     INTEGER, DIMENSION(:), ALLOCATABLE             ::COL,col_real
     INTEGER                                        ::NNZ
    END TYPE MTRR

    TYPE(MTRR),DIMENSION(:), ALLOCATABLE            ::MTR_R

end module

