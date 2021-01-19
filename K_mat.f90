module k_mat
contains
!----------------------------------------------------------------------------------------------------------------------------------------------
!     FUNCTION k_mat3D for calculation of K matrix - and initial assembly of sparse matrix (linear and nonlinear analysis and 3D integral)
!----------------------------------------------------------------------------------------------------------------------------------------------
   	subroutine k_mat3D
    use var_inputdata
    use var_analysis
    use struc_3D
    use integ_3D
    use mod_gauss
    use read_input_data
    use omp_lib
    
    implicit none
    integer::i14,j14,bne,agpt,bgpt
    integer::cntt,count,count1,count2,count1_var,count2_var,i,j,t,s,ele,div,cs,i3,j3,k3,elbm,i1,i2,nd_var,ndt,nds,ls,lt,ls1,lt1,ii,jj,ndnm,ndbm,i15,j15,ndnm1,c3,c1
    real*8::alpha,beta,eta,eta_phy,wt,wth,wtk,Term0
    real*8::Gitx,Gity,Gitz,Gjsx,Gjsy,Gjsz,Gitx_temp,Gity_temp,Gitz_temp,Git_temp(3),Term_temp(3),TermCAB(3),TermFDE(3),Git(3),Gjs(3)
    real*8,allocatable,dimension(:)::k_nu
    allocate(k_nu(9))
    
  !  !$omp parallel private (nd_var,count,bne,ndt,nds,lt,ls,count1_var,count2_var,count1,count2,k_nu,eta,eta_phy,wt,alpha,beta,wth,wtk,Term0,TermFDE,TermCAB,Git,Gitx,Gjs,Gjsx,cntt)    
  !  !$omp do
    do elbm=1,ne
        nd_var=(elbm-1)*(nne-1)
        count=3*sum(maxexp_nodes(1:nd_var))
        do ele=1,ane
            if (cs_ele_type(ele,2).eq.3) then
                agpt=agp_tri
                bgpt=1
            elseif (cs_ele_type(ele,2).eq.4) then
                agpt=agp_quad
                bgpt=agp_quad
            endif
            if (elem_beam_coll(ele,2).gt.0.and.elem_beam_coll(ele,2).le.elbm) cycle
            bne=(elbm-1)*ane+ele                                                    ! brick element number
            
            do i=1,nne
                ndt=(elbm-1)*(nne-1)+i
            do j=1,nne
                nds=(elbm-1)*(nne-1)+j
                
            do t=1,mexp
            lt=elemcon_replace(ele,t+1,ndt)
            if (lt.eq.0) cycle
                
            do s=1,mexp
            ls=elemcon_replace(ele,s+1,nds)
            if (ls.eq.0) cycle
                
            count1_var=sum(maxexp_nodes(nd_var+1:nds-1))
            count2_var=sum(maxexp_nodes(nd_var+1:ndt-1))
                        
            count1=(3*count1_var)+(3*ls-2)+count
            count2=(3*count2_var)+(3*lt-2)+count
                    
            if (count2.ge.count1) then
                
                k_nu=0.0d0
                    
                do i3=1,ngp
                eta=gpt_bm(i3)
                eta_phy=((y((nne-1)*elbm+1)-y((nne-1)*(elbm-1)+1))/2.0d0)*(eta+1.0d0)+y((nne-1)*(elbm-1)+1)         ! y-loc in physical csys
                wt=wpt_bm(i3)
                call expfun_beam(eta)
                do j3=1,agpt
                do k3=1,bgpt
                    
                if (cs_ele_type(ele,2).eq.3) then
                    alpha=agpt_cs_t(j3)
                    beta=bgpt_cs_t(j3)
                    wth=wpt_cs_t(j3)
                    wtk=1.0d0  
                elseif (cs_ele_type(ele,2).eq.4) then
                    alpha=gpt_cs(j3)
                    beta=gpt_cs(k3)
                    wth=wpt_cs(j3)
                    wtk=wpt_cs(k3)
                endif

                call brickshapefun(bne,alpha,eta,beta)
                call expfun_cs(ele,alpha,beta)
                call mat_stiff(ele,eta_phy)
                        
                Term0=wt*wth*wtk*det_jac_3D    
                    
                TermCAB(1)=fa(t)*sf(i)
                TermCAB(2)=f(t)*dsf(i)
                TermCAB(3)=fb(t)*sf(i)
                    
                TermFDE(1)=fa(s)*sf(j)
                TermFDE(2)=f(s)*dsf(j)
                TermFDE(3)=fb(s)*sf(j)
                    
                Git=matmul(jac3D_inv,TermCAB)
                Gjs=matmul(jac3D_inv,TermFDE)
                Gitx=Git(1); Gity=Git(2); Gitz=Git(3)
                Gjsx=Gjs(1); Gjsy=Gjs(2); Gjsz=Gjs(3)
                    
                k_nu(1)=k_nu(1)+(trcmat(1,1)*Gitx*Gjsx+trcmat(1,6)*Gity*Gjsx+trcmat(1,5)*Gitz*Gjsx+trcmat(1,6)*Gitx*Gjsy+trcmat(6,6)*Gity*Gjsy &
                                +trcmat(5,6)*Gitz*Gjsy+trcmat(1,5)*Gitx*Gjsz+trcmat(5,6)*Gity*Gjsz+trcmat(5,5)*Gitz*Gjsz)*(Term0)
                    
                k_nu(2)=k_nu(2)+(trcmat(1,6)*Gitx*Gjsx+trcmat(1,2)*Gity*Gjsx+trcmat(1,4)*Gitz*Gjsx+trcmat(6,6)*Gitx*Gjsy+trcmat(2,6)*Gity*Gjsy &
                                +trcmat(4,6)*Gitz*Gjsy+trcmat(5,6)*Gitx*Gjsz+trcmat(2,5)*Gity*Gjsz+trcmat(4,5)*Gitz*Gjsz)*(Term0)
                    
                k_nu(3)=k_nu(3)+(trcmat(1,5)*Gitx*Gjsx+trcmat(1,4)*Gity*Gjsx+trcmat(1,3)*Gitz*Gjsx+trcmat(5,6)*Gitx*Gjsy+trcmat(4,6)*Gity*Gjsy &
                                +trcmat(3,6)*Gitz*Gjsy+trcmat(5,5)*Gitx*Gjsz+trcmat(4,5)*Gity*Gjsz+trcmat(3,5)*Gitz*Gjsz)*(Term0)
                    
                k_nu(4)=k_nu(4)+(trcmat(1,6)*Gitx*Gjsx+trcmat(6,6)*Gity*Gjsx+trcmat(5,6)*Gitz*Gjsx+trcmat(1,2)*Gitx*Gjsy+trcmat(2,6)*Gity*Gjsy &
                                +trcmat(2,5)*Gitz*Gjsy+trcmat(1,4)*Gitx*Gjsz+trcmat(4,6)*Gity*Gjsz+trcmat(4,5)*Gitz*Gjsz)*(Term0)
                    
                k_nu(5)=k_nu(5)+(trcmat(6,6)*Gitx*Gjsx+trcmat(2,6)*Gity*Gjsx+trcmat(4,6)*Gitz*Gjsx+trcmat(2,6)*Gitx*Gjsy+trcmat(2,2)*Gity*Gjsy &
                                +trcmat(2,4)*Gitz*Gjsy+trcmat(4,6)*Gitx*Gjsz+trcmat(2,4)*Gity*Gjsz+trcmat(4,4)*Gitz*Gjsz)*(Term0)
                    
                k_nu(6)=k_nu(6)+(trcmat(5,6)*Gitx*Gjsx+trcmat(4,6)*Gity*Gjsx+trcmat(3,6)*Gitz*Gjsx+trcmat(2,5)*Gitx*Gjsy+trcmat(2,4)*Gity*Gjsy &
                                +trcmat(2,3)*Gitz*Gjsy+trcmat(4,5)*Gitx*Gjsz+trcmat(4,4)*Gity*Gjsz+trcmat(3,4)*Gitz*Gjsz)*(Term0)
                    
                k_nu(7)=k_nu(7)+(trcmat(1,5)*Gitx*Gjsx+trcmat(5,6)*Gity*Gjsx+trcmat(5,5)*Gitz*Gjsx+trcmat(1,4)*Gitx*Gjsy+trcmat(4,6)*Gity*Gjsy &
                                +trcmat(4,5)*Gitz*Gjsy+trcmat(1,3)*Gitx*Gjsz+trcmat(3,6)*Gity*Gjsz+trcmat(3,5)*Gitz*Gjsz)*(Term0)
                    
                k_nu(8)=k_nu(8)+(trcmat(5,6)*Gitx*Gjsx+trcmat(2,5)*Gity*Gjsx+trcmat(4,5)*Gitz*Gjsx+trcmat(4,6)*Gitx*Gjsy+trcmat(2,4)*Gity*Gjsy &
                                +trcmat(4,4)*Gitz*Gjsy+trcmat(3,6)*Gitx*Gjsz+trcmat(2,3)*Gity*Gjsz+trcmat(3,4)*Gitz*Gjsz)*(Term0) 
                    
                k_nu(9)=k_nu(9)+(trcmat(5,5)*Gitx*Gjsx+trcmat(4,5)*Gity*Gjsx+trcmat(3,5)*Gitz*Gjsx+trcmat(4,5)*Gitx*Gjsy+trcmat(4,4)*Gity*Gjsy &
                                +trcmat(3,4)*Gitz*Gjsy+trcmat(3,5)*Gitx*Gjsz+trcmat(3,4)*Gity*Gjsz+trcmat(3,3)*Gitz*Gjsz)*(Term0)
                                       
                enddo
                enddo
                enddo
                    
                cntt=1
                do ii=count1,count1+2
                    do jj=count2,count2+2
                        if (jj.ge.ii) then
                            call ASSEMBLY(ii,jj,k_nu(cntt))
						endif
                        cntt=cntt+1
                    enddo
				enddo 
				
	!!!!!!!!!!!!!!! for getting the condition number of the stiffness matrix !!!!!!!!!!!!!!!!!!!
	!			 cntt=1
    !            do ii=3*maxexp*(j-1)+3*ls-2,3*maxexp*(j-1)+3*ls-2+2
    !                do jj=3*maxexp*(i-1)+3*lt-2,3*maxexp*(i-1)+3*lt-2+2
    !                    if (jj.ge.ii) then
	!						kmat(band_width+ii-jj,count+jj)=kmat(band_width+ii-jj,count+jj)+k_nu(cntt)
    !                    endif
    !                    cntt=cntt+1
    !                enddo
    !            enddo
				
            endif

            enddo
            enddo
            enddo
            enddo
        enddo
        !write(*,*)'K matrix computation and assembly',elbm    
    enddo
  !  !$omp end do
  !  !$omp end parallel
    
    deallocate(k_nu)
    return
    end subroutine 
    
    
!------------------------------------------------------------------------
!   FUNCTION ASSEMBLY FOR initial assembly of sparse matrix
!------------------------------------------------------------------------
    SUBROUTINE ASSEMBLY(IR,IC,VAL)
      use var_inputdata
      use var_analysis

      IMPLICIT NONE
      INTEGER               ::  I,NZ,MAXCOL
      real*8                ::  VAL
      INTEGER               ::  IC,IR
      INTEGER               ::  IA,IB,IAB,A,B,AB,SO

      NZ=MTR_R(IR)%NNZ
      MAXCOL=MTR_R(IR)%COL(NZ)
      
      IF(IC.GT.MAXCOL) THEN
         MTR_R(IR)%VAL(NZ+1)=MTR_R(IR)%VAL(NZ+1)+VAL
         MTR_R(IR)%COL(NZ+1)=IC
         MTR_R(IR)%NNZ=MTR_R(IR)%NNZ+1
      ELSE
        IA=1
        IB=MTR_R(IR)%NNZ
        IAB=INT((IA+IB)/2)
        A=MTR_R(IR)%COL(IA)
        B=MTR_R(IR)%COL(IB)
        AB=MTR_R(IR)%COL(IAB)

10      IF((IB-IA).LT.10) GOTO 20
        
        IF (IC.GE.A.AND.IC.LT.AB) THEN
            IA=IA
            IB=IAB
            IAB=INT((IA+IB)/2)
            A=MTR_R(IR)%COL(IA)
            B=MTR_R(IR)%COL(IB)
            AB=MTR_R(IR)%COL(IAB)
            GOTO 10
        ELSEIF (IC.GE.AB.AND.IC.LE.B) THEN
            IA=IAB
            IB=IB
            IAB=INT((IA+IB)/2)
            A=MTR_R(IR)%COL(IA)
            B=MTR_R(IR)%COL(IB)
            AB=MTR_R(IR)%COL(IAB)
            GOTO 10
        ENDIF

20      DO I=IA,IB
            IF (IC.EQ.MTR_R(IR)%COL(I))THEN
                MTR_R(IR)%VAL(I)=MTR_R(IR)%VAL(I)+VAL
                EXIT
            ELSEIF (IC.LT.MTR_R(IR)%COL(I) )THEN
                MTR_R(IR)%VAL(I+1:NZ+1)=MTR_R(IR)%VAL(I:NZ)
                MTR_R(IR)%COL(I+1:NZ+1)=MTR_R(IR)%COL(I:NZ)
                MTR_R(IR)%VAL(I)=0.0D0
                MTR_R(IR)%VAL(I)=MTR_R(IR)%VAL(I)+VAL
                MTR_R(IR)%COL(I)=IC
                MTR_R(IR)%NNZ=MTR_R(IR)%NNZ+1
                EXIT
            ENDIF
        ENDDO
      ENDIF

      SO=SIZE(MTR_R(IR)%VAL)
      IF(SO-MTR_R(IR)%NNZ<1) THEN
          ALLOCATE (MTR_R(IR)%TEMP(SO))

          MTR_R(IR)%TEMP=0.0d0
          MTR_R(IR)%TEMP=MTR_R(IR)%VAL
          DEALLOCATE (MTR_R(IR)%VAL)
          ALLOCATE  (MTR_R(IR)%VAL(SO*2))
          MTR_R(IR)%VAL=0.0D0
          MTR_R(IR)%VAL(1:SO)=MTR_R(IR)%TEMP(1:SO)

          MTR_R(IR)%TEMP=0.0d0
          MTR_R(IR)%TEMP=MTR_R(IR)%COL
          DEALLOCATE (MTR_R(IR)%COL)
          ALLOCATE  (MTR_R(IR)%COL(SO*2))
          MTR_R(IR)%COL=0
          MTR_R(IR)%COL(1:SO)=INT(MTR_R(IR)%TEMP(1:SO))
          DEALLOCATE (MTR_R(IR)%TEMP)
      ENDIF
      
     return
    end subroutine
    
!----------------------------------------------------------------
!   FUNCTION kmat_sparse FOR FORMATION OF global sparse MATRIX
!----------------------------------------------------------------
    SUBROUTINE  kmat_sparse
      use var_inputdata
      use var_analysis
      
      IMPLICIT NONE
      INTEGER                :: I,J,JJ
      REAL*8                 :: MFACTOR

      tot_nnz=0

      DO I=1,tdof
        tot_nnz=tot_nnz+MTR_R(I)%NNZ
      ENDDO

      ALLOCATE(k_val(tot_nnz))
      ALLOCATE(col_idx(tot_nnz))
      ALLOCATE(row_ptr(tdof+1))
      
      col_idx=0
      k_val=0.0D0
      row_ptr=0
      diag_idx=0
      !K_MTR%NR=tdof
      !K_MTR%NC=tdof

      JJ=1
      row_ptr(1)=1
      DO I=1,tdof
          row_ptr(I+1)=row_ptr(I)+MTR_R(I)%NNZ
          DO J=1,MTR_R(I)%NNZ
              !if (J.eq.1) then
              !    diag_idx(I)=JJ
              !endif
            k_val(JJ)=MTR_R(I)%VAL(J)
            col_idx(JJ)=MTR_R(I)%COL(J)
            JJ=JJ+1
          ENDDO
      ENDDO

      MFACTOR=tot_nnz/(tdof**2.0D0)*100.0D0
      
      
      DO i=1,tot_nnz
        write(13,*)col_idx(i),k_val(i)
      ENDDO
      
      write(*,*)tot_nnz,tdof,MFACTOR
      return
    end subroutine
    
end module