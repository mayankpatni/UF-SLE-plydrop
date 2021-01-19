module mod_math
contains

!---------------------------------------------------------------------------------------
!      FUNCTION phy_to_norCS_1 FOR transforming physical coordinate to normal coordinate
!---------------------------------------------------------------------------------------
	subroutine phy_to_norCS_1
    use var_inputdata
    use var_analysis
    
        implicit none
        integer::normcnt,phycnt,ii1,ii2,ii3,ii4,nod,i
        real*8::norm1,norm2,dnorm,xphynode,zphynode,alpha1,beta1,ft1(4),norm1_beta,norm2_beta,dnorm_beta
        
        tpts(1)=(phynormesh+1)*(phynormesh+1)                       ! quad
        tpts(2)=(phynormesh+1)*(phynormesh+2)/2                     ! tri
        
        allocate(normcs(tpts(1),2),normcs_tri(tpts(2),2),phycs(count(cs_ele_type(:,2).eq.3)*tpts(2)+count(cs_ele_type(:,2).eq.4)*tpts(1),3))
        phycs=0.0
        
        norm1=-1.0d0
        norm2=1.0d0
        dnorm=(norm2-norm1)/phynormesh
        
        norm1_beta=0.0d0
        norm2_beta=sqrt(3.0d0)
        dnorm_beta=(norm2_beta-norm1_beta)/phynormesh

        !!!!!!! quad !!!!!!!!
        normcnt=0
        do ii1=1,phynormesh+1
            do ii2=1,phynormesh+1
                normcnt=normcnt+1
                normcs(normcnt,1)=norm1+dnorm*(ii2-1)
                normcs(normcnt,2)=norm1+dnorm*(ii1-1)
            enddo
        enddo
        
        !!!!!!!! tri !!!!!!!!!
        i=0
        normcnt=0
        do ii1=1,phynormesh+1
            do ii2=1,phynormesh+1-i
                normcnt=normcnt+1
                normcs_tri(normcnt,1)=norm1+dnorm*(ii2-1)+dnorm*(ii1-1)/2
                normcs_tri(normcnt,2)=norm1_beta+dnorm_beta*(ii1-1)
            enddo
            i=i+1
        enddo

        phycnt=0
        do ii3=1,ane
            if (cs_ele_type(ii3,2).eq.3) then                ! tri
                do ii1=1,tpts(2)
                    alpha1=normcs_tri(ii1,1)
                    beta1=normcs_tri(ii1,2)
                    ft1(1)=0.5d0*(1.0d0-alpha1-beta1/sqrt(3.0d0))   
                    ft1(2)=0.5d0*(1.0d0+alpha1-beta1/sqrt(3.0d0)) 
                    ft1(3)=beta1/sqrt(3.0d0) 
                    phycnt=phycnt+1
                    phycs(phycnt,1)=ii3
                
                    do ii2=1,3
                        nod=elemcon(ii3,ii2+1)
                        xphynode=cscord(nod,2)
                        zphynode=cscord(nod,3)
                        phycs(phycnt,2)=phycs(phycnt,2)+xphynode*ft1(ii2)
                        phycs(phycnt,3)=phycs(phycnt,3)+zphynode*ft1(ii2)
                    enddo
                enddo
            elseif (cs_ele_type(ii3,2).eq.4) then            ! quad
                do ii1=1,tpts(1)
                    alpha1=normcs(ii1,1)
                    beta1=normcs(ii1,2)
                    ft1(1)=0.25d0*(1.0d0-alpha1)*(1.0d0-beta1)   
                    ft1(2)=0.25d0*(1.0d0+alpha1)*(1.0d0-beta1)
                    ft1(3)=0.25d0*(1.0d0+alpha1)*(1.0d0+beta1)
                    ft1(4)=0.25d0*(1.0d0-alpha1)*(1.0d0+beta1) 
                    phycnt=phycnt+1
                    phycs(phycnt,1)=ii3
                
                    do ii2=1,4
                        nod=elemcon(ii3,ii2+1)
                        xphynode=cscord(nod,2)
                        zphynode=cscord(nod,3)
                        phycs(phycnt,2)=phycs(phycnt,2)+xphynode*ft1(ii2)
                        phycs(phycnt,3)=phycs(phycnt,3)+zphynode*ft1(ii2)
                    enddo
                enddo
            endif
		enddo
        
		!do i=1,tpts(2)
		!	write(11,*)normcs_tri(i,1),normcs_tri(i,2)
		!	write(13,*)phycs(i,2),phycs(i,3)
		!enddo
		
      return 
    end subroutine

!------------------------------------------------------------------------------------------------------------
!      FUNCTION phy_to_norCS_2_load FOR transforming physical coordinate to normal coordinate- for load vector
!-----------------------------------------------------------------------------------------------------------
	subroutine phy_to_norCS_2_load(xp,zp,alpha,beta,ele_cs)
    use var_inputdata
        implicit none
        integer::phycnt,ii1,minlc(1),ele_cs,minloc1,minloc2,counter
        real*8::alpha,beta,mini,xp,zp
        real*8,allocatable,dimension(:)::phycs1
        allocate(phycs1(count(cs_ele_type(:,2).eq.3)*tpts(2)+count(cs_ele_type(:,2).eq.4)*tpts(1)))
        
        phycnt=count(cs_ele_type(:,2).eq.3)*tpts(2)+count(cs_ele_type(:,2).eq.4)*tpts(1)
        phycs1=0.0
        
        do ii1=1,phycnt
            phycs1(ii1)=((phycs(ii1,2)-xp)**2)+((phycs(ii1,3)-zp)**2)
        enddo
        
        mini=minval(phycs1)
        minlc=minloc(phycs1)
        minloc1=minlc(1)
        
        if (mini.le.0.001) then
            ele_cs=phycs(minloc1,1)
            counter=count(cs_ele_type(1:ele_cs-1,2).eq.3)*tpts(2)+count(cs_ele_type(1:ele_cs-1,2).eq.4)*tpts(1)
            minloc2=minloc1-counter
            if (cs_ele_type(ele_cs,2).eq.3)then
                alpha=normcs_tri(minloc2,1)
                beta=normcs_tri(minloc2,2)
            elseif (cs_ele_type(ele_cs,2).eq.4)then
                alpha=normcs(minloc2,1)
                beta=normcs(minloc2,2)
            endif
        else
            ele_cs=0
            alpha=0
            beta=0
        endif
        
        deallocate(phycs1)
      return 
    end subroutine

    
!-----------------------------------------------------------------------------------------------------------
!      FUNCTION phy_to_norCS_2_post FOR transforming physical coordinate to normal coordinate - post processing
!-----------------------------------------------------------------------------------------------------------
	subroutine phy_to_norCS_2_post(xp,zp,alpha_p,beta_p,ele_cs_p,alpbet_cnt)
    use var_inputdata
        implicit none
        integer::phycnt,ii1,minlocat,ii2,alpbet_cnt,counter
        real*8::mini,xp,zp
        real*8,allocatable,dimension(:)::phycs1,alpha_p,beta_p
        integer,allocatable,dimension(:)::ele_cs_p
        allocate(phycs1(count(cs_ele_type(:,2).eq.3)*tpts(2)+count(cs_ele_type(:,2).eq.4)*tpts(1)))
        
        phycnt=count(cs_ele_type(:,2).eq.3)*tpts(2)+count(cs_ele_type(:,2).eq.4)*tpts(1)
        phycs1=0.0
        ele_cs_p=0
        alpha_p=0.0
        beta_p=0.0
        
        do ii1=1,phycnt
            phycs1(ii1)=((phycs(ii1,2)-xp)**2)+((phycs(ii1,3)-zp)**2)
        enddo
        
        mini=minval(phycs1)
        
        alpbet_cnt=0
        do ii2=1,phycnt
            if (phycs1(ii2).eq.mini) then
                alpbet_cnt=alpbet_cnt+1
                ele_cs_p(alpbet_cnt)=phycs(ii2,1)
                counter=count(cs_ele_type(1:ele_cs_p(alpbet_cnt)-1,2).eq.3)*tpts(2)+count(cs_ele_type(1:ele_cs_p(alpbet_cnt)-1,2).eq.4)*tpts(1)
                minlocat=ii2-counter
                if (cs_ele_type(ele_cs_p(alpbet_cnt),2).eq.3)then
                    alpha_p(alpbet_cnt)=normcs_tri(minlocat,1)
                    beta_p(alpbet_cnt)=normcs_tri(minlocat,2)
                elseif (cs_ele_type(ele_cs_p(alpbet_cnt),2).eq.4)then
                    alpha_p(alpbet_cnt)=normcs(minlocat,1)
                    beta_p(alpbet_cnt)=normcs(minlocat,2)
                endif
            endif
        enddo
        
        deallocate(phycs1)
      return 
    end subroutine
    
    
!-----------------------------------------------------------------------------------------------------------------------------
!      FUNCTION phy_to_norCS_2 FOR transforming z physical coordinate to beta normal coordinate - post processing- recovery
!-----------------------------------------------------------------------------------------------------------------------------
	!subroutine z_to_beta_recovery(xm1,zp,alpha,beta,ele_cs)
 !   use var_inputdata
 !       implicit none
 !       integer::phycnt,ii1,minlocat,ii2,ele_cs
 !       real*8::mini,xm1,zp,alpha,beta
 !       real*8,allocatable,dimension(:)::phycs1
 !       allocate(phycs1(tpts*ane))
 !       
 !       phycnt=tpts*ane
 !       phycs1=0.0
 !
 !       do ii1=1,phycnt
 !           phycs1(ii1)=((phycs(ii1,2)-xm1)**2)+((phycs(ii1,3)-zp)**2)
 !       enddo
 !       
 !       mini=minval(phycs1)
 !       
 !       do ii2=1,phycnt
 !           if (phycs1(ii2).eq.mini.and.phycs(ii2,1).eq.ele_cs) then
 !               minlocat=ii2-(ele_cs-1)*tpts
 !               !alpha=normcs(minlocat,1)
 !               beta=normcs(minlocat,2)
 !           endif
 !       enddo
 !       
 !       deallocate(phycs1)
 !     return 
 !   end subroutine

!-----------------------------------------------------------------------------------------
!      FUNCTION phy_to_norBeam FOR transforming physical coordinate to normal coordinate
!-----------------------------------------------------------------------------------------
	subroutine phy_to_norBeam(yp,eta,ele_beam)
    use var_inputdata
        implicit none
        integer::ii1,ele_beam,ynod1,ynod2
        real*8::eta,yp1,yp2,yp
        logical::FLAG
        ele_beam=0
        
        FLAG=.TRUE.
        do ii1=1,ne
            ynod1=((nne-1)*(ii1-1))+1
            ynod2=((nne-1)*ii1)+1
        
            if (yp.le.y(ynod2).and.FLAG .eq. .TRUE.) then
                ele_beam=ii1
                yp1=y(ynod1)
                yp2=y(ynod2)
                eta=-1.0d0+(2.0d0*(yp-yp1)/(yp2-yp1))
                FLAG=.FALSE.
            endif
        enddo
        
      return 
    end subroutine
 
!--------------------------------------
!       matrix print
!--------------------------------------
		SUBROUTINE PRINTF(A,N,M) 
		IMPLICIT NONE
		INTEGER:: N,I,J, M
		real*8, DIMENSION (N,M):: A
       
	   DO I = 1, N
	  
        WRITE(11,1003)  A(I,:)
	1003 format(//,24(en18.7,1x))
       END DO
	   
       WRITE(11,*)       
        END SUBROUTINE 
        
!--------------------------------------
!       integer matrix print
!--------------------------------------
		SUBROUTINE PRINTI(A,N,M) 
		IMPLICIT NONE
		INTEGER:: N,I,J, M
		integer, DIMENSION (N,M):: A
       
	   DO I = 1, N
	  
        WRITE(11,1001)  A(I,:)
	1001 format(//,24(I4,1x))
       END DO
	   
       WRITE(11,*)       
        END SUBROUTINE 


! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------
   INTEGER FUNCTION  FindMinimum(x, Start, endd)
      IMPLICIT  NONE
      real*8, DIMENSION(1:), INTENT(IN)  :: x
      INTEGER, INTENT(IN)                :: Start, endd
      real*8                             :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		                ! assume the first is the min
      Location = Start			                ! record its position
      DO i = Start+1, endd		                ! start with next elements
         IF (x(i) < Minimum) THEN	            !   if x(i) less than the min?
            Minimum  = x(i)		                !      Yes, a new minimum found
            Location = i                        !      record its position
         END IF
      END DO
      FindMinimum = Location        	        ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real*8, INTENT(INOUT)     :: a, b
      real*8                    :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      real*8, DIMENSION(1:), INTENT(INOUT)  :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1			                            ! except for the last
         Location = FindMinimum(x, i, Size)	                ! find min from this to last
         CALL  Swap(x(i), x(Location))	                    ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
   
   
!-----------------------------------------------------------------------------------------
!      FUNCTION rot_mat FOR rotation matrix
!-----------------------------------------------------------------------------------------
	subroutine rot_mat(alp,bet,gam,rmat)
    use var_inputdata
        implicit none
        real*8::alp,bet,gam
        real*8,allocatable,dimension(:,:)::rmat,rmat_x,rmat_y,rmat_z
        
        allocate(rmat_x(3,3),rmat_y(3,3),rmat_z(3,3))
        rmat_x=0.0
        rmat_y=0.0
        rmat_z=0.0
        
        rmat_x(1,1)=1.0
        rmat_x(2,2)=cosd(alp)
        rmat_x(2,3)=-sind(alp)
        rmat_x(3,2)=sind(alp)
        rmat_x(3,3)=cosd(alp)
        
        rmat_y(2,2)=1.0
        rmat_y(1,1)=cosd(bet)
        rmat_y(3,1)=-sind(bet)
        rmat_y(1,3)=sind(bet)
        rmat_y(3,3)=cosd(bet)
        
        rmat_z(3,3)=1.0
        rmat_z(1,1)=cosd(gam)
        rmat_z(1,2)=-sind(gam)
        rmat_z(2,1)=sind(gam)
        rmat_z(2,2)=cosd(gam)
        
        rmat=matmul(rmat_x,matmul(rmat_y,rmat_z))
        
        deallocate(rmat_x,rmat_y,rmat_z)
      return 
    end subroutine
    
!-----------------------------------------------------------------------------------------
!      FUNCTION inverse FOR matrix inversion
!-----------------------------------------------------------------------------------------  
  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real*8::a(n,n), c(n,n)
real*8::L(n,n), U(n,n), b(n), d(n), x(n)
real*8::coeff
integer::i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
  end subroutine inverse
  
!-----------------------------------------------------------------------------------------
!      FUNCTION matinv3 FOR matrix inversion of 3x3 matrix
!-----------------------------------------------------------------------------------------    
    subroutine matinv3(A,B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real*8 :: A(3,3), B(3,3)  !! Matrix
    real*8 :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
	end subroutine matinv3
        
!-----------------------------------------------------------------------------------------
!      FUNCTION cond_number FOR evaluating condition number of Kmatrix 
!-----------------------------------------------------------------------------------------
	subroutine cond_number
    use var_inputdata
    use var_analysis
        implicit none
          CHARACTER::UPLO,norm
          INTEGER ::INFO, KD, LDAB, N
          real*8::RCOND,ANORM
          INTEGER,allocatable::IWORK(:)
          real*8,allocatable::AB(:,:), WORK(:)
          
          real*8:: dlansb
          EXTERNAL dlansb

          EXTERNAL DPBTRF, DPBCON
          
      KD=band_width-1
      LDAB=KD+1
      UPLO='U'
      norm='1'
      N=tdof1
      
      allocate(AB(LDAB,N),IWORK(N),WORK(3*N))
      AB(:,:)=kmat(:,3*maxexp+1:tdof)
      
         ANORM=dlansb(norm,UPLO,N,KD,AB,LDAB,WORK)
       call DPBTRF(UPLO,N,KD,AB,LDAB,INFO)
       call DPBCON(UPLO,N,KD,AB,LDAB,ANORM,RCOND,WORK,IWORK,INFO)
       
       write(11,*)'condition number',rcond
       
       deallocate(AB,IWORK,WORK)
        
      return 
	end subroutine
	
end module
