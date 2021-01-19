module integ_2D_cs
contains

!--------------------------------------------------------------------------------------------------
!     FUNCTION cuf_2D_int for expansion function and their derivatives and surface integral for SLE
!--------------------------------------------------------------------------------------------------
   	subroutine cuf_2D_int(ele,cs)
    use var_inputdata
    use var_analysis
    use mod_gauss
    
        implicit none
        integer::i5,j5,i3,j3,ele,cs,cc,dd,hh
        real*8::alpha,beta,wth,wtk,sljacdet
        ALLOCATE(ft(mexp),fs(mexp),ftx(mexp),fsx(mexp),ftz(mexp),fsz(mexp),ftxx(mexp),ftxz(mexp),ftzz(mexp),ftxxx(mexp),ftxxz(mexp))
        allocate(ft_test(mexp,3))
        
        !ftfs=0.0d0
        !ftxfs=0.0d0
        !ftfsx=0.0d0
        !ftxfsx=0.0d0
        !ftzfs=0.0d0
        !ftfsz=0.0d0
        !ftzfsz=0.0d0
        !ftxfsz=0.0d0
        !ftzfsx=0.0d0
        ftfs_test=0.0d0
        
        allocate(gpt(agp),wpt(agp))   
        call gauss(agp)
             
	  	do i3=1,agp
            do j3=1,agp
    	        alpha=gpt(i3)
                beta=gpt(j3)
                wth=wpt(i3)
                wtk=wpt(j3)
     	                
                call slexpfun(ele,alpha,beta,sljacdet,cs)
                ft_test(:,1)=ftx(:)
                ft_test(:,2)=ft(:)
                ft_test(:,3)=ftz(:)
                            
                do i5=1,mexp
                    do j5=1,mexp
                        do cc=1,3
                            do dd=1,3 
                                hh=9*(mexp*(i5-1)+(j5-1))+3*(cc-1)+dd                                   
                                ftfs_test(hh)=ftfs_test(hh)+(sljacdet*wth*wtk*ft_test(i5,cc)*ft_test(j5,dd))
                            enddo
                        enddo
                        !ftfs(i5,j5)=ftfs(i5,j5)+(sljacdet*wth*wtk*ft(i5)*fs(j5))
                        !ftxfs(i5,j5)=ftxfs(i5,j5)+(sljacdet*wth*wtk*ftx(i5)*fs(j5))                            
                        !ftfsx(i5,j5)=ftfsx(i5,j5)+(sljacdet*wth*wtk*ft(i5)*fsx(j5))
                        !ftxfsx(i5,j5)=ftxfsx(i5,j5)+(sljacdet*wth*wtk*ftx(i5)*fsx(j5))
                        !ftzfs(i5,j5)=ftzfs(i5,j5)+(sljacdet*wth*wtk*ftz(i5)*fs(j5))
                        !ftfsz(i5,j5)=ftfsz(i5,j5)+(sljacdet*wth*wtk*ft(i5)*fsz(j5))
                        !ftzfsz(i5,j5)=ftzfsz(i5,j5)+(sljacdet*wth*wtk*ftz(i5)*fsz(j5))
                        !ftxfsz(i5,j5)=ftxfsz(i5,j5)+(sljacdet*wth*wtk*ftx(i5)*fsz(j5))
                        !ftzfsx(i5,j5)=ftzfsx(i5,j5)+(sljacdet*wth*wtk*ftz(i5)*fsx(j5))         
                    enddo
                enddo
            enddo
        enddo
            
    DEALLOCATE (ft,fs,ftx,fsx,ftz,fsz,ftxx,ftxz,ftzz,ftxxx,ftxxz,ft_test)
    deallocate(gpt,wpt)
    return
    end subroutine

!--------------------------------------------------------------------------------------------------
!     FUNCTION cuf_2D_int_nl for expansion function and their derivatives and surface integral for SLE
!--------------------------------------------------------------------------------------------------
   	subroutine cuf_2D_int_nl(ele,cs)
    use var_inputdata
    use var_analysis
    use mod_gauss
    
        implicit none
        integer::i5,j5,k5,l5,i3,j3,ele,cs,aa,bb,cc,dd,hh
        real*8::alpha,beta,wth,wtk,sljacdet
        ALLOCATE(ft(mexp),fs(mexp),ftx(mexp),fsx(mexp),ftz(mexp),fsz(mexp),ftxx(mexp),ftxz(mexp),ftzz(mexp),ftxxx(mexp),ftxxz(mexp))
        allocate(ft_test(mexp,3))
  
        ftfsfpfq_test=0.0d0; ftfsfp_test=0.0d0; ftfs_test=0.0d0
        
        allocate(gpt(agp),wpt(agp))   
        call gauss(agp)
             
	  	do i3=1,agp
            do j3=1,agp
    	        alpha=gpt(i3)
                beta=gpt(j3)
                wth=wpt(i3)
                wtk=wpt(j3)
     	                
                call slexpfun(ele,alpha,beta,sljacdet,cs)
                ft_test(:,1)=ftx(:)
                ft_test(:,2)=ft(:)
                ft_test(:,3)=ftz(:)
            
                do i5=1,mexp
                    do j5=1,mexp
                        do cc=1,3
                            do dd=1,3 
                                hh=9*(mexp*(i5-1)+(j5-1))+3*(cc-1)+dd                                   
                                ftfs_test(hh)=ftfs_test(hh)+(sljacdet*wth*wtk*ft_test(i5,cc)*ft_test(j5,dd))
                            enddo
                        enddo
                    enddo
                enddo
                
                do i5=1,mexp
                    do j5=1,mexp
                        do k5=1,mexp
                            do bb=1,3
                                do cc=1,3
                                    do dd=1,3 
                                        hh=27*(mexp*mexp*(i5-1)+mexp*(j5-1)+(k5-1))+9*(bb-1)+3*(cc-1)+dd                             
                                        ftfsfp_test(hh)=ftfsfp_test(hh)+(sljacdet*wth*wtk*ft_test(i5,bb)*ft_test(j5,cc)*ft_test(k5,dd))
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
                
                do i5=1,mexp
                    do j5=1,mexp
                        do k5=1,mexp
                           do l5=1,mexp
                               do aa=1,3
                                   do bb=1,3
                                       do cc=1,3
                                           do dd=1,3
                                               hh=81*(mexp*mexp*mexp*(i5-1)+mexp*mexp*(j5-1)+mexp*(k5-1)+l5-1)+27*(aa-1)+9*(bb-1)+3*(cc-1)+dd                                   
                                               ftfsfpfq_test(hh)=ftfsfpfq_test(hh)+(sljacdet*wth*wtk*ft_test(i5,aa)*ft_test(j5,bb)*ft_test(k5,cc)*ft_test(l5,dd))
                                           enddo
                                       enddo
                                   enddo
                               enddo
                           enddo
                        enddo
                    enddo
                enddo
                
            enddo
        enddo
            
    DEALLOCATE (ft,fs,ftx,fsx,ftz,fsz,ftxx,ftxz,ftzz,ftxxx,ftxxz,ft_test)
    deallocate(gpt,wpt)
    return
    end subroutine

!-----------------------------------------------------------------
!     FUNCTION slexpfun for evaluating expansion function for SLE 
!-----------------------------------------------------------------
   	subroutine slexpfun(ele,alpha,beta,sljacdet,cs)
    use var_inputdata
    use var_analysis
        implicit none
        integer::i1,ele,nodenum,cs
        real*8::alpha,beta
        real*8::sljacdet,xa,xb,za,zb,xab,zab,afun,bfun,cfun,dfun,efun,ffun,pfun,qfun,rfun,sfun,aafun,bbfun,ccfun,ddfun
        real*8,ALLOCATABLE, DIMENSION (:)::fta,ftb,ftaa,ftbb,ftab,ftaaa,ftaab,ftaba,ftabb,ftbba,ftbbb
        ALLOCATE(fta(mexp),ftb(mexp),ftaa(mexp),ftbb(mexp),ftab(mexp),ftaaa(mexp),ftaab(mexp),ftaba(mexp),ftabb(mexp),ftbba(mexp),ftbbb(mexp))
        xa=0.0d0
        xb=0.0d0
        za=0.0d0
        zb=0.0d0
        xab=0.0d0
        zab=0.0d0
        
          ft(1)=0.25d0*(1.0d0-alpha)*(1.0d0-beta)   
          ft(2)=0.25d0*(1.0d0+alpha)*(1.0d0-beta)
          ft(3)=0.25d0*(1.0d0+alpha)*(1.0d0+beta)
          ft(4)=0.25d0*(1.0d0-alpha)*(1.0d0+beta) 
          
          fta(1)=-0.25d0*(1.0d0-beta)   
          fta(2)=0.25d0*(1.0d0-beta)
          fta(3)=0.25d0*(1.0d0+beta)
          fta(4)=-0.25d0*(1.0d0+beta)
          
          ftb(1)=-0.25d0*(1.0d0-alpha)  
          ftb(2)=-0.25d0*(1.0d0+alpha)
          ftb(3)=0.25d0*(1.0d0+alpha)
          ftb(4)=0.25d0*(1.0d0-alpha)
          
          ftaa(1)=0.0; ftaa(2)=0.0; ftaa(3)=0.0; ftaa(4)=0.0    
          ftbb(1)=0.0; ftbb(2)=0.0; ftbb(3)=0.0; ftbb(4)=0.0 
          ftab(1)=0.25d0; ftab(2)=-0.25d0; ftab(3)=0.25d0; ftab(4)=-0.25d0
          
          ftaaa(1)=0.0; ftaaa(2)=0.0; ftaaa(3)=0.0; ftaaa(4)=0.0
          ftaab(1)=0.0; ftaab(2)=0.0; ftaab(3)=0.0; ftaab(4)=0.0
          
          ftbba(1)=0.0; ftbba(2)=0.0; ftbba(3)=0.0; ftbba(4)=0.0
          ftbbb(1)=0.0; ftbbb(2)=0.0; ftbbb(3)=0.0; ftbbb(4)=0.0
          
          ftaba(1)=0.0; ftaba(2)=0.0; ftaba(3)=0.0; ftaba(4)=0.0
          ftabb(1)=0.0; ftabb(2)=0.0; ftabb(3)=0.0; ftabb(4)=0.0
          
          if (mexp.eq.4) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
          
          ft(5)=0.5d0*(1.0d0-beta)*((alpha**2)-1.0d0)
          ft(6)=0.5d0*(1.0d0+alpha)*((beta**2)-1.0d0)
          ft(7)=0.5d0*(1.0d0+beta)*((alpha**2)-1.0d0)
          ft(8)=0.5d0*(1.0d0-alpha)*((beta**2)-1.0d0)
          
          fta(5)=0.5d0*(1.0d0-beta)*(2.0*alpha)
          fta(6)=0.5d0*((beta**2)-1.0d0)
          fta(7)=0.5d0*(1.0d0+beta)*(2.0*alpha)
          fta(8)=-0.5d0*((beta**2)-1.0d0)
          
          ftb(5)=-0.5d0*((alpha**2)-1.0d0)
          ftb(6)=0.5d0*(1.0d0+alpha)*(2.0*beta)
          ftb(7)=0.5d0*((alpha**2)-1.0d0)
          ftb(8)=0.5d0*(1.0d0-alpha)*(2.0d0*beta)
          
          ftaa(5)=0.5d0*(1.0d0-beta)*(2.0); ftaa(6)=0.0d0; ftaa(7)=0.5d0*(1.0d0+beta)*(2.0); ftaa(8)=0.0d0
          ftbb(5)=0.0d0; ftbb(6)=0.5d0*(1.0d0+alpha)*(2.0); ftbb(7)=0.0d0; ftbb(8)=0.5d0*(1.0d0-alpha)*(2.0)
          ftab(5)=-0.5d0*(2.0*alpha); ftab(6)=0.5d0*(2.0*beta); ftab(7)=0.5d0*(2.0*alpha); ftab(8)=-0.5d0*(2.0*beta)
          
          ftaaa(5)=0.0d0; ftaaa(6)=0.0d0; ftaaa(7)=0.0d0; ftaaa(8)=0.0d0
          ftbba(5)=0.0d0; ftbba(6)=1.0d0; ftbba(7)=0.0d0; ftbba(8)=-1.0d0
          ftaba(5)=-1.0d0; ftaba(6)=0.0d0; ftaba(7)=1.0d0; ftaba(8)=0.0d0
          
          ftaab(5)=-1.0d0; ftaab(6)=0.0d0; ftaab(7)=1.0d0; ftaab(8)=0.0d0
          ftbbb(5)=0.0d0; ftbbb(6)=0.0d0; ftbbb(7)=0.0d0; ftbbb(8)=0.0d0
          ftabb(5)=0.0d0; ftabb(6)=1.0d0; ftabb(7)=0.0d0; ftabb(8)=-1.0d0
          
          if (mexp.eq.8) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(9)=0.5d0*(1.0d0-beta)*((alpha**3)-alpha)
          ft(10)=0.5d0*(1.0d0+alpha)*((beta**3)-beta)
          !ft(11)=0.5d0*(1.0d0+beta)*((alpha**3)-alpha)
          !ft(12)=0.5d0*(1.0d0-alpha)*((beta**3)-beta)
          ft(11)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**3)+alpha)
          ft(12)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**3)+beta)
          
          fta(9)=0.5d0*(1.0d0-beta)*(3.0d0*(alpha**2)-1.0d0)
          fta(10)=0.5d0*((beta**3)-beta)
          !fta(11)=0.5d0*(1.0d0+beta)*(3.0d0*(alpha**2)-1.0d0)
          !fta(12)=-0.5d0*((beta**3)-beta)
          fta(11)=0.5d0*(1.0d0+beta)*(-3.0d0*(alpha**2)+1.0d0)
          fta(12)=-0.5d0*(-1.0d0*(beta**3)+beta)
          
          ftb(9)=-0.5d0*((alpha**3)-alpha)
          ftb(10)=0.5d0*(1.0d0+alpha)*(3.0d0*(beta**2)-1.0d0)
          !ftb(11)=0.5d0*((alpha**3)-alpha)
          !ftb(12)=0.5d0*(1.0d0-alpha)*(3.0d0*(beta**2)-1.0d0)
          ftb(11)=0.5d0*(-1.0d0*(alpha**3)+alpha)
          ftb(12)=0.5d0*(1.0d0-alpha)*(-3.0d0*(beta**2)+1.0d0)
          
          ftaa(9)=0.5d0*(1.0d0-beta)*(6.0*alpha); ftaa(10)=0.0d0; ftaa(11)=0.5d0*(1.0d0+beta)*(-6.0d0*alpha); ftaa(12)=0.0d0
          ftbb(9)=0.0d0; ftbb(10)=0.5d0*(1.0d0+alpha)*(6.0*beta); ftbb(11)=0.0d0; ftbb(12)=0.5d0*(1.0d0-alpha)*(-6.0d0*beta)
          ftab(9)=-0.5d0*(3.0d0*(alpha**2)-1.0d0); ftab(10)=0.5d0*(3.0d0*(beta**2)-1.0d0); ftab(11)=0.5d0*(-3.0d0*(alpha**2)+1.0d0); ftab(12)=-0.5d0*(-3.0d0*(beta**2)+1.0d0)
          
          ftaaa(9)=0.5d0*(1.0d0-beta)*(6.0); ftaaa(10)=0.0d0; ftaaa(11)=0.5d0*(1.0d0+beta)*(-6.0d0); ftaaa(12)=0.0d0
          ftbba(9)=0.0d0; ftbba(10)=0.5d0*(6.0*beta); ftbba(11)=0.0d0; ftbba(12)=-0.5d0*(-6.0d0*beta)
          ftaba(9)=-0.5d0*(6.0*alpha); ftaba(10)=0.0d0; ftaba(11)=0.5d0*(-6.0d0*alpha); ftaba(12)=0.0d0
          
          ftaab(9)=-0.5d0*(6.0*alpha); ftaab(10)=0.0d0; ftaab(11)=0.5d0*(-6.0d0*alpha); ftaab(12)=0.0d0
          ftbbb(9)=0.0d0; ftbbb(10)=0.5d0*(1.0d0+alpha)*(6.0); ftbbb(11)=0.0d0; ftbbb(12)=0.5d0*(1.0d0-alpha)*(-6.0d0)
          ftabb(9)=0.0d0; ftabb(10)=0.5d0*(6.0*beta); ftabb(11)=0.0d0; ftabb(12)=-0.5d0*(-6.0d0*beta)
          
          if (mexp.eq.12) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(13)=0.5d0*(1.0d0-beta)*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ft(14)=0.5d0*(1.0d0+alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(15)=0.5d0*(1.0d0+beta)*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ft(16)=0.5d0*(1.0d0-alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(17)=((alpha**2)-1.0d0)*((beta**2)-1.0d0) 
          
          fta(13)=0.5d0*(1.0d0-beta)*(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)
          fta(14)=0.5d0*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(15)=0.5d0*(1.0d0+beta)*(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)
          fta(16)=-0.5d0*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(17)=(2.0*alpha)*((beta**2)-1.0d0)
          
          ftb(13)=-0.5d0*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ftb(14)=0.5d0*(1.0d0+alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta)
          ftb(15)=0.5d0*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ftb(16)=0.5d0*(1.0d0-alpha)*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(17)=(2.0*beta)*((alpha**2)-1.0d0)
          
          ftaa(13)=0.5d0*(1.0d0-beta)*(12.0*(alpha**2)-(20.0d0/9.0d0)); ftaa(14)=0.0d0; ftaa(15)=0.5d0*(1.0d0+beta)*(12.0*(alpha**2)-(20.0d0/9.0d0)); ftaa(16)=0.0d0
          ftaa(17)=2.0d0*((beta**2)-1.0d0)
          ftbb(13)=0.0d0; ftbb(14)=0.5d0*(1.0d0+alpha)*(12.0*(beta**2)-(20.0d0/9.0d0)); ftbb(15)=0.0d0; ftbb(16)=0.5d0*(1.0d0-alpha)*(12.0*(beta**2)-(20.0d0/9.0d0))
          ftbb(17)=2.0d0*((alpha**2)-1.0d0)
          ftab(13)=-0.5d0*(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha); ftab(14)=0.5d0*(4.0*(beta**3)-(20.0d0/9.0d0)*beta); ftab(15)=0.5d0*(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha); 
          ftab(16)=-0.5d0*(4.0*(beta**3)-(20.0d0/9.0d0)*beta); ftab(17)=(2.0*alpha)*(2.0*beta)
          
          ftaaa(13)=0.5d0*(1.0d0-beta)*(24.0*alpha); ftaaa(14)=0.0d0; ftaaa(15)=0.5d0*(1.0d0+beta)*(24.0*alpha); ftaaa(16)=0.0d0
          ftaaa(17)=0.0d0
          ftbba(13)=0.0d0; ftbba(14)=0.5d0*(12.0*(beta**2)-(20.0d0/9.0d0)); ftbba(15)=0.0d0; ftbba(16)=-0.5d0*(12.0*(beta**2)-(20.0d0/9.0d0))
          ftbba(17)=2.0d0*(2.0*alpha)
          ftaba(13)=-0.5d0*(12.0*(alpha**2)-(20.0d0/9.0d0)); ftaba(14)=0.0d0; ftaba(15)=0.5d0*(12.0*(alpha**2)-(20.0d0/9.0d0)); 
          ftaba(16)=0.0d0; ftaba(17)=(2.0)*(2.0*beta)
          
          ftaab(13)=-0.5d0*(12.0*(alpha**2)-(20.0d0/9.0d0)); ftaab(14)=0.0d0; ftaab(15)=0.5d0*(12.0*(alpha**2)-(20.0d0/9.0d0)); ftaab(16)=0.0d0
          ftaab(17)=2.0d0*(2.0*beta)
          ftbbb(13)=0.0d0; ftbbb(14)=0.5d0*(1.0d0+alpha)*(24.0*beta); ftbbb(15)=0.0d0; ftbbb(16)=0.5d0*(1.0d0-alpha)*(24.0*beta)
          ftbbb(17)=0.0d0
          ftabb(13)=0.0d0; ftabb(14)=0.5d0*(12.0*(beta**2)-(20.0d0/9.0d0)); ftabb(15)=0.0d0; 
          ftabb(16)=-0.5d0*(12.0*(beta**2)-(20.0d0/9.0d0)); ftabb(17)=(2.0*alpha)*(2.0)
          
          if (mexp.eq.17) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(18)=0.5d0*(1.0d0-beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+(1.0d0/4.0d0)*alpha)
          ft(19)=0.5d0*(1.0d0+alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          !ft(20)=0.5d0*(1.0+beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+(1.0d0/4.0)*alpha)
          !ft(21)=0.5d0*(1.0d0-alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(20)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**5)+(5.0d0/4.0d0)*(alpha**3)-(1.0d0/4.0d0)*alpha)
          ft(21)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**5)+(5.0d0/4.0d0)*(beta**3)-(1.0d0/4.0d0)*beta)
          ft(22)=((alpha**2)-1.0d0)*((beta**3)-beta)
          ft(23)=((beta**2)-1.0d0)*((alpha**3)-alpha)  
          
          fta(18)=0.5d0*(1.0d0-beta)*(5.0d0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25d0)
          fta(19)=0.5d0*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(0.25d0*beta))
          !fta(20)=0.5d0*(1.0d0+beta)*(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)
          !fta(21)=-0.5d0*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(20)=0.5d0*(1.0d0+beta)*(-5.0d0*(alpha**4)+(15.0d0/4.0d0)*(alpha**2)-0.25d0)
          fta(21)=-0.5d0*(-1.0d0*(beta**5)+(5.0d0/4.0d0)*(beta**3)-(1.0d0/4.0d0)*beta)
          fta(22)=(2.0*alpha)*((beta**3)-beta)
          fta(23)=((beta**2)-1.0d0)*(3.0*(alpha**2)-1.0d0)
          
          ftb(18)=-0.5d0*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25d0*alpha)
          ftb(19)=0.5d0*(1.0d0+alpha)*(5.0d0*(beta**4)-(15.0/4.0)*(beta**2)+0.25d0)
          !ftb(20)=0.5d0*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          !ftb(21)=0.5d0*(1.0d0-alpha)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25d0)
          ftb(20)=0.5d0*(-1.0d0*(alpha**5)+(5.0d0/4.0d0)*(alpha**3)-(1.0d0/4.0d0)*alpha)
          ftb(21)=0.5d0*(1.0d0-alpha)*(-5.0d0*(beta**4)+(15.0d0/4.0d0)*(beta**2)-0.25d0)
          ftb(22)=((alpha**2)-1.0d0)*(3.0*(beta**2)-1.0d0)
          ftb(23)=(2.0*beta)*((alpha**3)-alpha)
          
          ftaa(18)=0.5d0*(1.0d0-beta)*(20.0d0*(alpha**3)-(30.0d0/4.0d0)*alpha)
          ftaa(19)=0.0d0
          ftaa(20)=0.5d0*(1.0d0+beta)*(-20.0d0*(alpha**3)+(30.0d0/4.0d0)*alpha)
          ftaa(21)=0.0d0
          ftaa(22)=2.0d0*((beta**3)-beta)
          ftaa(23)=((beta**2)-1.0d0)*(6.0d0*alpha)
          
          ftbb(18)=0.0d0
          ftbb(19)=0.5d0*(1.0d0+alpha)*(20.0d0*(beta**3)-(30.0d0/4.0d0)*beta)
          ftbb(20)=0.0d0
          ftbb(21)=0.5d0*(1.0d0-alpha)*(-20.0d0*(beta**3)+(30.0d0/4.0d0)*beta)
          ftbb(22)=((alpha**2)-1.0d0)*(6.0d0*beta)
          ftbb(23)=2.0d0*((alpha**3)-alpha)
          
          ftab(18)=-0.5d0*(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25d0)
          ftab(19)=0.5d0*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25d0)
          ftab(20)=0.5d0*(-5.0d0*(alpha**4)+(15.0d0/4.0d0)*(alpha**2)-0.25d0)
          ftab(21)=-0.5d0*(-5.0d0*(beta**4)+(15.0d0/4.0d0)*(beta**2)-0.25d0)
          ftab(22)=(2.0d0*alpha)*(3.0d0*(beta**2)-1.0d0)
          ftab(23)=(2.0d0*beta)*(3.0d0*(alpha**2)-1.0d0)
          
          ftaaa(18)=0.5d0*(1.0d0-beta)*(60.0d0*(alpha**2)-(30.0d0/4.0d0)); ftaaa(19)=0.0d0; ftaaa(20)=0.5d0*(1.0d0+beta)*(-60.0d0*(alpha**2)+(30.0d0/4.0d0))
          ftaaa(21)=0.0d0; ftaaa(22)=0.0d0; ftaaa(23)=((beta**2)-1.0d0)*(6.0d0)
          ftbba(18)=0.0d0; ftbba(19)=0.5d0*(20.0*(beta**3)-(30.0d0/4.0d0)*beta); ftbba(20)=0.0d0
          ftbba(21)=-0.5d0*(-20.0d0*(beta**3)+(30.0d0/4.0d0)*beta); ftbba(22)=(2.0*alpha)*(6.0*beta); ftbba(23)=2.0*(3.0*(alpha**2)-1.0d0)
          ftaba(18)=-0.5d0*(20.0d0*(alpha**3)-(30.0d0/4.0d0)*(alpha)); ftaba(19)=0.0d0; ftaba(20)=0.5d0*(-20.0d0*(alpha**3)+(30.0d0/4.0d0)*(alpha))
          ftaba(21)=0.0d0; ftaba(22)=(2.0d0)*(3.0d0*(beta**2)-1.0d0); ftaba(23)=(2.0d0*beta)*(6.0d0*alpha)
          
          ftaab(18)=-0.5d0*(20.0*(alpha**3)-(30.0d0/4.0d0)*alpha); ftaab(19)=0.0d0; ftaab(20)=0.5d0*(-20.0d0*(alpha**3)+(30.0d0/4.0d0)*alpha)
          ftaab(21)=0.0d0; ftaab(22)=2.0d0*(3.0*(beta**2)-1.0d0); ftaab(23)=(2.0d0*beta)*(6.0d0*alpha)
          ftbbb(18)=0.0d0; ftbbb(19)=0.5d0*(1.0d0+alpha)*(60.0*(beta**2)-(30.0/4.0)); ftbbb(20)=0.0d0
          ftbbb(21)=0.5d0*(1.0d0-alpha)*(-60.0d0*(beta**2)+(30.0d0/4.0d0)); ftbbb(22)=((alpha**2)-1.0d0)*(6.0d0); ftbbb(23)=0.0d0
          ftabb(18)=0.0d0; ftabb(19)=0.5d0*(20.0*(beta**3)-(30.0d0/4.0d0)*beta); ftabb(20)=0.0d0
          ftabb(21)=-0.5d0*(-20.0d0*(beta**3)+(30.0d0/4.0d0)*beta); ftabb(22)=(2.0d0*alpha)*(6.0d0*beta); ftabb(23)=(2.0)*(3.0*(alpha**2)-1.0d0)
          
          if (mexp.eq.23) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(24)=0.5d0*(1.0d0-beta)*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ft(25)=0.5d0*(1.0d0+alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(26)=0.5d0*(1.0d0+beta)*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ft(27)=0.5d0*(1.0d0-alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(28)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**2)-1.0d0)
          ft(29)=((alpha**3)-alpha)*((beta**3)-beta)
          ft(30)=((alpha**2)-1.0d0)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0)) 
          
          fta(24)=0.5d0*(1.0d0-beta)*(6.0*(alpha**5)-(28.0d0/5.0d0)*(alpha**3)+(518.0d0/625.0d0)*alpha)
          fta(25)=0.5d0*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(26)=0.5d0*(1.0d0+beta)*(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)
          fta(27)=-0.5d0*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(28)=(4.0d0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**2)-1.0d0)
          fta(29)=(3.0d0*(alpha**2)-1.0d0)*((beta**3)-beta)
          fta(30)=(2.0d0*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          
          ftb(24)=-0.5d0*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ftb(25)=0.5d0*(1.0d0+alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(26)=0.5d0*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ftb(27)=0.5d0*(1.0d0-alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(28)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(2.0*beta)
          ftb(29)=((alpha**3)-alpha)*(3.0*(beta**2)-1.0d0)
          ftb(30)=((alpha**2)-1.0d0)*(4.0*(beta**3)-(20.0/9.0)*beta)
          
          ftaa(24)=0.5d0*(1.0d0-beta)*(30.0*(alpha**4)-(84.0d0/5.0d0)*(alpha**2)+(518.0d0/625.0d0))
          ftaa(25)=0.0d0
          ftaa(26)=0.5d0*(1.0d0+beta)*(30.0*(alpha**4)-(84.0d0/5.0d0)*(alpha**2)+(518.0d0/625.0d0))
          ftaa(27)=0.0d0
          ftaa(28)=(12.0*(alpha**2)-(20.0d0/9.0d0))*((beta**2)-1.0d0)
          ftaa(29)=(6.0*alpha)*((beta**3)-beta)
          ftaa(30)=2.0*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          
          ftbb(24)=0.0d0
          ftbb(25)=0.5d0*(1.0d0+alpha)*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0))
          ftbb(26)=0.0d0
          ftbb(27)=0.5d0*(1.0d0-alpha)*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0))
          ftbb(28)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*2.0
          ftbb(29)=((alpha**3)-alpha)*(6.0*beta)
          ftbb(30)=((alpha**2)-1.0d0)*(12.0*(beta**2)-(20.0/9.0))
          
          ftab(24)=-0.5d0*(6.0*(alpha**5)-(28.0d0/5.0d0)*(alpha**3)+(518.0d0/625.0d0)*alpha)
          ftab(25)=0.5d0*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftab(26)=0.5d0*(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)
          ftab(27)=-0.5d0*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftab(28)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*(2.0*beta)
          ftab(29)=(3.0*(alpha**2)-1.0d0)*(3.0*(beta**2)-1.0d0)
          ftab(30)=(2.0*alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta)
          
          ftaaa(24)=0.5d0*(1.0d0-beta)*(120.0*(alpha**3)-(168.0d0/5.0d0)*(alpha)); ftaaa(25)=0.0d0
          ftaaa(26)=0.5d0*(1.0d0+beta)*(120.0*(alpha**3)-(168.0d0/5.0d0)*(alpha)); ftaaa(27)=0.0d0
          ftaaa(28)=(24.0*(alpha))*((beta**2)-1.0d0); ftaaa(29)=(6.0)*((beta**3)-beta); ftaaa(30)=0.0d0
          
          ftbba(24)=0.0d0; ftbba(25)=0.5d0*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0))
          ftbba(26)=0.0d0; ftbba(27)=-0.5d0*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0))
          ftbba(28)=(4.0*(alpha**3)-(20.0d0/9.0d0)*(alpha))*2.0; ftbba(29)=(3.0*(alpha**2)-1.0d0)*(6.0*beta)
          ftbba(30)=(2.0*alpha)*(12.0*(beta**2)-(20.0/9.0))
          
          ftaba(24)=-0.5d0*(30.0*(alpha**4)-(84.0d0/5.0d0)*(alpha**2)+(518.0d0/625.0d0)); ftaba(25)=0.0d0
          ftaba(26)=0.5d0*(30.0*(alpha**4)-(84.0/5.0)*(alpha**2)+(518.0d0/625.0d0)); ftaba(27)=0.0d0
          ftaba(28)=(12.0*(alpha**2)-(20.0d0/9.0d0))*(2.0*beta); ftaba(29)=(6.0*alpha)*(3.0*(beta**2)-1.0d0); ftaba(30)=(2.0)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta)
          
          ftaab(24)=-0.5d0*(30.0*(alpha**4)-(84.0d0/5.0d0)*(alpha**2)+(518.0d0/625.0d0)); ftaab(25)=0.0d0
          ftaab(26)=0.5d0*(30.0*(alpha**4)-(84.0d0/5.0d0)*(alpha**2)+(518.0d0/625.0d0)); ftaab(27)=0.0d0
          ftaab(28)=(12.0*(alpha**2)-(20.0d0/9.0d0))*(2.0*beta); ftaab(29)=(6.0*alpha)*(3.0*(beta**2)-1.0d0)
          ftaab(30)=2.0*(4.0*(beta**3)-(20.0d0/9.0d0)*(beta))
          
          ftbbb(24)=0.0d0; ftbbb(25)=0.5d0*(1.0d0+alpha)*(120.0*(beta**3)-(168.0d0/5.0d0)*(beta)); ftbbb(26)=0.0d0
          ftbbb(27)=0.5d0*(1.0d0-alpha)*(120.0*(beta**3)-(168.0d0/5.0d0)*(beta)); ftbbb(28)=0.0d0
          ftbbb(29)=((alpha**3)-alpha)*(6.0); ftbbb(30)=((alpha**2)-1.0d0)*(24.0*beta)
          
          ftabb(24)=0.0; ftabb(25)=0.5d0*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0)); ftabb(26)=0.0
          ftabb(27)=-0.5d0*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0)); ftabb(28)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*(2.0)
          ftabb(29)=(3.0*(alpha**2)-1.0d0)*(6.0*beta); ftabb(30)=(2.0*alpha)*(12.0*(beta**2)-(20.0d0/9.0d0))
          
          if (mexp.eq.30) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(31)=0.5d0*(1.0d0-beta)*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          ft(32)=0.5d0*(1.0d0+alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          !ft(33)=0.5d0*(1.0d0+beta)*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          !ft(34)=0.5d0*(1.0d0-alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(33)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**7)+(14.0d0/9.0d0)*(alpha**5)-(49.0d0/81.0d0)*(alpha**3)+(4.0d0/81.0d0)*alpha)
          ft(34)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**7)+(14.0d0/9.0d0)*(beta**5)-(49.0d0/81.0d0)*(beta**3)+(4.0d0/81.0d0)*beta) 
          ft(35)=((alpha**2)-1.0d0)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+0.25*beta)
          ft(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**3)-beta) 
          ft(37)=((alpha**3)-alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(38)=((beta**2)-1.0d0)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          
          fta(31)=0.5d0*(1.0d0-beta)*(7.0d0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))
          fta(32)=0.5d0*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          !fta(33)=0.5d0*(1.0+beta)*(7.0*(alpha**6)-(70.0/9.0)*(alpha**4)+(147.0/81.0d0)*(alpha**2)-(4.0/81.0d0))
          !fta(34)=-0.5d0*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          fta(33)=0.5d0*(1.0d0+beta)*(-7.0d0*(alpha**6)+(70.0/9.0)*(alpha**4)-(147.0/81.0d0)*(alpha**2)+(4.0/81.0d0))
          fta(34)=-0.5d0*(-1.0d0*(beta**7)+(14.0d0/9.0d0)*(beta**5)-(49.0d0/81.0d0)*(beta**3)+(4.0d0/81.0d0)*beta)
          
          fta(35)=(2.0*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+0.25*beta)
          fta(36)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**3)-beta)
          fta(37)=(3.0*(alpha**2)-1.0d0)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(38)=((beta**2)-1.0d0)*(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)
          
          ftb(31)=-0.5d0*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          ftb(32)=0.5d0*(1.0d0+alpha)*(7.0d0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0))   
          !ftb(33)=0.5d0*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          !ftb(34)=0.5d0*(1.0-alpha)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftb(33)=0.5d0*(-1.0d0*(alpha**7)+(14.0d0/9.0d0)*(alpha**5)-(49.0d0/81.0d0)*(alpha**3)+(4.0d0/81.0d0)*alpha)
          ftb(34)=0.5d0*(1.0d0-alpha)*(-7.0d0*(beta**6)+(70.0d0/9.0d0)*(beta**4)-(147.0d0/81.0d0)*(beta**2)+(4.0d0/81.0d0)) 
          
          ftb(35)=((alpha**2)-1.0d0)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25)  
          ftb(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(3.0*(beta**2)-1.0d0) 
          ftb(37)=((alpha**3)-alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta)   
          ftb(38)=(2.0*beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          
          ftaa(31)=0.5d0*(1.0d0-beta)*(42.0d0*(alpha**5)-(280.0d0/9.0d0)*(alpha**3)+(294.0d0/81.0d0)*alpha)
          ftaa(32)=0.0d0
          ftaa(33)=0.5d0*(1.0d0+beta)*(-42.0d0*(alpha**5)+(280.0d0/9.0d0)*(alpha**3)-(294.0d0/81.0d0)*alpha)
          ftaa(34)=0.0d0
          ftaa(35)=2.0d0*((beta**5)-(5.0d0/4.0d0)*(beta**3)+0.25*beta)
          ftaa(36)=(12.0d0*(alpha**2)-(20.0d0/9.0d0))*((beta**3)-beta)
          ftaa(37)=(6.0d0*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ftaa(38)=((beta**2)-1.0d0)*(20.0d0*(alpha**3)-(30.0d0/4.0d0)*alpha)
          
          ftbb(31)=0.0d0
          ftbb(32)=0.5d0*(1.0d0+alpha)*(42.0d0*(beta**5)-(280.0d0/9.0d0)*(beta**3)+(294.0d0/81.0d0)*beta)   
          ftbb(33)=0.0d0
          ftbb(34)=0.5d0*(1.0d0-alpha)*(-42.0d0*(beta**5)+(280.0d0/9.0d0)*(beta**3)-(294.0d0/81.0d0)*beta)  
          ftbb(35)=((alpha**2)-1.0d0)*(20.0*(beta**3)-(30.0d0/4.0d0)*beta)  
          ftbb(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(6.0*beta) 
          ftbb(37)=((alpha**3)-alpha)*(12.0*(beta**2)-(20.0d0/9.0d0))   
          ftbb(38)=2.0d0*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          
          ftab(31)=-0.5d0*(7.0d0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))
          ftab(32)=0.5d0*(7.0d0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0))
          ftab(33)=0.5d0*(-7.0d0*(alpha**6)+(70.0d0/9.0d0)*(alpha**4)-(147.0d0/81.0d0)*(alpha**2)+(4.0d0/81.0d0))
          ftab(34)=-0.5d0*(-7.0d0*(beta**6)+(70.0d0/9.0d0)*(beta**4)-(147.0d0/81.0d0)*(beta**2)+(4.0d0/81.0d0))
          ftab(35)=(2.0d0*alpha)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25d0)
          ftab(36)=(4.0d0*(alpha**3)-(20.0d0/9.0d0)*alpha)*(3.0d0*(beta**2)-1.0d0)
          ftab(37)=(3.0d0*(alpha**2)-1.0d0)*(4.0d0*(beta**3)-(20.0d0/9.0d0)*beta)
          ftab(38)=(2.0d0*beta)*(5.0d0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25d0)
          
          
          
          
          ftaaa(31)=0.5d0*(1.0d0-beta)*(210.0*(alpha**4)-(840.0d0/9.0d0)*(alpha**2)+(294.0d0/81.0d0)); ftaaa(32)=0.0d0
          ftaaa(33)=0.5d0*(1.0d0+beta)*(-210.0d0*(alpha**4)+(840.0d0/9.0d0)*(alpha**2)-(294.0d0/81.0d0)); ftaaa(34)=0.0d0
          ftaaa(35)=0.0d0; ftaaa(36)=(24.0d0*alpha)*((beta**3)-beta); ftaaa(37)=(6.0d0)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ftaaa(38)=((beta**2)-1.0d0)*(60.0*(alpha**2)-(30.0d0/4.0d0))
          
          ftbba(31)=0.0d0; ftbba(32)=0.5d0*(42.0d0*(beta**5)-(280.0d0/9.0d0)*(beta**3)+(294.0d0/81.0d0)*beta)   
          ftbba(33)=0.0d0; ftbba(34)=-0.5d0*(-42.0d0*(beta**5)+(280.0d0/9.0d0)*(beta**3)-(294.0d0/81.0d0)*beta)  
          ftbba(35)=(2.0d0*alpha)*(20.0*(beta**3)-(30.0d0/4.0d0)*beta); ftbba(36)=(4.0d0*(alpha**3)-(20.0d0/9.0d0)*(alpha))*(6.0d0*beta) 
          ftbba(37)=(3.0d0*(alpha**2)-1.0d0)*(12.0*(beta**2)-(20.0d0/9.0d0)); ftbba(38)=2.0*(5.0d0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)
          
          ftaba(31)=-0.5d0*(42.0*(alpha**5)-(280.0d0/9.0d0)*(alpha**3)+(294.0d0/81.0d0)*alpha); ftaba(32)=0.0d0
          ftaba(33)=0.5d0*(-42.0*(alpha**5)+(280.0d0/9.0d0)*(alpha**3)-(294.0d0/81.0d0)*alpha); ftaba(34)=0.0d0
          ftaba(35)=(2.0d0)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25); ftaba(36)=(12.0*(alpha**2)-(20.0d0/9.0d0))*(3.0*(beta**2)-1.0d0)
          ftaba(37)=(6.0*alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta); ftaba(38)=(2.0*beta)*(20.0*(alpha**3)-(30.0d0/4.0d0)*(alpha))
          
          ftaab(31)=-0.5d0*(42.0*(alpha**5)-(280.0d0/9.0d0)*(alpha**3)+(294.0d0/81.0d0)*alpha)
          ftaab(32)=0.0d0
          ftaab(33)=0.5d0*(-42.0*(alpha**5)+(280.0d0/9.0d0)*(alpha**3)-(294.0d0/81.0d0)*alpha)
          ftaab(34)=0.0d0
          ftaab(35)=2.0*(5.0d0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25)
          ftaab(36)=(12.0*(alpha**2)-(20.0d0/9.0d0))*(3.0d0*(beta**2)-1.0d0)
          ftaab(37)=(6.0*alpha)*(4.0d0*(beta**3)-(20.0d0/9.0d0)*(beta))
          ftaab(38)=(2.0d0*beta)*(20.0*(alpha**3)-(30.0d0/4.0d0)*alpha)
          
          ftbbb(31)=0.0d0; ftbbb(32)=0.5d0*(1.0d0+alpha)*(210.0*(beta**4)-(840.0d0/9.0d0)*(beta**2)+(294.0d0/81.0d0))  
          ftbbb(33)=0.0d0; ftbbb(34)=0.5d0*(1.0d0-alpha)*(-210.0*(beta**4)+(840.0d0/9.0d0)*(beta**2)-(294.0d0/81.0d0))  
          ftbbb(35)=((alpha**2)-1.0d0)*(60.0*(beta**2)-(30.0d0/4.0d0)); ftbbb(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(6.0d0) 
          ftbbb(37)=((alpha**3)-alpha)*(24.0d0*beta); ftbbb(38)=0.0d0
          
          ftabb(31)=0.0d0; ftabb(32)=0.5d0*(42.0*(beta**5)-(280.0d0/9.0d0)*(beta**3)+(294.0d0/81.0d0)*beta)
          ftabb(33)=0.0d0; ftabb(34)=-0.5d0*(-42.0*(beta**5)+(280.0d0/9.0d0)*(beta**3)-(294.0d0/81.0d0)*beta)
          ftabb(35)=(2.0*alpha)*(20.0*(beta**3)-(30.0d0/4.0d0)*(beta)); ftabb(36)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*(6.0d0*beta)
          ftabb(37)=(3.0*(alpha**2)-1.0d0)*(12.0*(beta**2)-(20.0d0/9.0d0)); ftabb(38)=(2.0d0)*(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)
          
          if (mexp.eq.38) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(39)=0.5d0*(1.0d0-beta)*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ft(40)=0.5d0*(1.0d0+alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(41)=0.5d0*(1.0d0+beta)*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ft(42)=0.5d0*(1.0d0-alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(43)=((alpha**2)-1.0d0)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(44)=((alpha**3)-alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(45)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(46)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**3)-beta)
          ft(47)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**2)-1.0d0)
          
          fta(39)=0.5d0*(1.0d0-beta)*(8.0d0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))
          fta(40)=0.5d0*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          fta(41)=0.5d0*(1.0d0+beta)*(8.0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))
          fta(42)=-0.5d0*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          fta(43)=(2.0*alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(44)=(3.0*(alpha**2)-1.0d0)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(45)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(46)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*((beta**3)-beta)
          fta(47)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*((beta**2)-1.0d0)
          
          ftb(39)=-0.5d0*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ftb(40)=0.5d0*(1.0d0+alpha)*(8.0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftb(41)=0.5d0*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ftb(42)=0.5d0*(1.0d0-alpha)*(8.0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftb(43)=((alpha**2)-1.0d0)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(44)=((alpha**3)-alpha)*(5.0*(beta**4)-(15.0/4.0)*(beta**2)+0.25)
          ftb(45)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(46)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(3.0*(beta**2)-1.0d0)
          ftb(47)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(2.0*beta)
          
          ftaa(39)=0.5d0*(1.0d0-beta)*(56.0d0*(alpha**6)-(360.0d0/7.0d0)*(alpha**4)+(3384.0d0/343.0d0)*(alpha**2)-(25832.0d0/117649.0d0))
          ftaa(40)=0.0d0
          ftaa(41)=0.5d0*(1.0d0+beta)*(56.0d0*(alpha**6)-(360.0d0/7.0d0)*(alpha**4)+(3384.0d0/343.0d0)*(alpha**2)-(25832.0d0/117649.0d0))
          ftaa(42)=0.0d0
          ftaa(43)=2.0*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ftaa(44)=(6.0*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ftaa(45)=(12.0*(alpha**2)-(20.0d0/9.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ftaa(46)=(20.0*(alpha**3)-(30.0d0/4.0d0)*alpha)*((beta**3)-beta)
          ftaa(47)=(30.0*(alpha**4)-(84.0/5.0)*(alpha**2)+(518.0d0/625.0d0))*((beta**2)-1.0d0)
          
          ftbb(39)=0.0d0
          ftbb(40)=0.5d0*(1.0d0+alpha)*(56.0*(beta**6)-(360.0d0/7.0d0)*(beta**4)+(3384.0d0/343.0d0)*(beta**2)-(25832.0d0/117649.0d0))
          ftbb(41)=0.0d0
          ftbb(42)=0.5d0*(1.0d0-alpha)*(56.0*(beta**6)-(360.0d0/7.0d0)*(beta**4)+(3384.0d0/343.0d0)*(beta**2)-(25832.0d0/117649.0d0))
          ftbb(43)=((alpha**2)-1.0d0)*(30.0*(beta**4)-(84.0/5.0)*(beta**2)+(518.0d0/625.0d0))
          ftbb(44)=((alpha**3)-alpha)*(20.0*(beta**3)-(30.0/4.0)*beta)
          ftbb(45)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(12.0*(beta**2)-(20.0/9.0))
          ftbb(46)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(6.0*beta)
          ftbb(47)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*2.0
          
          ftab(39)=-0.5d0*(8.0d0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))
          ftab(40)=0.5d0*(8.0d0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftab(41)=0.5d0*(8.0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))
          ftab(42)=-0.5d0*(8.0d0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftab(43)=(2.0*alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftab(44)=(3.0*(alpha**2)-1.0d0)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+(1.0d0/4.0d0))
          ftab(45)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*(beta))
          ftab(46)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*(3.0*(beta**2)-1.0d0)
          ftab(47)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*(2.0*beta)
          
          if (mexp.eq.47) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(48)=0.5d0*(1.0d0-beta)*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          ft(49)=0.5d0*(1.0d0+alpha)*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          !ft(50)=0.5d0*(1.0d0+beta)*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          !ft(51)=0.5d0*(1.0d0-alpha)*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          ft(50)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**9)+(15.0d0/8.0d0)*(alpha**7)-(273.0d0/256.0d0)*(alpha**5)+(205.0d0/1024.0d0)*(alpha**3)-(9.0d0/1024.0d0)*alpha)
          ft(51)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**9)+(15.0d0/8.0d0)*(beta**7)-(273.0d0/256.0d0)*(beta**5)+(205.0d0/1024.0d0)*(beta**3)-(9.0d0/1024.0d0)*beta)
          ft(52)=((alpha**2)-1.0d0)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(53)=((alpha**3)-alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(54)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(55)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(56)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**3)-beta)
          ft(57)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*((beta**2)-1.0d0)
          
          fta(48)=0.5d0*(1.0d0-beta)*(9.0d0*(alpha**8)-(105.0d0/8.0d0)*(alpha**6)+(1365.0d0/256.0d0)*(alpha**4)-(615.0d0/1024.0d0)*(alpha**2)+(9.0d0/1024.0d0))
          fta(49)=0.5d0*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          !fta(50)=0.5d0*(1.0d0+beta)*(9.0d0*(alpha**8)-(105.0d0/8.0d0)*(alpha**6)+(1365.0d0/256.0d0)*(alpha**4)-(615.0d0/1024.0d0)*(alpha**2)+(9.0d0/1024.0d0))
          !fta(51)=-0.5d0*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          fta(50)=0.5d0*(1.0d0+beta)*(-9.0d0*(alpha**8)+(105.0d0/8.0d0)*(alpha**6)-(1365.0d0/256.0d0)*(alpha**4)+(615.0d0/1024.0d0)*(alpha**2)-(9.0d0/1024.0d0))
          fta(51)=-0.5d0*(-1.0d0*(beta**9)+(15.0d0/8.0d0)*(beta**7)-(273.0d0/256.0d0)*(beta**5)+(205.0d0/1024.0d0)*(beta**3)-(9.0d0/1024.0d0)*beta)
          fta(52)=(2.0*alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          fta(53)=(3.0*(alpha**2)-1.0d0)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(54)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(55)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(56)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*((beta**3)-beta)
          fta(57)=(7.0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))*((beta**2)-1.0d0)
          
          ftb(48)=-0.5d0*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          ftb(49)=0.5d0*(1.0d0+alpha)*(9.0d0*(beta**8)-(105.0d0/8.0d0)*(beta**6)+(1365.0d0/256.0d0)*(beta**4)-(615.0d0/1024.0d0)*(beta**2)+(9.0d0/1024.0d0))
          !ftb(50)=0.5d0*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          !ftb(51)=0.5d0*(1.0d0-alpha)*(9.0d0*(beta**8)-(105.0d0/8.0d0)*(beta**6)+(1365.0d0/256.0d0)*(beta**4)-(615.0d0/1024.0d0)*(beta**2)+(9.0d0/1024.0d0))
          ftb(50)=0.5d0*(-1.0d0*(alpha**9)+(15.0d0/8.0d0)*(alpha**7)-(273.0d0/256.0d0)*(alpha**5)+(205.0d0/1024.0d0)*(alpha**3)-(9.0d0/1024.0d0)*alpha)
          ftb(51)=0.5d0*(1.0d0-alpha)*(-9.0d0*(beta**8)+(105.0d0/8.0d0)*(beta**6)-(1365.0d0/256.0d0)*(beta**4)+(615.0d0/1024.0d0)*(beta**2)-(9.0d0/1024.0d0))
          ftb(52)=((alpha**2)-1.0d0)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftb(53)=((alpha**3)-alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(54)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(5.0*(beta**4)-(15.0/4.0)*(beta**2)+0.25)
          ftb(55)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(56)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(3.0*(beta**2)-1.0d0)
          ftb(57)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*(2.0d0*beta)
          
          ftaa(48)=0.5d0*(1.0d0-beta)*(72.0d0*(alpha**7)-(630.0d0/8.0d0)*(alpha**5)+(5460.0d0/256.0d0)*(alpha**3)-(1230.0d0/1024.0d0)*(alpha))
          ftaa(49)=0.0d0
          ftaa(50)=0.5d0*(1.0d0+beta)*(-72.0d0*(alpha**7)+(630.0d0/8.0d0)*(alpha**5)-(5460.0d0/256.0d0)*(alpha**3)+(1230.0d0/1024.0d0)*(alpha))
          ftaa(51)=0.0d0
          ftaa(52)=(2.0d0)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ftaa(53)=(6.0d0*alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ftaa(54)=(12.0d0*(alpha**2)-(20.0d0/9.0d0))*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ftaa(55)=(20.0*(alpha**3)-(30.0d0/4.0d0)*(alpha))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ftaa(56)=(30.0*(alpha**4)-(84.0d0/5.0d0)*(alpha**2)+(518.0d0/625.0d0))*((beta**3)-beta)
          ftaa(57)=(42.0*(alpha**5)-(280.0d0/9.0d0)*(alpha**3)+(294.0d0/81.0d0)*(alpha))*((beta**2)-1.0d0)
          
          ftbb(48)=0.0d0
          ftbb(49)=0.5d0*(1.0d0+alpha)*(72.0d0*(beta**7)-(630.0d0/8.0d0)*(beta**5)+(5460.0d0/256.0d0)*(beta**3)-(1230.0d0/1024.0d0)*(beta))
          ftbb(50)=0.0d0
          ftbb(51)=0.5d0*(1.0d0-alpha)*(-72.0d0*(beta**7)+(630.0d0/8.0d0)*(beta**5)-(5460.0d0/256.0d0)*(beta**3)+(1230.0d0/1024.0d0)*(beta))
          ftbb(52)=((alpha**2)-1.0d0)*(42.0*(beta**5)-(280.0d0/9.0d0)*(beta**3)+(294.0d0/81.0d0)*(beta)) 
          ftbb(53)=((alpha**3)-alpha)*(30.0*(beta**4)-(84.0d0/5.0d0)*(beta**2)+(518.0d0/625.0d0))
          ftbb(54)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(20.0*(beta**3)-(30.0/4.0)*(beta))
          ftbb(55)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(12.0*(beta**2)-(20.0/9.0))
          ftbb(56)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(6.0d0*beta)
          ftbb(57)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*(2.0d0)
          
          ftab(48)=-0.5d0*(9.0d0*(alpha**8)-(105.0d0/8.0d0)*(alpha**6)+(1365.0d0/256.0d0)*(alpha**4)-(615.0d0/1024.0d0)*(alpha**2)+(9.0d0/1024.0d0))
          ftab(49)=0.5d0*(9.0d0*(beta**8)-(105.0d0/8.0d0)*(beta**6)+(1365.0d0/256.0d0)*(beta**4)-(615.0d0/1024.0d0)*(beta**2)+(9.0d0/1024.0d0))
          ftab(50)=0.5d0*(-9.0d0*(alpha**8)+(105.0d0/8.0d0)*(alpha**6)-(1365.0d0/256.0d0)*(alpha**4)+(615.0d0/1024.0d0)*(alpha**2)-(9.0d0/1024.0d0))
          ftab(51)=-0.5d0*(-9.0d0*(beta**8)+(105.0d0/8.0d0)*(beta**6)-(1365.0d0/256.0d0)*(beta**4)+(615.0d0/1024.0d0)*(beta**2)-(9.0d0/1024.0d0))
          ftab(52)=(2.0*alpha)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftab(53)=(3.0*(alpha**2)-1.0d0)*(6.0d0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*(beta))
          ftab(54)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*(5.0d0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+(1.0d0/4.0d0))
          ftab(55)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*(4.0d0*(beta**3)-(20.0d0/9.0d0)*(beta))
          ftab(56)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*(3.0d0*(beta**2)-1.0d0)
          ftab(57)=(7.0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))*(2.0d0*beta)
          
          if (mexp.eq.57) then
              GOTO 2002
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          
          ft(58)=0.5*(1.0d0-beta)*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ft(59)=0.5*(1.0d0+alpha)*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          ft(60)=0.5*(1.0d0+beta)*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ft(61)=0.5*(1.0d0-alpha)*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          ft(62)=((alpha**2)-1.0d0)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(63)=((alpha**3)-alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(64)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(65)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(66)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(67)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*((beta**3)-beta)
          ft(68)=((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))*((beta**2)-1.0d0)
          
          fta(58)=0.5*(1.0d0-beta)*(10.0d0*(alpha**9)-(440.0d0/27.0d0)*(alpha**7)+(5825.0d0/729.0d0)*(alpha**5)-(691240.0d0/531441.0d0)*(alpha**3)+(234938.0d0/4782969.0d0)*alpha)
          fta(59)=0.5*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          fta(60)=0.5*(1.0d0+beta)*(10.0d0*(alpha**9)-(440.0d0/27.0d0)*(alpha**7)+(5825.0d0/729.0d0)*(alpha**5)-(691240.0d0/531441.0d0)*(alpha**3)+(234938.0d0/4782969.0d0)*alpha)
          fta(61)=-0.5*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          fta(62)=(2.0*alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          fta(63)=(3.0*(alpha**2)-1.0d0)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          fta(64)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(65)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(66)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(67)=(7.0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))*((beta**3)-beta)
          fta(68)=(8.0d0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))*((beta**2)-1.0d0)
          
          ftb(58)=-0.5*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ftb(59)=0.5*(1.0d0+alpha)*(10.0d0*(beta**9)-(440.0d0/27.0d0)*(beta**7)+(5825.0d0/729.0d0)*(beta**5)-(691240.0d0/531441.0d0)*(beta**3)+(234938.0d0/4782969.0d0)*beta)
          ftb(60)=0.5*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ftb(61)=0.5*(1.0d0-alpha)*(10.0d0*(beta**9)-(440.0d0/27.0d0)*(beta**7)+(5825.0d0/729.0d0)*(beta**5)-(691240.0d0/531441.0d0)*(beta**3)+(234938.0d0/4782969.0d0)*beta)
          ftb(62)=((alpha**2)-1.0d0)*(8.0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftb(63)=((alpha**3)-alpha)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftb(64)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(65)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(5.0*(beta**4)-(15.0/4.0)*(beta**2)+0.25)
          ftb(66)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(67)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*(3.0*(beta**2)-1.0d0)
          ftb(68)=((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))*(2.0*beta)
          
          if (mexp.eq.68) then
              GOTO 2002
          endif
    
2002   	do i1=1,4
            nodenum=elemcon(ele,i1+1,cs)
            xa=xa+fta(i1)*cscord(nodenum,2,cs)
            xb=xb+ftb(i1)*cscord(nodenum,2,cs)
            za=za+fta(i1)*cscord(nodenum,3,cs)
            zb=zb+ftb(i1)*cscord(nodenum,3,cs)
            
            xab=xab+ftab(i1)*cscord(nodenum,2,cs)
            zab=zab+ftab(i1)*cscord(nodenum,3,cs)
        enddo
   
    sljacdet=xa*zb-za*xb

    do i1=1,mexp
        ftx(i1)=(zb*fta(i1)-za*ftb(i1))/sljacdet
        ftz(i1)=(xa*ftb(i1)-xb*fta(i1))/sljacdet
        
        afun=(sljacdet*ftbb(i1)*(za**2))+(2.0*ftb(i1)*za*(zb*xa*zab-za*zb*xab))
        bfun=zb*(sljacdet*zb*ftaa(i1)-2.0*sljacdet*za*ftab(i1)+2.0*xab*za*zb*fta(i1)-2.0*zab*za*xb*fta(i1))
        ftxx(i1)=(afun+bfun)/(sljacdet**3)
        
        cfun=(sljacdet*ftbb(i1)*(xa**2))+(2.0*ftb(i1)*xa*(xb*xa*zab-za*xb*xab))
        dfun=xb*(sljacdet*xb*ftaa(i1)-2.0*sljacdet*xa*ftab(i1)+2.0*xab*xa*zb*fta(i1)-2.0*zab*xa*xb*fta(i1))
        ftzz(i1)=(cfun+dfun)/(sljacdet**3)
        
        efun=(sljacdet*ftab(i1)*xb*za)-(sljacdet*fta(i1)*xb*zab)-(sljacdet*ftaa(i1)*xb*zb)+(xb*(xab*za-xa*zab)*(ftb(i1)*za-fta(i1)*zb))
        ffun=(sljacdet*ftbb(i1)*xa*za)+(sljacdet*ftb(i1)*xa*zab)-(sljacdet*ftab(i1)*xa*zb)+(xa*(xb*zab-xab*zb)*(ftb(i1)*za-fta(i1)*zb))
        ftxz(i1)=(efun-ffun)/(sljacdet**3)
        
        pfun=(za**2)*(ftbb(i1)*(xa*zab-za*xab)+sljacdet*ftbba(i1))+2.0d0*za*zab*xa*(ftab(i1)*xb+ftb(i1)*xab)-2.0d0*za*za*xab*(ftab(i1)*zb+ftb(i1)*zab)
        qfun=(xa*zab-za*xab)*(ftaa(i1)*zb*zb-2.0d0*za*zb*ftab(i1))+zb*sljacdet*(ftaaa(i1)*zb+2.0d0*ftaa(i1)*zab)-2.0d0*sljacdet*za*(ftaba(i1)*zb+ftab(i1)*zab) &
              +2.0d0*xab*za*zb*(ftaa(i1)*zb+2.0d0*fta(i1)*zab)-2.0d0*za*zab*(ftaa(i1)*xb*zb+fta(i1)*xab*zb+fta(i1)*xb*zab)
        rfun=ftbb(i1)*za*za*(xab*zb-zab*xb)+sljacdet*ftbbb(i1)*za*za+2.0d0*ftbb(i1)*za*zab*(sljacdet+xa*xb)+2.0d0*ftb(i1)*zab*(zab*xa*xb+za*xab) &
             -2.0d0*za*zb*xab*(ftbb(i1)*za+2.0d0*ftb(i1)*zab)
        sfun=zb*((xab*zb-zab*xb)*(ftaa(i1)*zb-2.0d0*ftab(i1)*za)+sljacdet*(ftaab(i1)*zb-2.0d0*ftabb(i1)*za-2.0d0*ftab(i1)*zab) &
             +(ftab(i1)*za+fta(i1)*zab)*(2.0d0*xab*zb-2.0d0*xb*zab))
        
        aafun=zb*(sljacdet*(pfun+qfun)-3.0d0*xa*zab*(afun+bfun))
        bbfun=za*(sljacdet*(rfun+sfun)+3.0d0*xb*zab*(afun+bfun)-6.0d0*xab*zb*(afun+bfun)) 
        ccfun=-1.0d0*xb*(sljacdet*(pfun+qfun)+3.0d0*xab*za*(afun+bfun))
        ddfun=xa*(sljacdet*(rfun+sfun)-3.0d0*xab*zb*(afun+bfun)+6.0d0*xb*zab*(afun+bfun))
        
        ftxxx(i1)=(aafun-bbfun)/(sljacdet**5)
        ftxxz(i1)=(ccfun+ddfun)/(sljacdet**5)
    enddo
     
     fs(:)=ft(:)
     fsx(:)=ftx(:)
     fsz(:)=ftz(:)
     
     
     DEALLOCATE(fta,ftb,ftaa,ftbb,ftab,ftaaa,ftaab,ftaba,ftabb,ftbba,ftbbb)
     return
    end subroutine
    
    
!-----------------------------------------------------------------
!     FUNCTION slexpfun1 for evaluating expansion function for SLE 
!-----------------------------------------------------------------
   	subroutine slexpfun1(alpha,beta)
    use var_inputdata
    use var_analysis
        implicit none
        real*8::alpha,beta
        
          ft(1)=0.25d0*(1.0d0-alpha)*(1.0d0-beta)   
          ft(2)=0.25d0*(1.0d0+alpha)*(1.0d0-beta)
          ft(3)=0.25d0*(1.0d0+alpha)*(1.0d0+beta)
          ft(4)=0.25d0*(1.0d0-alpha)*(1.0d0+beta) 
          if (mexp.eq.4) then
              GOTO 2001
          endif
          
          ft(5)=0.5d0*(1.0d0-beta)*((alpha**2)-1.0d0)
          ft(6)=0.5d0*(1.0d0+alpha)*((beta**2)-1.0d0)
          ft(7)=0.5d0*(1.0d0+beta)*((alpha**2)-1.0d0)
          ft(8)=0.5d0*(1.0d0-alpha)*((beta**2)-1.0d0)
          if (mexp.eq.8) then
              GOTO 2001
          endif
              
          ft(9)=0.5d0*(1.0d0-beta)*((alpha**3)-alpha)
          ft(10)=0.5d0*(1.0d0+alpha)*((beta**3)-beta)
          !ft(11)=0.5d0*(1.0d0+beta)*((alpha**3)-alpha)
          !ft(12)=0.5d0*(1.0d0-alpha)*((beta**3)-beta)
          ft(11)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**3)+alpha)
          ft(12)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**3)+beta)
          if (mexp.eq.12) then
              GOTO 2001
          endif
          
          ft(13)=0.5d0*(1.0d0-beta)*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ft(14)=0.5d0*(1.0d0+alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(15)=0.5d0*(1.0d0+beta)*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ft(16)=0.5d0*(1.0d0-alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(17)=((alpha**2)-1.0d0)*((beta**2)-1.0d0) 
          if (mexp.eq.17) then
              GOTO 2001
          endif
          
          ft(18)=0.5*(1.0d0-beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+(1.0d0/4.0d0)*alpha)
          ft(19)=0.5*(1.0d0+alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          !ft(20)=0.5*(1.0+beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+(1.0/4.0)*alpha)
          !ft(21)=0.5*(1.0d0-alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(20)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**5)+(5.0d0/4.0d0)*(alpha**3)-(1.0d0/4.0d0)*alpha)
          ft(21)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**5)+(5.0d0/4.0d0)*(beta**3)-(1.0d0/4.0d0)*beta)
          ft(22)=((alpha**2)-1.0d0)*((beta**3)-beta)
          ft(23)=((beta**2)-1.0d0)*((alpha**3)-alpha)  
          if (mexp.eq.23) then
              GOTO 2001
          endif
          
          ft(24)=0.5*(1.0d0-beta)*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ft(25)=0.5*(1.0d0+alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(26)=0.5*(1.0d0+beta)*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ft(27)=0.5*(1.0d0-alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(28)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**2)-1.0d0)
          ft(29)=((alpha**3)-alpha)*((beta**3)-beta)
          ft(30)=((alpha**2)-1.0d0)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0)) 
          if (mexp.eq.30) then
              GOTO 2001
          endif
          
          ft(31)=0.5d0*(1.0d0-beta)*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          ft(32)=0.5d0*(1.0d0+alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          !ft(33)=0.5d0*(1.0d0+beta)*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          !ft(34)=0.5d0*(1.0d0-alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(33)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**7)+(14.0d0/9.0d0)*(alpha**5)-(49.0d0/81.0d0)*(alpha**3)+(4.0d0/81.0d0)*alpha)
          ft(34)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**7)+(14.0d0/9.0d0)*(beta**5)-(49.0d0/81.0d0)*(beta**3)+(4.0d0/81.0d0)*beta) 
          ft(35)=((alpha**2)-1.0d0)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+0.25*beta)
          ft(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**3)-beta) 
          ft(37)=((alpha**3)-alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(38)=((beta**2)-1.0d0)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          if (mexp.eq.38) then
              GOTO 2001
          endif
          
          ft(39)=0.5*(1.0d0-beta)*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ft(40)=0.5*(1.0d0+alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(41)=0.5*(1.0d0+beta)*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ft(42)=0.5*(1.0d0-alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(43)=((alpha**2)-1.0d0)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(44)=((alpha**3)-alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(45)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(46)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**3)-beta)
          ft(47)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**2)-1.0d0)
          if (mexp.eq.47) then
              GOTO 2001
          endif
          
          ft(48)=0.5d0*(1.0d0-beta)*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          ft(49)=0.5d0*(1.0d0+alpha)*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          !ft(50)=0.5d0*(1.0d0+beta)*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          !ft(51)=0.5d0*(1.0d0-alpha)*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          ft(50)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**9)+(15.0d0/8.0d0)*(alpha**7)-(273.0d0/256.0d0)*(alpha**5)+(205.0d0/1024.0d0)*(alpha**3)-(9.0d0/1024.0d0)*alpha)
          ft(51)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**9)+(15.0d0/8.0d0)*(beta**7)-(273.0d0/256.0d0)*(beta**5)+(205.0d0/1024.0d0)*(beta**3)-(9.0d0/1024.0d0)*beta)
          ft(52)=((alpha**2)-1.0d0)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(53)=((alpha**3)-alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(54)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(55)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(56)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**3)-beta)
          ft(57)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*((beta**2)-1.0d0)
          if (mexp.eq.57) then
              GOTO 2001
          endif
          
          ft(58)=0.5*(1.0d0-beta)*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ft(59)=0.5*(1.0d0+alpha)*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          ft(60)=0.5*(1.0d0+beta)*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ft(61)=0.5*(1.0d0-alpha)*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          ft(62)=((alpha**2)-1.0d0)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(63)=((alpha**3)-alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(64)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(65)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(66)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(67)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*((beta**3)-beta)
          ft(68)=((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))*((beta**2)-1.0d0)
          if (mexp.eq.68) then
              GOTO 2001
          endif
     
2001     return
    end subroutine
    
!-----------------------------------------------------------------
!     FUNCTION slexpfun_ab for evaluating expansion function for SLE 
!-----------------------------------------------------------------
   	subroutine slexpfun_ab(ele,cs,alpha,beta,xa,xb,za,zb)
    use var_inputdata
    use var_analysis
        implicit none
        integer::ele,cs
        real*8::alpha,beta
        real*8::xa,xb,za,zb
        real*8,ALLOCATABLE, DIMENSION (:)::fta,ftb
        
        allocate(fta(mexp),ftb(mexp))
        
        xa=0.0d0
        xb=0.0d0
        za=0.0d0
        zb=0.0d0
        
          ft(1)=0.25d0*(1.0d0-alpha)*(1.0d0-beta)   
          ft(2)=0.25d0*(1.0d0+alpha)*(1.0d0-beta)
          ft(3)=0.25d0*(1.0d0+alpha)*(1.0d0+beta)
          ft(4)=0.25d0*(1.0d0-alpha)*(1.0d0+beta) 
          
          fta(1)=-0.25d0*(1.0d0-beta)   
          fta(2)=0.25d0*(1.0d0-beta)
          fta(3)=0.25d0*(1.0d0+beta)
          fta(4)=-0.25d0*(1.0d0+beta)
          
          ftb(1)=-0.25d0*(1.0d0-alpha)  
          ftb(2)=-0.25d0*(1.0d0+alpha)
          ftb(3)=0.25d0*(1.0d0+alpha)
          ftb(4)=0.25d0*(1.0d0-alpha)
          
          if (mexp.eq.4) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
          
          ft(5)=0.5d0*(1.0d0-beta)*((alpha**2)-1.0d0)
          ft(6)=0.5d0*(1.0d0+alpha)*((beta**2)-1.0d0)
          ft(7)=0.5d0*(1.0d0+beta)*((alpha**2)-1.0d0)
          ft(8)=0.5d0*(1.0d0-alpha)*((beta**2)-1.0d0)
          
          fta(5)=0.5d0*(1.0d0-beta)*(2.0*alpha)
          fta(6)=0.5d0*((beta**2)-1.0d0)
          fta(7)=0.5d0*(1.0d0+beta)*(2.0*alpha)
          fta(8)=-0.5d0*((beta**2)-1.0d0)
          
          ftb(5)=-0.5d0*((alpha**2)-1.0d0)
          ftb(6)=0.5d0*(1.0d0+alpha)*(2.0*beta)
          ftb(7)=0.5d0*((alpha**2)-1.0d0)
          ftb(8)=0.5d0*(1.0d0-alpha)*(2.0d0*beta)
          
          if (mexp.eq.8) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(9)=0.5d0*(1.0d0-beta)*((alpha**3)-alpha)
          ft(10)=0.5d0*(1.0d0+alpha)*((beta**3)-beta)
          !ft(11)=0.5d0*(1.0d0+beta)*((alpha**3)-alpha)
          !ft(12)=0.5d0*(1.0d0-alpha)*((beta**3)-beta)
          ft(11)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**3)+alpha)
          ft(12)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**3)+beta)
          
          fta(9)=0.5d0*(1.0d0-beta)*(3.0d0*(alpha**2)-1.0d0)
          fta(10)=0.5d0*((beta**3)-beta)
          !fta(11)=0.5d0*(1.0d0+beta)*(3.0d0*(alpha**2)-1.0d0)
          !fta(12)=-0.5d0*((beta**3)-beta)
          fta(11)=0.5d0*(1.0d0+beta)*(-3.0d0*(alpha**2)+1.0d0)
          fta(12)=-0.5d0*(-1.0d0*(beta**3)+beta)
          
          ftb(9)=-0.5d0*((alpha**3)-alpha)
          ftb(10)=0.5d0*(1.0d0+alpha)*(3.0d0*(beta**2)-1.0d0)
          !ftb(11)=0.5d0*((alpha**3)-alpha)
          !ftb(12)=0.5d0*(1.0d0-alpha)*(3.0d0*(beta**2)-1.0d0)
          ftb(11)=0.5d0*(-1.0d0*(alpha**3)+alpha)
          ftb(12)=0.5d0*(1.0d0-alpha)*(-3.0d0*(beta**2)+1.0d0)
          
          if (mexp.eq.12) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(13)=0.5d0*(1.0d0-beta)*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ft(14)=0.5d0*(1.0d0+alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(15)=0.5d0*(1.0d0+beta)*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ft(16)=0.5d0*(1.0d0-alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(17)=((alpha**2)-1.0d0)*((beta**2)-1.0d0) 
          
          fta(13)=0.5d0*(1.0d0-beta)*(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)
          fta(14)=0.5d0*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(15)=0.5d0*(1.0d0+beta)*(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)
          fta(16)=-0.5d0*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(17)=(2.0*alpha)*((beta**2)-1.0d0)
          
          ftb(13)=-0.5d0*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ftb(14)=0.5d0*(1.0d0+alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta)
          ftb(15)=0.5d0*((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))
          ftb(16)=0.5d0*(1.0d0-alpha)*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(17)=(2.0*beta)*((alpha**2)-1.0d0)
          
          if (mexp.eq.17) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(18)=0.5d0*(1.0d0-beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+(1.0d0/4.0d0)*alpha)
          ft(19)=0.5d0*(1.0d0+alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          !ft(20)=0.5d0*(1.0+beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+(1.0d0/4.0)*alpha)
          !ft(21)=0.5d0*(1.0d0-alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(20)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**5)+(5.0d0/4.0d0)*(alpha**3)-(1.0d0/4.0d0)*alpha)
          ft(21)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**5)+(5.0d0/4.0d0)*(beta**3)-(1.0d0/4.0d0)*beta)
          ft(22)=((alpha**2)-1.0d0)*((beta**3)-beta)
          ft(23)=((beta**2)-1.0d0)*((alpha**3)-alpha)  
          
          fta(18)=0.5d0*(1.0d0-beta)*(5.0d0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25d0)
          fta(19)=0.5d0*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(0.25d0*beta))
          !fta(20)=0.5d0*(1.0d0+beta)*(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)
          !fta(21)=-0.5d0*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(20)=0.5d0*(1.0d0+beta)*(-5.0d0*(alpha**4)+(15.0d0/4.0d0)*(alpha**2)-0.25d0)
          fta(21)=-0.5d0*(-1.0d0*(beta**5)+(5.0d0/4.0d0)*(beta**3)-(1.0d0/4.0d0)*beta)
          fta(22)=(2.0*alpha)*((beta**3)-beta)
          fta(23)=((beta**2)-1.0d0)*(3.0*(alpha**2)-1.0d0)
          
          ftb(18)=-0.5d0*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25d0*alpha)
          ftb(19)=0.5d0*(1.0d0+alpha)*(5.0d0*(beta**4)-(15.0/4.0)*(beta**2)+0.25d0)
          !ftb(20)=0.5d0*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          !ftb(21)=0.5d0*(1.0d0-alpha)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25d0)
          ftb(20)=0.5d0*(-1.0d0*(alpha**5)+(5.0d0/4.0d0)*(alpha**3)-(1.0d0/4.0d0)*alpha)
          ftb(21)=0.5d0*(1.0d0-alpha)*(-5.0d0*(beta**4)+(15.0d0/4.0d0)*(beta**2)-0.25d0)
          ftb(22)=((alpha**2)-1.0d0)*(3.0*(beta**2)-1.0d0)
          ftb(23)=(2.0*beta)*((alpha**3)-alpha)
          
          if (mexp.eq.23) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(24)=0.5d0*(1.0d0-beta)*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ft(25)=0.5d0*(1.0d0+alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(26)=0.5d0*(1.0d0+beta)*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ft(27)=0.5d0*(1.0d0-alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(28)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**2)-1.0d0)
          ft(29)=((alpha**3)-alpha)*((beta**3)-beta)
          ft(30)=((alpha**2)-1.0d0)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0)) 
          
          fta(24)=0.5d0*(1.0d0-beta)*(6.0*(alpha**5)-(28.0d0/5.0d0)*(alpha**3)+(518.0d0/625.0d0)*alpha)
          fta(25)=0.5d0*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(26)=0.5d0*(1.0d0+beta)*(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)
          fta(27)=-0.5d0*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(28)=(4.0d0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**2)-1.0d0)
          fta(29)=(3.0d0*(alpha**2)-1.0d0)*((beta**3)-beta)
          fta(30)=(2.0d0*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          
          ftb(24)=-0.5d0*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ftb(25)=0.5d0*(1.0d0+alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(26)=0.5d0*((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))
          ftb(27)=0.5d0*(1.0d0-alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(28)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(2.0*beta)
          ftb(29)=((alpha**3)-alpha)*(3.0*(beta**2)-1.0d0)
          ftb(30)=((alpha**2)-1.0d0)*(4.0*(beta**3)-(20.0/9.0)*beta)
          
          if (mexp.eq.30) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(31)=0.5d0*(1.0d0-beta)*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          ft(32)=0.5d0*(1.0d0+alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          !ft(33)=0.5d0*(1.0d0+beta)*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          !ft(34)=0.5d0*(1.0d0-alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(33)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**7)+(14.0d0/9.0d0)*(alpha**5)-(49.0d0/81.0d0)*(alpha**3)+(4.0d0/81.0d0)*alpha)
          ft(34)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**7)+(14.0d0/9.0d0)*(beta**5)-(49.0d0/81.0d0)*(beta**3)+(4.0d0/81.0d0)*beta) 
          ft(35)=((alpha**2)-1.0d0)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+0.25*beta)
          ft(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**3)-beta) 
          ft(37)=((alpha**3)-alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(38)=((beta**2)-1.0d0)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          
          fta(31)=0.5d0*(1.0d0-beta)*(7.0d0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))
          fta(32)=0.5d0*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          !fta(33)=0.5d0*(1.0+beta)*(7.0*(alpha**6)-(70.0/9.0)*(alpha**4)+(147.0/81.0d0)*(alpha**2)-(4.0/81.0d0))
          !fta(34)=-0.5d0*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          fta(33)=0.5d0*(1.0d0+beta)*(-7.0d0*(alpha**6)+(70.0/9.0)*(alpha**4)-(147.0/81.0d0)*(alpha**2)+(4.0/81.0d0))
          fta(34)=-0.5d0*(-1.0d0*(beta**7)+(14.0d0/9.0d0)*(beta**5)-(49.0d0/81.0d0)*(beta**3)+(4.0d0/81.0d0)*beta)
          
          fta(35)=(2.0*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+0.25*beta)
          fta(36)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**3)-beta)
          fta(37)=(3.0*(alpha**2)-1.0d0)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(38)=((beta**2)-1.0d0)*(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)
          
          ftb(31)=-0.5d0*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          ftb(32)=0.5d0*(1.0d0+alpha)*(7.0d0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0))   
          !ftb(33)=0.5d0*((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)
          !ftb(34)=0.5d0*(1.0-alpha)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftb(33)=0.5d0*(-1.0d0*(alpha**7)+(14.0d0/9.0d0)*(alpha**5)-(49.0d0/81.0d0)*(alpha**3)+(4.0d0/81.0d0)*alpha)
          ftb(34)=0.5d0*(1.0d0-alpha)*(-7.0d0*(beta**6)+(70.0d0/9.0d0)*(beta**4)-(147.0d0/81.0d0)*(beta**2)+(4.0d0/81.0d0)) 
          
          ftb(35)=((alpha**2)-1.0d0)*(5.0*(beta**4)-(15.0d0/4.0d0)*(beta**2)+0.25)  
          ftb(36)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(3.0*(beta**2)-1.0d0) 
          ftb(37)=((alpha**3)-alpha)*(4.0*(beta**3)-(20.0d0/9.0d0)*beta)   
          ftb(38)=(2.0*beta)*((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)
          
          if (mexp.eq.38) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(39)=0.5d0*(1.0d0-beta)*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ft(40)=0.5d0*(1.0d0+alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(41)=0.5d0*(1.0d0+beta)*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ft(42)=0.5d0*(1.0d0-alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(43)=((alpha**2)-1.0d0)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(44)=((alpha**3)-alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(45)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(46)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**3)-beta)
          ft(47)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**2)-1.0d0)
          
          fta(39)=0.5d0*(1.0d0-beta)*(8.0d0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))
          fta(40)=0.5d0*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          fta(41)=0.5d0*(1.0d0+beta)*(8.0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))
          fta(42)=-0.5d0*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          fta(43)=(2.0*alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(44)=(3.0*(alpha**2)-1.0d0)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(45)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(46)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*((beta**3)-beta)
          fta(47)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*((beta**2)-1.0d0)
          
          ftb(39)=-0.5d0*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ftb(40)=0.5d0*(1.0d0+alpha)*(8.0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftb(41)=0.5d0*((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))
          ftb(42)=0.5d0*(1.0d0-alpha)*(8.0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftb(43)=((alpha**2)-1.0d0)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(44)=((alpha**3)-alpha)*(5.0*(beta**4)-(15.0/4.0)*(beta**2)+0.25)
          ftb(45)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(46)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(3.0*(beta**2)-1.0d0)
          ftb(47)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(2.0*beta)
          
          if (mexp.eq.47) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          
          ft(48)=0.5d0*(1.0d0-beta)*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          ft(49)=0.5d0*(1.0d0+alpha)*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          !ft(50)=0.5d0*(1.0d0+beta)*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          !ft(51)=0.5d0*(1.0d0-alpha)*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          ft(50)=0.5d0*(1.0d0+beta)*(-1.0d0*(alpha**9)+(15.0d0/8.0d0)*(alpha**7)-(273.0d0/256.0d0)*(alpha**5)+(205.0d0/1024.0d0)*(alpha**3)-(9.0d0/1024.0d0)*alpha)
          ft(51)=0.5d0*(1.0d0-alpha)*(-1.0d0*(beta**9)+(15.0d0/8.0d0)*(beta**7)-(273.0d0/256.0d0)*(beta**5)+(205.0d0/1024.0d0)*(beta**3)-(9.0d0/1024.0d0)*beta)
          ft(52)=((alpha**2)-1.0d0)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(53)=((alpha**3)-alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(54)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(55)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(56)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**3)-beta)
          ft(57)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*((beta**2)-1.0d0)
          
          fta(48)=0.5d0*(1.0d0-beta)*(9.0d0*(alpha**8)-(105.0d0/8.0d0)*(alpha**6)+(1365.0d0/256.0d0)*(alpha**4)-(615.0d0/1024.0d0)*(alpha**2)+(9.0d0/1024.0d0))
          fta(49)=0.5d0*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          !fta(50)=0.5d0*(1.0d0+beta)*(9.0d0*(alpha**8)-(105.0d0/8.0d0)*(alpha**6)+(1365.0d0/256.0d0)*(alpha**4)-(615.0d0/1024.0d0)*(alpha**2)+(9.0d0/1024.0d0))
          !fta(51)=-0.5d0*((beta**9)-(15.0d0/8.0d0)*(beta**7)+(273.0d0/256.0d0)*(beta**5)-(205.0d0/1024.0d0)*(beta**3)+(9.0d0/1024.0d0)*beta)
          fta(50)=0.5d0*(1.0d0+beta)*(-9.0d0*(alpha**8)+(105.0d0/8.0d0)*(alpha**6)-(1365.0d0/256.0d0)*(alpha**4)+(615.0d0/1024.0d0)*(alpha**2)-(9.0d0/1024.0d0))
          fta(51)=-0.5d0*(-1.0d0*(beta**9)+(15.0d0/8.0d0)*(beta**7)-(273.0d0/256.0d0)*(beta**5)+(205.0d0/1024.0d0)*(beta**3)-(9.0d0/1024.0d0)*beta)
          fta(52)=(2.0*alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          fta(53)=(3.0*(alpha**2)-1.0d0)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(54)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(55)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(56)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*((beta**3)-beta)
          fta(57)=(7.0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))*((beta**2)-1.0d0)
          
          ftb(48)=-0.5d0*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          ftb(49)=0.5d0*(1.0d0+alpha)*(9.0d0*(beta**8)-(105.0d0/8.0d0)*(beta**6)+(1365.0d0/256.0d0)*(beta**4)-(615.0d0/1024.0d0)*(beta**2)+(9.0d0/1024.0d0))
          !ftb(50)=0.5d0*((alpha**9)-(15.0d0/8.0d0)*(alpha**7)+(273.0d0/256.0d0)*(alpha**5)-(205.0d0/1024.0d0)*(alpha**3)+(9.0d0/1024.0d0)*alpha)
          !ftb(51)=0.5d0*(1.0d0-alpha)*(9.0d0*(beta**8)-(105.0d0/8.0d0)*(beta**6)+(1365.0d0/256.0d0)*(beta**4)-(615.0d0/1024.0d0)*(beta**2)+(9.0d0/1024.0d0))
          ftb(50)=0.5d0*(-1.0d0*(alpha**9)+(15.0d0/8.0d0)*(alpha**7)-(273.0d0/256.0d0)*(alpha**5)+(205.0d0/1024.0d0)*(alpha**3)-(9.0d0/1024.0d0)*alpha)
          ftb(51)=0.5d0*(1.0d0-alpha)*(-9.0d0*(beta**8)+(105.0d0/8.0d0)*(beta**6)-(1365.0d0/256.0d0)*(beta**4)+(615.0d0/1024.0d0)*(beta**2)-(9.0d0/1024.0d0))
          ftb(52)=((alpha**2)-1.0d0)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftb(53)=((alpha**3)-alpha)*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(54)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(5.0*(beta**4)-(15.0/4.0)*(beta**2)+0.25)
          ftb(55)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(56)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(3.0*(beta**2)-1.0d0)
          ftb(57)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*(2.0d0*beta)
          
          if (mexp.eq.57) then
              GOTO 2004
          endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          
          ft(58)=0.5*(1.0d0-beta)*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ft(59)=0.5*(1.0d0+alpha)*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          ft(60)=0.5*(1.0d0+beta)*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ft(61)=0.5*(1.0d0-alpha)*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          ft(62)=((alpha**2)-1.0d0)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          ft(63)=((alpha**3)-alpha)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          ft(64)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          ft(65)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          ft(66)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          ft(67)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*((beta**3)-beta)
          ft(68)=((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))*((beta**2)-1.0d0)
          
          fta(58)=0.5*(1.0d0-beta)*(10.0d0*(alpha**9)-(440.0d0/27.0d0)*(alpha**7)+(5825.0d0/729.0d0)*(alpha**5)-(691240.0d0/531441.0d0)*(alpha**3)+(234938.0d0/4782969.0d0)*alpha)
          fta(59)=0.5*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          fta(60)=0.5*(1.0d0+beta)*(10.0d0*(alpha**9)-(440.0d0/27.0d0)*(alpha**7)+(5825.0d0/729.0d0)*(alpha**5)-(691240.0d0/531441.0d0)*(alpha**3)+(234938.0d0/4782969.0d0)*alpha)
          fta(61)=-0.5*((beta**10)-(55.0d0/27.0d0)*(beta**8)+(2926.0d0/2187.0d0)*(beta**6)-(172810.0d0/531441.0d0)*(beta**4)+(117469.0d0/4782969.0d0)*(beta**2)-(1225.0d0/4782969.0d0))
          fta(62)=(2.0*alpha)*((beta**8)-(12.0d0/7.0d0)*(beta**6)+(282.0d0/343.0d0)*(beta**4)-(12916.0d0/117649.0d0)*(beta**2)+(225.0d0/117649.0d0))
          fta(63)=(3.0*(alpha**2)-1.0d0)*((beta**7)-(14.0d0/9.0d0)*(beta**5)+(49.0d0/81.0d0)*(beta**3)-(4.0d0/81.0d0)*beta)
          fta(64)=(4.0*(alpha**3)-(20.0d0/9.0d0)*alpha)*((beta**6)-(7.0d0/5.0d0)*(beta**4)+(259.0d0/625.0d0)*(beta**2)-(9.0d0/625.0d0))
          fta(65)=(5.0*(alpha**4)-(15.0d0/4.0d0)*(alpha**2)+0.25)*((beta**5)-(5.0d0/4.0d0)*(beta**3)+(1.0d0/4.0d0)*beta)
          fta(66)=(6.0*(alpha**5)-(28.0/5.0)*(alpha**3)+(518.0d0/625.0d0)*alpha)*((beta**4)-(10.0d0/9.0d0)*(beta**2)+(1.0d0/9.0d0))
          fta(67)=(7.0*(alpha**6)-(70.0d0/9.0d0)*(alpha**4)+(147.0d0/81.0d0)*(alpha**2)-(4.0d0/81.0d0))*((beta**3)-beta)
          fta(68)=(8.0d0*(alpha**7)-(72.0d0/7.0d0)*(alpha**5)+(1128.0d0/343.0d0)*(alpha**3)-(25832.0d0/117649.0d0)*(alpha))*((beta**2)-1.0d0)
          
          ftb(58)=-0.5*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ftb(59)=0.5*(1.0d0+alpha)*(10.0d0*(beta**9)-(440.0d0/27.0d0)*(beta**7)+(5825.0d0/729.0d0)*(beta**5)-(691240.0d0/531441.0d0)*(beta**3)+(234938.0d0/4782969.0d0)*beta)
          ftb(60)=0.5*((alpha**10)-(55.0d0/27.0d0)*(alpha**8)+(2926.0d0/2187.0d0)*(alpha**6)-(172810.0d0/531441.0d0)*(alpha**4)+(117469.0d0/4782969.0d0)*(alpha**2)-(1225.0d0/4782969.0d0))
          ftb(61)=0.5*(1.0d0-alpha)*(10.0d0*(beta**9)-(440.0d0/27.0d0)*(beta**7)+(5825.0d0/729.0d0)*(beta**5)-(691240.0d0/531441.0d0)*(beta**3)+(234938.0d0/4782969.0d0)*beta)
          ftb(62)=((alpha**2)-1.0d0)*(8.0*(beta**7)-(72.0d0/7.0d0)*(beta**5)+(1128.0d0/343.0d0)*(beta**3)-(25832.0d0/117649.0d0)*(beta))
          ftb(63)=((alpha**3)-alpha)*(7.0*(beta**6)-(70.0d0/9.0d0)*(beta**4)+(147.0d0/81.0d0)*(beta**2)-(4.0d0/81.0d0)) 
          ftb(64)=((alpha**4)-(10.0d0/9.0d0)*(alpha**2)+(1.0d0/9.0d0))*(6.0*(beta**5)-(28.0d0/5.0d0)*(beta**3)+(518.0d0/625.0d0)*beta)
          ftb(65)=((alpha**5)-(5.0d0/4.0d0)*(alpha**3)+0.25*alpha)*(5.0*(beta**4)-(15.0/4.0)*(beta**2)+0.25)
          ftb(66)=((alpha**6)-(7.0d0/5.0d0)*(alpha**4)+(259.0d0/625.0d0)*(alpha**2)-(9.0d0/625.0d0))*(4.0*(beta**3)-(20.0/9.0)*beta)
          ftb(67)=((alpha**7)-(14.0d0/9.0d0)*(alpha**5)+(49.0d0/81.0d0)*(alpha**3)-(4.0d0/81.0d0)*alpha)*(3.0*(beta**2)-1.0d0)
          ftb(68)=((alpha**8)-(12.0d0/7.0d0)*(alpha**6)+(282.0d0/343.0d0)*(alpha**4)-(12916.0d0/117649.0d0)*(alpha**2)+(225.0d0/117649.0d0))*(2.0*beta)
          
          if (mexp.eq.68) then
              GOTO 2004
          endif
    
2004    xa=fta(1)*cscord(elemcon(ele,2,cs),2,cs)+fta(2)*cscord(elemcon(ele,3,cs),2,cs)+fta(3)*cscord(elemcon(ele,4,cs),2,cs)+fta(4)*cscord(elemcon(ele,5,cs),2,cs)
        xb=ftb(1)*cscord(elemcon(ele,2,cs),2,cs)+ftb(2)*cscord(elemcon(ele,3,cs),2,cs)+ftb(3)*cscord(elemcon(ele,4,cs),2,cs)+ftb(4)*cscord(elemcon(ele,5,cs),2,cs)
        za=fta(1)*cscord(elemcon(ele,2,cs),3,cs)+fta(2)*cscord(elemcon(ele,3,cs),3,cs)+fta(3)*cscord(elemcon(ele,4,cs),3,cs)+fta(4)*cscord(elemcon(ele,5,cs),3,cs)
        zb=ftb(1)*cscord(elemcon(ele,2,cs),3,cs)+ftb(2)*cscord(elemcon(ele,3,cs),3,cs)+ftb(3)*cscord(elemcon(ele,4,cs),3,cs)+ftb(4)*cscord(elemcon(ele,5,cs),3,cs)
        
        deallocate(fta,ftb)
     
     return
    end subroutine 

    
end module
