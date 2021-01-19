module mod_gauss
contains

    subroutine gauss_cs(a)
    use var_analysis
    
        implicit none
        integer::a
        
	        if (a==2)then
     	        gpt_cs(1)=-1.0/sqrt(3.0)
     	        gpt_cs(2)=1.0/sqrt(3.0)
     	        wpt_cs(1)=1.0
     	        wpt_cs(2)=1.0
 	        
 	        elseif (a==3)THEN
     	        gpt_cs(1)=-0.774596669241483
     	        gpt_cs(2)=0.0
     	        gpt_cs(3)=0.774596669241483
     	        wpt_cs(1)=0.555555555555555
     	        wpt_cs(2)=0.888888888888888
                wpt_cs(3)=0.555555555555555
               
            elseif (a==4)THEN
     	        gpt_cs(1)=-0.861136311594053
     	        gpt_cs(2)=-0.339981043584856
     	        gpt_cs(3)=0.339981043584856
                gpt_cs(4)=0.861136311594053
     	        wpt_cs(1)=0.347854845137454
     	        wpt_cs(2)=0.652145154862546
     	        wpt_cs(3)=0.652145154862546
                wpt_cs(4)=0.347854845137454
        
            elseif (a==5)THEN
                gpt_cs(1)=-0.906179845939
                gpt_cs(2)=-0.538469310106
		        gpt_cs(3)=0.0
		        gpt_cs(4)=0.538469310106
		        gpt_cs(5)=0.906179845939
     	        wpt_cs(1)=0.236926885056
		        wpt_cs(2)=0.478628670499
		        wpt_cs(3)=0.568888888889
		        wpt_cs(4)=0.478628670499
		        wpt_cs(5)=0.236926885056
                
          elseif (a==6)THEN
                gpt_cs(1)=-0.9324695142031520
                gpt_cs(2)=-0.6612093864662645
		        gpt_cs(3)=-0.2386191860831969
		        gpt_cs(4)=0.2386191860831969
		        gpt_cs(5)=0.6612093864662645
                gpt_cs(6)=0.9324695142031520
                
     	        wpt_cs(1)=0.1713244923791703
		        wpt_cs(2)=0.3607615730481386
		        wpt_cs(3)=0.4679139345726910
		        wpt_cs(4)=0.4679139345726910
		        wpt_cs(5)=0.3607615730481386
                wpt_cs(6)=0.1713244923791703
                
            elseif (a==7)THEN
                gpt_cs(1)=-0.9491079123427585
                gpt_cs(2)=-0.7415311855993944
		        gpt_cs(3)=-0.4058451513773971
		        gpt_cs(4)=0.0
		        gpt_cs(5)=0.4058451513773971
                gpt_cs(6)=0.7415311855993944
		        gpt_cs(7)=0.9491079123427585
                
     	        wpt_cs(1)=0.1294849661688696
		        wpt_cs(2)=0.2797053914892766
		        wpt_cs(3)=0.3818300505051189
		        wpt_cs(4)=0.4179591836734693
		        wpt_cs(5)=0.3818300505051189
                wpt_cs(6)=0.2797053914892766
		        wpt_cs(7)=0.1294849661688696
                
            elseif (a==8)THEN
                gpt_cs(1)=-0.9602898564975362
                gpt_cs(2)=-0.7966664774136267
		        gpt_cs(3)=-0.5255324099163289
		        gpt_cs(4)=-0.1834346424956498
		        gpt_cs(5)=0.1834346424956498
                gpt_cs(6)=0.5255324099163289
		        gpt_cs(7)=0.7966664774136267
                gpt_cs(8)=0.9602898564975362
                
     	        wpt_cs(1)=0.1012285362903762
		        wpt_cs(2)=0.2223810344533744
		        wpt_cs(3)=0.3137066458778872
		        wpt_cs(4)=0.3626837833783619
		        wpt_cs(5)=0.3626837833783619
                wpt_cs(6)=0.3137066458778872
		        wpt_cs(7)=0.2223810344533744
                wpt_cs(8)=0.1012285362903762
                
            elseif (a==9)THEN
                gpt_cs(1)=-0.9681602395076260
                gpt_cs(2)=-0.8360311073266357
		        gpt_cs(3)=-0.6133714327005903
		        gpt_cs(4)=-0.3242534234038089
		        gpt_cs(5)=0.0000000000000000
                gpt_cs(6)=0.3242534234038089
		        gpt_cs(7)=0.6133714327005903
                gpt_cs(8)=0.8360311073266357
		        gpt_cs(9)=0.9681602395076260
                
     	        wpt_cs(1)=0.0812743883615744
		        wpt_cs(2)=0.1806481606948574
		        wpt_cs(3)=0.2606106964029354
		        wpt_cs(4)=0.3123470770400028
		        wpt_cs(5)=0.3302393550012597
                wpt_cs(6)=0.3123470770400028
		        wpt_cs(7)=0.2606106964029354
                wpt_cs(8)=0.1806481606948574
		        wpt_cs(9)=0.0812743883615744
                
            elseif (a==10)THEN
                gpt_cs(1)=-0.9739065285171717
                gpt_cs(2)=-0.8650633666889845
		        gpt_cs(3)=-0.6794095682990244
		        gpt_cs(4)=-0.4333953941292471
		        gpt_cs(5)=-0.1488743389816312
                gpt_cs(6)=0.1488743389816312
		        gpt_cs(7)=0.4333953941292471
                gpt_cs(8)=0.6794095682990244
		        gpt_cs(9)=0.8650633666889845
                gpt_cs(10)=0.9739065285171717
                
     	        wpt_cs(1)=0.0666713443086881
		        wpt_cs(2)=0.1494513491505805
		        wpt_cs(3)=0.2190863625159820
		        wpt_cs(4)=0.2692667193099963
		        wpt_cs(5)=0.2955242247147528
                wpt_cs(6)=0.2955242247147528
		        wpt_cs(7)=0.2692667193099963
                wpt_cs(8)=0.2190863625159820
		        wpt_cs(9)=0.1494513491505805
                wpt_cs(10)=0.0666713443086881
            endif
                
	return
    end subroutine   
    
    subroutine gauss_cs_tri(a)
    use var_analysis
    
    implicit none
    integer::a,jj
    real*8::n1,n2,n3
    
    !!!!! gauss points taken from Higher-Order Finite Element Methods - Solin et. al. Chapter 4
    if (a==1) then
        agpt_cs_t_test(1)=-1.0d0/3.0d0
        
        bgpt_cs_t_test(1)=-1.0d0/3.0d0
        
        wpt_cs_t_test(1)=2.0d0
        
    elseif (a==3)then
     	agpt_cs_t_test(1)=-2.0d0/3.0d0
        agpt_cs_t_test(2)=-2.0d0/3.0d0
        agpt_cs_t_test(3)=1.0d0/3.0d0
                
        bgpt_cs_t_test(1)=-2.0d0/3.0d0
        bgpt_cs_t_test(2)=1.0d0/3.0d0
        bgpt_cs_t_test(3)=-2.0d0/3.0d0
                
        wpt_cs_t_test(1)=2.0d0/3.0d0
        wpt_cs_t_test(2)=2.0d0/3.0d0
        wpt_cs_t_test(3)=2.0d0/3.0d0
        
    elseif (a==4)then
     	agpt_cs_t_test(1)=-1.0d0/3.0d0
        agpt_cs_t_test(2)=-3.0d0/5.0d0
        agpt_cs_t_test(3)=-3.0d0/5.0d0
        agpt_cs_t_test(4)=1.0d0/5.0d0
                
        bgpt_cs_t_test(1)=-1.0d0/3.0d0
        bgpt_cs_t_test(2)=-3.0d0/5.0d0
        bgpt_cs_t_test(3)=1.0d0/5.0d0
        bgpt_cs_t_test(4)=-3.0d0/5.0d0
                
        wpt_cs_t_test(1)=-1.1250d0
        wpt_cs_t_test(2)=1.041666666666667
        wpt_cs_t_test(3)=1.041666666666667
        wpt_cs_t_test(4)=1.041666666666667
        
    elseif (a==6)then
     	agpt_cs_t_test(1)=-0.108103018168070
        agpt_cs_t_test(2)=-0.108103018168070
        agpt_cs_t_test(3)=-0.783793963663860
        agpt_cs_t_test(4)=-0.816847572980458
        agpt_cs_t_test(5)=-0.816847572980458
        agpt_cs_t_test(6)=0.633695145960918
                
        bgpt_cs_t_test(1)=-0.108103018168070
        bgpt_cs_t_test(2)=-0.783793963663860
        bgpt_cs_t_test(3)=-0.108103018168070
        bgpt_cs_t_test(4)=-0.816847572980458
        bgpt_cs_t_test(5)=0.633695145960918
        bgpt_cs_t_test(6)=-0.816847572980458
                
        wpt_cs_t_test(1)=0.446763179356022
        wpt_cs_t_test(2)=0.446763179356022
        wpt_cs_t_test(3)=0.446763179356022
        wpt_cs_t_test(4)=0.219903487310644
        wpt_cs_t_test(5)=0.219903487310644
        wpt_cs_t_test(6)=0.219903487310644
        
    elseif (a==7)then
     	agpt_cs_t_test(1)=-0.333333333333333
        agpt_cs_t_test(2)=-0.059715871789770
        agpt_cs_t_test(3)=-0.059715871789770
        agpt_cs_t_test(4)=-0.880568256420460
        agpt_cs_t_test(5)=-0.797426985353088
        agpt_cs_t_test(6)=-0.797426985353088
        agpt_cs_t_test(7)=0.594853970706174
                
        bgpt_cs_t_test(1)=-0.333333333333333
        bgpt_cs_t_test(2)=-0.059715871789770
        bgpt_cs_t_test(3)=-0.880568256420460
        bgpt_cs_t_test(4)=-0.059715871789770
        bgpt_cs_t_test(5)=-0.797426985353088
        bgpt_cs_t_test(6)=0.594853970706174
        bgpt_cs_t_test(7)=-0.797426985353088

        wpt_cs_t_test(1)=0.450d0
        wpt_cs_t_test(2)=0.264788305577012
        wpt_cs_t_test(3)=0.264788305577012
        wpt_cs_t_test(4)=0.264788305577012
        wpt_cs_t_test(5)=0.251878361089654
        wpt_cs_t_test(6)=0.251878361089654
        wpt_cs_t_test(7)=0.251878361089654
        
    elseif (a==12)then
     	agpt_cs_t_test(1)=-0.2521397644872685
        agpt_cs_t_test(2)=0.0d0
        agpt_cs_t_test(3)=0.2521397644872685
        agpt_cs_t_test(4)=-0.810732956525494
        agpt_cs_t_test(5)=0.0d0
        agpt_cs_t_test(6)=0.810732956525494
        agpt_cs_t_test(7)=-0.326150048087615
        agpt_cs_t_test(8)=0.326150048087615
        agpt_cs_t_test(9)=-0.583357449276582
        agpt_cs_t_test(10)=0.583357449276582
        agpt_cs_t_test(11)=-0.257207401188967
        agpt_cs_t_test(12)=0.257207401188967
                
        bgpt_cs_t_test(1)=0.249286745170910*sqrt(3.0d0)
        bgpt_cs_t_test(2)=0.501426509658179*sqrt(3.0d0) 
        bgpt_cs_t_test(3)=0.249286745170910*sqrt(3.0d0) 
        bgpt_cs_t_test(4)=0.063089014491502*sqrt(3.0d0) 
        bgpt_cs_t_test(5)=0.873821971016996*sqrt(3.0d0) 
        bgpt_cs_t_test(6)=0.063089014491502*sqrt(3.0d0) 
        bgpt_cs_t_test(7)=0.053145049844817*sqrt(3.0d0) 
        bgpt_cs_t_test(8)=0.053145049844817*sqrt(3.0d0) 
        bgpt_cs_t_test(9)=0.310352451033784*sqrt(3.0d0) 
        bgpt_cs_t_test(10)=0.310352451033784*sqrt(3.0d0) 
        bgpt_cs_t_test(11)=0.636502499121399*sqrt(3.0d0) 
        bgpt_cs_t_test(12)=0.636502499121399*sqrt(3.0d0) 
                
        wpt_cs_t_test(1)=0.116786275726379*sqrt(3.0d0)
        wpt_cs_t_test(2)=0.116786275726379*sqrt(3.0d0)
        wpt_cs_t_test(3)=0.116786275726379*sqrt(3.0d0) 
        wpt_cs_t_test(4)=0.050844906370207*sqrt(3.0d0) 
        wpt_cs_t_test(5)=0.050844906370207*sqrt(3.0d0) 
        wpt_cs_t_test(6)=0.050844906370207*sqrt(3.0d0) 
        wpt_cs_t_test(7)=0.082851075618374*sqrt(3.0d0) 
        wpt_cs_t_test(8)=0.082851075618374*sqrt(3.0d0) 
        wpt_cs_t_test(9)=0.082851075618374*sqrt(3.0d0) 
        wpt_cs_t_test(10)=0.082851075618374*sqrt(3.0d0) 
        wpt_cs_t_test(11)=0.082851075618374*sqrt(3.0d0) 
        wpt_cs_t_test(12)=0.082851075618374*sqrt(3.0d0) 
		
		!agpt_cs_t_test(1)=-0.501426509658180
  !      agpt_cs_t_test(2)=-0.501426509658180
  !      agpt_cs_t_test(3)=0.002853019316358
  !      agpt_cs_t_test(4)=-0.873821971016996
  !      agpt_cs_t_test(5)=-0.873821971016996
  !      agpt_cs_t_test(6)=0.747643942033992
  !      agpt_cs_t_test(7)=-0.379295097932432
  !      agpt_cs_t_test(8)=0.273004998242798
  !      agpt_cs_t_test(9)=-0.893709900310366
  !      agpt_cs_t_test(10)=-0.379295097932432
  !      agpt_cs_t_test(11)=0.273004998242798
  !      agpt_cs_t_test(12)=-0.893709900310366
  !              
  !      bgpt_cs_t_test(1)=-0.501426509658180
  !      bgpt_cs_t_test(2)=0.002853019316358
  !      bgpt_cs_t_test(3)=-0.501426509658180
  !      bgpt_cs_t_test(4)=-0.873821971016996
  !      bgpt_cs_t_test(5)=0.747643942033992
  !      bgpt_cs_t_test(6)=-0.873821971016996
  !      bgpt_cs_t_test(7)=0.273004998242798
  !      bgpt_cs_t_test(8)=-0.893709900310366
  !      bgpt_cs_t_test(9)=-0.379295097932432
  !      bgpt_cs_t_test(10)=-0.893709900310366
  !      bgpt_cs_t_test(11)=-0.379295097932432
  !      bgpt_cs_t_test(12)=0.273004998242798
  !              
  !      wpt_cs_t_test(1)=0.233572551452758
  !      wpt_cs_t_test(2)=0.233572551452758
  !      wpt_cs_t_test(3)=0.233572551452758
  !      wpt_cs_t_test(4)=0.101689812740414
  !      wpt_cs_t_test(5)=0.101689812740414
  !      wpt_cs_t_test(6)=0.101689812740414
  !      wpt_cs_t_test(7)=0.165702151236748
  !      wpt_cs_t_test(8)=0.165702151236748
  !      wpt_cs_t_test(9)=0.165702151236748
  !      wpt_cs_t_test(10)=0.165702151236748
  !      wpt_cs_t_test(11)=0.165702151236748
  !      wpt_cs_t_test(12)=0.165702151236748
        
	elseif (a==13)then
     	!agpt_cs_t_test(1)=-0.333333333333333
      !  agpt_cs_t_test(2)=-0.479308067841920
      !  agpt_cs_t_test(3)=-0.479308067841920
      !  agpt_cs_t_test(4)=-0.041383864316160
      !  agpt_cs_t_test(5)=-0.869739794195568
      !  agpt_cs_t_test(6)=-0.869739794195568
      !  agpt_cs_t_test(7)=0.739479588391136
      !  agpt_cs_t_test(8)=-0.374269007990252
      !  agpt_cs_t_test(9)=0.276888377139620
      !  agpt_cs_t_test(10)=-0.902619369149368
      !  agpt_cs_t_test(11)=-0.374269007990252
      !  agpt_cs_t_test(12)=0.276888377139620
      !  agpt_cs_t_test(13)=-0.902619369149368
      !          
      !  bgpt_cs_t_test(1)=-0.333333333333333
      !  bgpt_cs_t_test(2)=-0.479308067841920
      !  bgpt_cs_t_test(3)=-0.041383864316160
      !  bgpt_cs_t_test(4)=-0.479308067841920
      !  bgpt_cs_t_test(5)=-0.869739794195568
      !  bgpt_cs_t_test(6)=0.739479588391136
      !  bgpt_cs_t_test(7)=-0.869739794195568
      !  bgpt_cs_t_test(8)=0.276888377139620
      !  bgpt_cs_t_test(9)=-0.902619369149368
      !  bgpt_cs_t_test(10)=-0.374269007990252
      !  bgpt_cs_t_test(11)=-0.902619369149368
      !  bgpt_cs_t_test(12)=-0.374269007990252
      !  bgpt_cs_t_test(13)=0.276888377139620
      !          
      !  wpt_cs_t_test(1)=-0.299140088935364
      !  wpt_cs_t_test(2)=0.351230514866416
      !  wpt_cs_t_test(3)=0.351230514866416
      !  wpt_cs_t_test(4)=0.351230514866416
      !  wpt_cs_t_test(5)=0.106694471217676
      !  wpt_cs_t_test(6)=0.106694471217676
      !  wpt_cs_t_test(7)=0.106694471217676
      !  wpt_cs_t_test(8)=0.154227521780514
      !  wpt_cs_t_test(9)=0.154227521780514
      !  wpt_cs_t_test(10)=0.154227521780514
      !  wpt_cs_t_test(11)=0.154227521780514
      !  wpt_cs_t_test(12)=0.154227521780514
      !  wpt_cs_t_test(13)=0.154227521780514
		
	elseif (a==16)then
     	agpt_cs_t_test(1)=0.0d0
        agpt_cs_t_test(2)=-0.377877764878169
        agpt_cs_t_test(3)=0.0d0
        agpt_cs_t_test(4)=0.377877764878169
        agpt_cs_t_test(5)=-0.48829207674472
        agpt_cs_t_test(6)=0.0d0
        agpt_cs_t_test(7)=0.48829207674472
        agpt_cs_t_test(8)=-0.848358315048907
        agpt_cs_t_test(9)=0.0d0
        agpt_cs_t_test(10)=0.848358315048907
        agpt_cs_t_test(11)=-0.465379563320766
        agpt_cs_t_test(12)=0.465379563320766
        agpt_cs_t_test(13)=-0.720097615545446
		agpt_cs_t_test(14)=0.720097615545446
        agpt_cs_t_test(15)=-0.25471805222468
        agpt_cs_t_test(16)=0.25471805222468
                
        bgpt_cs_t_test(1)=0.5773502691896252
        bgpt_cs_t_test(2)=0.7955180984628107
        bgpt_cs_t_test(3)=0.1410146106432558
        bgpt_cs_t_test(4)=0.7955180984628107
        bgpt_cs_t_test(5)=0.2954347072379003
        bgpt_cs_t_test(6)=1.141181393093077
        bgpt_cs_t_test(7)=0.2954347072379003
        bgpt_cs_t_test(8)=0.087550367626882
        bgpt_cs_t_test(9)=1.556950072315113
        bgpt_cs_t_test(10)=0.087550367626882
        bgpt_cs_t_test(11)=0.0145401809922787
        bgpt_cs_t_test(12)=0.0145401809922787
        bgpt_cs_t_test(13)=0.4557247890504072
		bgpt_cs_t_test(14)=0.4557247890504072
        bgpt_cs_t_test(15)=1.261785837526191
        bgpt_cs_t_test(16)=1.261785837526191
                
        wpt_cs_t_test(1)=0.144315607677787*sqrt(3.0d0)
        wpt_cs_t_test(2)=0.095091634267285*sqrt(3.0d0)
        wpt_cs_t_test(3)=0.095091634267285*sqrt(3.0d0)
        wpt_cs_t_test(4)=0.095091634267285*sqrt(3.0d0)
        wpt_cs_t_test(5)=0.103217370534718*sqrt(3.0d0)
        wpt_cs_t_test(6)=0.103217370534718*sqrt(3.0d0)
        wpt_cs_t_test(7)=0.103217370534718*sqrt(3.0d0)
        wpt_cs_t_test(8)=0.032458497623198*sqrt(3.0d0)
        wpt_cs_t_test(9)=0.032458497623198*sqrt(3.0d0)
        wpt_cs_t_test(10)=0.032458497623198*sqrt(3.0d0)
        wpt_cs_t_test(11)=0.027230314174435*sqrt(3.0d0)
        wpt_cs_t_test(12)=0.027230314174435*sqrt(3.0d0)
        wpt_cs_t_test(13)=0.027230314174435*sqrt(3.0d0)
		wpt_cs_t_test(14)=0.027230314174435*sqrt(3.0d0)
        wpt_cs_t_test(15)=0.027230314174435*sqrt(3.0d0)
        wpt_cs_t_test(16)=0.027230314174435*sqrt(3.0d0)
		
	elseif (a==19)then
     	agpt_cs_t_test(1)=0.0d0
        agpt_cs_t_test(2)=-0.469047557596212
        agpt_cs_t_test(3)=0.0d0
        agpt_cs_t_test(4)=0.469047557596212
        agpt_cs_t_test(5)=-0.311268774478811
        agpt_cs_t_test(6)=0.0d0
        agpt_cs_t_test(7)=0.311268774478811
        agpt_cs_t_test(8)=-0.435389393142901
        agpt_cs_t_test(9)=0.0d0
        agpt_cs_t_test(10)=0.435389393142901
        agpt_cs_t_test(11)=-0.865811459816641
        agpt_cs_t_test(12)=0.0d0
        agpt_cs_t_test(13)=0.865811459816641
		agpt_cs_t_test(14)=-0.519235609623732
        agpt_cs_t_test(15)=0.519235609623732
        agpt_cs_t_test(16)=-0.704360186729762
		agpt_cs_t_test(17)=0.704360186729762
        agpt_cs_t_test(18)=-0.18512457710603
        agpt_cs_t_test(19)=0.18512457710603
                
        bgpt_cs_t_test(1)=0.5773502691896252
        bgpt_cs_t_test(2)=0.8481550028305364
        bgpt_cs_t_test(3)=0.0357408019078062
        bgpt_cs_t_test(4)=0.8481550028305364
        bgpt_cs_t_test(5)=0.7570613799252922
        bgpt_cs_t_test(6)=0.2179280477182946
        bgpt_cs_t_test(7)=0.7570613799252922
        bgpt_cs_t_test(8)=0.3259780858562641
        bgpt_cs_t_test(9)=1.080094635856351
        bgpt_cs_t_test(10)=0.3259780858562641
        bgpt_cs_t_test(11)=0.0774737897970252
        bgpt_cs_t_test(12)=1.577103227974829
        bgpt_cs_t_test(13)=0.0774737897970252
		bgpt_cs_t_test(14)=0.0638060013489606
        bgpt_cs_t_test(15)=0.0638060013489606
        bgpt_cs_t_test(16)=0.3844511746263067
		bgpt_cs_t_test(17)=0.3844511746263067
        bgpt_cs_t_test(18)=1.28379363159361
        bgpt_cs_t_test(19)=1.28379363159361
                
        wpt_cs_t_test(1)=0.097135796282799*sqrt(3.0d0)
        wpt_cs_t_test(2)=0.031334700227139*sqrt(3.0d0)
        wpt_cs_t_test(3)=0.031334700227139*sqrt(3.0d0)
        wpt_cs_t_test(4)=0.031334700227139*sqrt(3.0d0)
        wpt_cs_t_test(5)=0.077827541004774*sqrt(3.0d0)
        wpt_cs_t_test(6)=0.077827541004774*sqrt(3.0d0)
        wpt_cs_t_test(7)=0.077827541004774*sqrt(3.0d0)
        wpt_cs_t_test(8)=0.079647738927210*sqrt(3.0d0)
        wpt_cs_t_test(9)=0.079647738927210*sqrt(3.0d0)
        wpt_cs_t_test(10)=0.079647738927210*sqrt(3.0d0)
        wpt_cs_t_test(11)=0.025577675658698*sqrt(3.0d0)
        wpt_cs_t_test(12)=0.025577675658698*sqrt(3.0d0)
        wpt_cs_t_test(13)=0.025577675658698*sqrt(3.0d0)
		wpt_cs_t_test(14)=0.043283539377289*sqrt(3.0d0)
        wpt_cs_t_test(15)=0.043283539377289*sqrt(3.0d0)
        wpt_cs_t_test(16)=0.043283539377289*sqrt(3.0d0)
		wpt_cs_t_test(17)=0.043283539377289*sqrt(3.0d0)
        wpt_cs_t_test(18)=0.043283539377289*sqrt(3.0d0)
        wpt_cs_t_test(19)=0.043283539377289*sqrt(3.0d0)
		
	elseif (a==25)then
     	agpt_cs_t_test(1)=0.0d0
        agpt_cs_t_test(2)=-0.456732900150971
        agpt_cs_t_test(3)=0.0d0
        agpt_cs_t_test(4)=0.456732900150971
        agpt_cs_t_test(5)=-0.671555273544889
        agpt_cs_t_test(6)=0.0d0
        agpt_cs_t_test(7)=0.671555273544889
        agpt_cs_t_test(8)=-0.242413103056878
        agpt_cs_t_test(9)=0.242413103056878
        agpt_cs_t_test(10)=-0.408645722406119
        agpt_cs_t_test(11)=0.408645722406119
        agpt_cs_t_test(12)=-0.166232619349241
        agpt_cs_t_test(13)=0.166232619349241
		agpt_cs_t_test(14)=-0.481651343957508
        agpt_cs_t_test(15)=0.481651343957508
        agpt_cs_t_test(16)=-0.703320369834725
		agpt_cs_t_test(17)=0.703320369834725
        agpt_cs_t_test(18)=-0.221669025877217
        agpt_cs_t_test(19)=0.221669025877217
		agpt_cs_t_test(20)=-0.856852682575301
        agpt_cs_t_test(21)=0.856852682575301
        agpt_cs_t_test(22)=-0.9141151181872
		agpt_cs_t_test(23)=0.9141151181872
        agpt_cs_t_test(24)=-0.0572624356119
        agpt_cs_t_test(25)=0.0572624356119
                
        bgpt_cs_t_test(1)=0.5773502691896252
        bgpt_cs_t_test(2)=0.8410451320395473
        bgpt_cs_t_test(3)=0.0499605434897809
        bgpt_cs_t_test(4)=0.8410451320395473
        bgpt_cs_t_test(5)=0.1896276512327713
        bgpt_cs_t_test(6)=1.352795505103335
        bgpt_cs_t_test(7)=0.1896276512327713
        bgpt_cs_t_test(8)=0.245444103825883
        bgpt_cs_t_test(9)=0.245444103825883
        bgpt_cs_t_test(10)=0.5333674464140256
        bgpt_cs_t_test(11)=0.5333674464140256
        bgpt_cs_t_test(12)=0.9532392573289687
        bgpt_cs_t_test(13)=0.9532392573289687
		bgpt_cs_t_test(14)=0.0433073925777868
        bgpt_cs_t_test(15)=0.0433073925777868
        bgpt_cs_t_test(16)=0.4272494078614268
		bgpt_cs_t_test(17)=0.4272494078614268
        bgpt_cs_t_test(18)=1.261494007129664
        bgpt_cs_t_test(19)=1.261494007129664
		bgpt_cs_t_test(20)=0.0165251770189535
        bgpt_cs_t_test(21)=0.0165251770189535
        bgpt_cs_t_test(22)=0.1157066248639074
		bgpt_cs_t_test(23)=0.1157066248639074
        bgpt_cs_t_test(24)=1.599819005686015
        bgpt_cs_t_test(25)=1.599819005686015
                
        wpt_cs_t_test(1)=0.090817990382754*sqrt(3.0d0)
        wpt_cs_t_test(2)=0.036725957756467*sqrt(3.0d0)
        wpt_cs_t_test(3)=0.036725957756467*sqrt(3.0d0)
        wpt_cs_t_test(4)=0.036725957756467*sqrt(3.0d0)
        wpt_cs_t_test(5)=0.045321059435528*sqrt(3.0d0)
        wpt_cs_t_test(6)=0.045321059435528*sqrt(3.0d0)
        wpt_cs_t_test(7)=0.045321059435528*sqrt(3.0d0)
        wpt_cs_t_test(8)=0.072757916845420*sqrt(3.0d0)
        wpt_cs_t_test(9)=0.072757916845420*sqrt(3.0d0)
        wpt_cs_t_test(10)=0.072757916845420*sqrt(3.0d0)
        wpt_cs_t_test(11)=0.072757916845420*sqrt(3.0d0)
        wpt_cs_t_test(12)=0.072757916845420*sqrt(3.0d0)
        wpt_cs_t_test(13)=0.072757916845420*sqrt(3.0d0)
		wpt_cs_t_test(14)=0.028327242531057*sqrt(3.0d0)
        wpt_cs_t_test(15)=0.028327242531057*sqrt(3.0d0)
        wpt_cs_t_test(16)=0.028327242531057*sqrt(3.0d0)
		wpt_cs_t_test(17)=0.028327242531057*sqrt(3.0d0)
        wpt_cs_t_test(18)=0.028327242531057*sqrt(3.0d0)
        wpt_cs_t_test(19)=0.028327242531057*sqrt(3.0d0)
		wpt_cs_t_test(20)=0.009421666963733*sqrt(3.0d0)
        wpt_cs_t_test(21)=0.009421666963733*sqrt(3.0d0)
        wpt_cs_t_test(22)=0.009421666963733*sqrt(3.0d0)
		wpt_cs_t_test(23)=0.009421666963733*sqrt(3.0d0)
        wpt_cs_t_test(24)=0.009421666963733*sqrt(3.0d0)
        wpt_cs_t_test(25)=0.009421666963733*sqrt(3.0d0)
		
	elseif (a==33)then
     	agpt_cs_t_test(1)=-0.464652169321415
        agpt_cs_t_test(2)=0.0d0
        agpt_cs_t_test(3)=0.464652169321415
        agpt_cs_t_test(4)=-0.319173176883382
        agpt_cs_t_test(5)=0.0d0
        agpt_cs_t_test(6)=0.319173176883382
        agpt_cs_t_test(7)=-0.186368844963652
        agpt_cs_t_test(8)=0.0d0
        agpt_cs_t_test(9)=0.186368844963652
        agpt_cs_t_test(10)=-0.617271563375242
        agpt_cs_t_test(11)=0.0d0
        agpt_cs_t_test(12)=0.617271563375242
        agpt_cs_t_test(13)=-0.93604794864037
		agpt_cs_t_test(14)=0.0d0
        agpt_cs_t_test(15)=0.93604794864037
        agpt_cs_t_test(16)=-0.333229966094282
		agpt_cs_t_test(17)=0.333229966094282
        agpt_cs_t_test(18)=-0.493599741245086
        agpt_cs_t_test(19)=0.493599741245086
		agpt_cs_t_test(20)=-0.160369775150816
        agpt_cs_t_test(21)=0.160369775150816
        agpt_cs_t_test(22)=-0.414510505797863
		agpt_cs_t_test(23)=0.414510505797863
        agpt_cs_t_test(24)=-0.672997754565546
        agpt_cs_t_test(25)=0.672997754565546
		agpt_cs_t_test(26)=-0.258487248767683
		agpt_cs_t_test(27)=0.258487248767683
        agpt_cs_t_test(28)=-0.741762117636476
        agpt_cs_t_test(29)=0.741762117636476
		agpt_cs_t_test(30)=-0.832279982995737
        agpt_cs_t_test(31)=0.832279982995737
        agpt_cs_t_test(32)=-0.090517865359267
		agpt_cs_t_test(33)=0.090517865359267
                
        bgpt_cs_t_test(1)=0.488217389773805*sqrt(3.0d0)
        bgpt_cs_t_test(2)=0.023565220452390*sqrt(3.0d0)
        bgpt_cs_t_test(3)=0.488217389773805*sqrt(3.0d0)
        bgpt_cs_t_test(4)=0.439724392294460*sqrt(3.0d0)
        bgpt_cs_t_test(5)=0.120551215411079*sqrt(3.0d0)
        bgpt_cs_t_test(6)=0.439724392294460*sqrt(3.0d0)
        bgpt_cs_t_test(7)=0.271210385012116*sqrt(3.0d0)
        bgpt_cs_t_test(8)=0.457579229975768*sqrt(3.0d0)
        bgpt_cs_t_test(9)=0.271210385012116*sqrt(3.0d0)
        bgpt_cs_t_test(10)=0.127576145541586*sqrt(3.0d0)
        bgpt_cs_t_test(11)=0.744847708916828*sqrt(3.0d0)
        bgpt_cs_t_test(12)=0.127576145541586*sqrt(3.0d0)
        bgpt_cs_t_test(13)=0.021317350453210*sqrt(3.0d0)
		bgpt_cs_t_test(14)=0.957365299093579*sqrt(3.0d0)
        bgpt_cs_t_test(15)=0.021317350453210*sqrt(3.0d0)
        bgpt_cs_t_test(16)=0.115343494534698*sqrt(3.0d0)
		bgpt_cs_t_test(17)=0.115343494534698*sqrt(3.0d0)
        bgpt_cs_t_test(18)=0.275713269685514*sqrt(3.0d0)
        bgpt_cs_t_test(19)=0.275713269685514*sqrt(3.0d0)
		bgpt_cs_t_test(20)=0.608943235779788*sqrt(3.0d0)
        bgpt_cs_t_test(21)=0.608943235779788*sqrt(3.0d0)
        bgpt_cs_t_test(22)=0.022838332222257*sqrt(3.0d0)
		bgpt_cs_t_test(23)=0.022838332222257*sqrt(3.0d0)
        bgpt_cs_t_test(24)=0.281325580989940*sqrt(3.0d0)
        bgpt_cs_t_test(25)=0.281325580989940*sqrt(3.0d0)
		bgpt_cs_t_test(26)=0.695836086787803*sqrt(3.0d0)
		bgpt_cs_t_test(27)=0.695836086787803*sqrt(3.0d0)
        bgpt_cs_t_test(28)=0.025734050548330*sqrt(3.0d0)
        bgpt_cs_t_test(29)=0.025734050548330*sqrt(3.0d0)
		bgpt_cs_t_test(30)=0.116251915907597*sqrt(3.0d0)
        bgpt_cs_t_test(31)=0.116251915907597*sqrt(3.0d0)
        bgpt_cs_t_test(32)=0.858014033544073*sqrt(3.0d0)
		bgpt_cs_t_test(33)=0.858014033544073*sqrt(3.0d0)
                
        wpt_cs_t_test(1)=0.025731066440455*sqrt(3.0d0)
        wpt_cs_t_test(2)=0.025731066440455*sqrt(3.0d0)
        wpt_cs_t_test(3)=0.025731066440455*sqrt(3.0d0)
        wpt_cs_t_test(4)=0.043692544538038*sqrt(3.0d0)
        wpt_cs_t_test(5)=0.043692544538038*sqrt(3.0d0)
        wpt_cs_t_test(6)=0.043692544538038*sqrt(3.0d0)
        wpt_cs_t_test(7)=0.062858224217885*sqrt(3.0d0)
        wpt_cs_t_test(8)=0.062858224217885*sqrt(3.0d0)
        wpt_cs_t_test(9)=0.062858224217885*sqrt(3.0d0)
        wpt_cs_t_test(10)=0.034796112930709*sqrt(3.0d0)
        wpt_cs_t_test(11)=0.034796112930709*sqrt(3.0d0)
        wpt_cs_t_test(12)=0.034796112930709*sqrt(3.0d0)
        wpt_cs_t_test(13)=0.006166261051559*sqrt(3.0d0)
		wpt_cs_t_test(14)=0.006166261051559*sqrt(3.0d0)
        wpt_cs_t_test(15)=0.006166261051559*sqrt(3.0d0)
        wpt_cs_t_test(16)=0.040371557766381*sqrt(3.0d0)
		wpt_cs_t_test(17)=0.040371557766381*sqrt(3.0d0)
        wpt_cs_t_test(18)=0.040371557766381*sqrt(3.0d0)
        wpt_cs_t_test(19)=0.040371557766381*sqrt(3.0d0)
		wpt_cs_t_test(20)=0.040371557766381*sqrt(3.0d0)
        wpt_cs_t_test(21)=0.040371557766381*sqrt(3.0d0)
        wpt_cs_t_test(22)=0.022356773202303*sqrt(3.0d0)
		wpt_cs_t_test(23)=0.022356773202303*sqrt(3.0d0)
        wpt_cs_t_test(24)=0.022356773202303*sqrt(3.0d0)
        wpt_cs_t_test(25)=0.022356773202303*sqrt(3.0d0)
		wpt_cs_t_test(26)=0.022356773202303*sqrt(3.0d0)
		wpt_cs_t_test(27)=0.022356773202303*sqrt(3.0d0)
        wpt_cs_t_test(28)=0.017316231108659*sqrt(3.0d0)
        wpt_cs_t_test(29)=0.017316231108659*sqrt(3.0d0)
		wpt_cs_t_test(30)=0.017316231108659*sqrt(3.0d0)
        wpt_cs_t_test(31)=0.017316231108659*sqrt(3.0d0)
        wpt_cs_t_test(32)=0.017316231108659*sqrt(3.0d0)
		wpt_cs_t_test(33)=0.017316231108659*sqrt(3.0d0)
    endif
    
 !   do jj=1,a
 !       n1=-0.5d0*(agpt_cs_t_test(jj)+bgpt_cs_t_test(jj))
 !       n2=0.5d0*(1.0d0+agpt_cs_t_test(jj))
 !       n3=0.5d0*(1.0d0+bgpt_cs_t_test(jj))
 !       agpt_cs_t(jj)=-1.0d0*n1+1.0d0*n2+0.0d0*n3
	!	if (abs(agpt_cs_t(jj)).le.1e-8) then
	!		agpt_cs_t(jj)=0.0d0
	!	endif
 !       bgpt_cs_t(jj)=0.0d0*n1+0.0d0*n2+sqrt(3.0d0)*n3
 !       wpt_cs_t(jj)=0.5d0*wpt_cs_t_test(jj)*sqrt(3.0d0)
 !      ! write(*,*)agpt_cs_t(jj),bgpt_cs_t(jj),wpt_cs_t(jj)
	!enddo
	
	do jj=1,a
        agpt_cs_t(jj)=agpt_cs_t_test(jj)
        bgpt_cs_t(jj)=bgpt_cs_t_test(jj)
        wpt_cs_t(jj)=wpt_cs_t_test(jj)
       ! write(*,*)agpt_cs_t(jj),bgpt_cs_t(jj),wpt_cs_t(jj)
    enddo

	return
    end subroutine   
    
 !   subroutine gauss_cs_tri_test(a)
 !   use var_analysis
 !   
 !   implicit none
 !   integer::a
 !   
 !   !!!!!!!!!!!!!!!!!
 !   integer::ii,jj
 !   real*8::n1,n2,n3,alp1,bet1,alp2,bet2,alp3,bet3
 !   real*8,allocatable,dimension(:)::aa,bb,alp,bet
 !   
 !   ii=6
 !   allocate(aa(ii),bb(ii),alp(ii),bet(ii))
 !   aa=0.0d0; bb=0.0d0; alp=0.0d0; bet=0.0d0
 !   
 !   alp1=-1.0d0; bet1=0.0d0
 !   alp2=1.0d0; bet2=0.0d0
 !   alp3=0.0d0; bet3=sqrt(3.0d0)
 !   
 !   aa(1)=-0.108103018168070;    bb(1)=-0.108103018168070
 !   aa(2)=-0.108103018168070;    bb(2)=-0.783793963663860
 !   aa(3)=-0.783793963663860;    bb(3)=-0.108103018168070
 !   
 !   aa(4)=-0.816847572980458;    bb(4)=-0.816847572980458
 !   aa(5)=-0.816847572980458;    bb(5)=0.633695145960918
 !   aa(6)=0.633695145960918;    bb(6)=-0.816847572980458
 !   
 !   do jj=1,ii
 !   n1=-0.5d0*(aa(jj)+bb(jj))
 !   n2=0.5d0*(1.0d0+aa(jj))
 !   n3=0.5d0*(1.0d0+bb(jj))
 !   alp(jj)=n1*alp1+n2*alp2+n3*alp3
 !   bet(jj)=n1*bet1+n2*bet2+n3*bet3
 !   write(*,*)alp(jj),bet(jj)
 !   enddo
 !
 !       !!!!!!!!!!!!!!!!!!
 !       
 !       
	!if (a==3)then
 !    	agpt_cs_t(1)=-0.5d0
 !       agpt_cs_t(2)=0.0d0
 !       agpt_cs_t(3)=0.5d0
 !               
 !       bgpt_cs_t(1)=sqrt(3.0d0)/6.0d0
 !       bgpt_cs_t(2)=2.0d0/sqrt(3.0d0)
 !       bgpt_cs_t(3)=sqrt(3.0d0)/6.0d0
 !               
 !       wpt_cs_t(1)=1.0d0/sqrt(3.0d0)
 !       wpt_cs_t(2)=1.0d0/sqrt(3.0d0)
 !       wpt_cs_t(3)=1.0d0/sqrt(3.0d0)
 !       
 !   elseif (a==6)then
 !    	agpt_cs_t(1)=-0.725271359
 !       agpt_cs_t(2)=0.725271359
 !       agpt_cs_t(3)=0.0d0
 !       agpt_cs_t(4)=0.0d0
 !       agpt_cs_t(5)=0.337845472
 !       agpt_cs_t(6)=-0.337845472
 !               
 !       bgpt_cs_t(1)=0.158614654
 !       bgpt_cs_t(2)=0.158614654
 !       bgpt_cs_t(3)=1.414821498
 !       bgpt_cs_t(4)=0.187239919
 !       bgpt_cs_t(5)=0.772405443
 !       bgpt_cs_t(6)=0.772405443
 !               
 !       wpt_cs_t(1)=0.109951743655322*sqrt(3.0d0)
 !       wpt_cs_t(2)=0.109951743655322*sqrt(3.0d0)
 !       wpt_cs_t(3)=0.109951743655322*sqrt(3.0d0)
 !       wpt_cs_t(4)=0.223381589678011*sqrt(3.0d0)
 !       wpt_cs_t(5)=0.223381589678011*sqrt(3.0d0)
 !       wpt_cs_t(6)=0.223381589678011*sqrt(3.0d0)
 !       
 !       elseif (a==7)then
 !    	agpt_cs_t(1)=0.0d0
 !       agpt_cs_t(2)=-0.696140478
 !       agpt_cs_t(3)=0.696140478
 !       agpt_cs_t(4)=0.0d0
 !       agpt_cs_t(5)=0.0d0
 !       agpt_cs_t(6)=0.410426192
 !       agpt_cs_t(7)=-0.410426192
 !               
 !       bgpt_cs_t(1)=1.0d0/sqrt(3.0d0)
 !       bgpt_cs_t(2)=0.175433376
 !       bgpt_cs_t(3)=0.175433376
 !       bgpt_cs_t(4)=1.381184054
 !       bgpt_cs_t(5)=0.103430923
 !       bgpt_cs_t(6)=0.814309941
 !       bgpt_cs_t(7)=0.814309941
 !               
 !       wpt_cs_t(1)=0.2250d0*sqrt(3.0d0)
 !       wpt_cs_t(2)=0.125939180544827*sqrt(3.0d0)
 !       wpt_cs_t(3)=0.125939180544827*sqrt(3.0d0)
 !       wpt_cs_t(4)=0.125939180544827*sqrt(3.0d0)
 !       wpt_cs_t(5)=0.132394152788506*sqrt(3.0d0)
 !       wpt_cs_t(6)=0.132394152788506*sqrt(3.0d0)
 !       wpt_cs_t(7)=0.132394152788506*sqrt(3.0d0)
 !       
 !       elseif (a==12)then
 !    	agpt_cs_t(1)=-0.810732956
 !       agpt_cs_t(2)=0.810732956
 !       agpt_cs_t(3)=0.0d0
 !       agpt_cs_t(4)=-0.252139764
 !       agpt_cs_t(5)=0.252139764
 !       agpt_cs_t(6)=0.0d0
 !       
 !       !agpt_cs_t(7)=
 !       !agpt_cs_t(8)=
 !       !agpt_cs_t(9)=
 !       !agpt_cs_t(10)=
 !       !agpt_cs_t(11)=
 !       !agpt_cs_t(12)=
 !               
 !       bgpt_cs_t(1)=0.109273378
 !       bgpt_cs_t(2)=0.109273378
 !       bgpt_cs_t(3)=1.513504051
 !       bgpt_cs_t(4)=0.431777308
 !       bgpt_cs_t(5)=0.431777308
 !       bgpt_cs_t(6)=0.86849619
 !       
 !       !bgpt_cs_t(7)=
 !       !bgpt_cs_t(8)=
 !       !bgpt_cs_t(9)=
 !       !bgpt_cs_t(10)=
 !       !bgpt_cs_t(11)=
 !       !bgpt_cs_t(12)=
 !               
 !       wpt_cs_t(1)=0.050844906370207*sqrt(3.0d0)
 !       wpt_cs_t(2)=0.050844906370207*sqrt(3.0d0)
 !       wpt_cs_t(3)=0.050844906370207*sqrt(3.0d0)
 !       wpt_cs_t(4)=0.116786275726379*sqrt(3.0d0)
 !       wpt_cs_t(5)=0.116786275726379*sqrt(3.0d0)
 !       wpt_cs_t(6)=0.116786275726379*sqrt(3.0d0)
 !       wpt_cs_t(7)=0.082851075618374*sqrt(3.0d0)
 !       wpt_cs_t(8)=0.082851075618374*sqrt(3.0d0)
 !       wpt_cs_t(9)=0.082851075618374*sqrt(3.0d0)
 !       wpt_cs_t(10)=0.082851075618374*sqrt(3.0d0)
 !       wpt_cs_t(11)=0.082851075618374*sqrt(3.0d0)
 !       wpt_cs_t(12)=0.082851075618374*sqrt(3.0d0)
 !      
 !   endif
 !   
	!return
 !   end subroutine   

    subroutine gauss_bm(a)
    use var_analysis
    
        implicit none
        integer::a
        
	        if (a==2)then
     	        gpt_bm(1)=-1.0/sqrt(3.0)
     	        gpt_bm(2)=1.0/sqrt(3.0)
     	        wpt_bm(1)=1.0
     	        wpt_bm(2)=1.0
 	        
 	        elseif (a==3)THEN
     	        gpt_bm(1)=-0.774596669241483
     	        gpt_bm(2)=0.0
     	        gpt_bm(3)=0.774596669241483
     	        wpt_bm(1)=0.555555555555555
     	        wpt_bm(2)=0.888888888888888
                wpt_bm(3)=0.555555555555555
               
            elseif (a==4)THEN
     	        gpt_bm(1)=-0.861136311594053
     	        gpt_bm(2)=-0.339981043584856
     	        gpt_bm(3)=0.339981043584856
                gpt_bm(4)=0.861136311594053
     	        wpt_bm(1)=0.347854845137454
     	        wpt_bm(2)=0.652145154862546
     	        wpt_bm(3)=0.652145154862546
                wpt_bm(4)=0.347854845137454
        
            elseif (a==5)THEN
                gpt_bm(1)=-0.906179845939
                gpt_bm(2)=-0.538469310106
		        gpt_bm(3)=0.0
		        gpt_bm(4)=0.538469310106
		        gpt_bm(5)=0.906179845939
     	        wpt_bm(1)=0.236926885056
		        wpt_bm(2)=0.478628670499
		        wpt_bm(3)=0.568888888889
		        wpt_bm(4)=0.478628670499
		        wpt_bm(5)=0.236926885056
                
          elseif (a==6)THEN
                gpt_bm(1)=-0.9324695142031520
                gpt_bm(2)=-0.6612093864662645
		        gpt_bm(3)=-0.2386191860831969
		        gpt_bm(4)=0.2386191860831969
		        gpt_bm(5)=0.6612093864662645
                gpt_bm(6)=0.9324695142031520
                
     	        wpt_bm(1)=0.1713244923791703
		        wpt_bm(2)=0.3607615730481386
		        wpt_bm(3)=0.4679139345726910
		        wpt_bm(4)=0.4679139345726910
		        wpt_bm(5)=0.3607615730481386
                wpt_bm(6)=0.1713244923791703
                
            elseif (a==7)THEN
                gpt_bm(1)=-0.9491079123427585
                gpt_bm(2)=-0.7415311855993944
		        gpt_bm(3)=-0.4058451513773971
		        gpt_bm(4)=0.0
		        gpt_bm(5)=0.4058451513773971
                gpt_bm(6)=0.7415311855993944
		        gpt_bm(7)=0.9491079123427585
                
     	        wpt_bm(1)=0.1294849661688696
		        wpt_bm(2)=0.2797053914892766
		        wpt_bm(3)=0.3818300505051189
		        wpt_bm(4)=0.4179591836734693
		        wpt_bm(5)=0.3818300505051189
                wpt_bm(6)=0.2797053914892766
		        wpt_bm(7)=0.1294849661688696
                
            elseif (a==8)THEN
                gpt_bm(1)=-0.9602898564975362
                gpt_bm(2)=-0.7966664774136267
		        gpt_bm(3)=-0.5255324099163289
		        gpt_bm(4)=-0.1834346424956498
		        gpt_bm(5)=0.1834346424956498
                gpt_bm(6)=0.5255324099163289
		        gpt_bm(7)=0.7966664774136267
                gpt_bm(8)=0.9602898564975362
                
     	        wpt_bm(1)=0.1012285362903762
		        wpt_bm(2)=0.2223810344533744
		        wpt_bm(3)=0.3137066458778872
		        wpt_bm(4)=0.3626837833783619
		        wpt_bm(5)=0.3626837833783619
                wpt_bm(6)=0.3137066458778872
		        wpt_bm(7)=0.2223810344533744
                wpt_bm(8)=0.1012285362903762
                
            elseif (a==9)THEN
                gpt_bm(1)=-0.9681602395076260
                gpt_bm(2)=-0.8360311073266357
		        gpt_bm(3)=-0.6133714327005903
		        gpt_bm(4)=-0.3242534234038089
		        gpt_bm(5)=0.0000000000000000
                gpt_bm(6)=0.3242534234038089
		        gpt_bm(7)=0.6133714327005903
                gpt_bm(8)=0.8360311073266357
		        gpt_bm(9)=0.9681602395076260
                
     	        wpt_bm(1)=0.0812743883615744
		        wpt_bm(2)=0.1806481606948574
		        wpt_bm(3)=0.2606106964029354
		        wpt_bm(4)=0.3123470770400028
		        wpt_bm(5)=0.3302393550012597
                wpt_bm(6)=0.3123470770400028
		        wpt_bm(7)=0.2606106964029354
                wpt_bm(8)=0.1806481606948574
		        wpt_bm(9)=0.0812743883615744
            endif
	return
    end subroutine

    end module

 !   subroutine gauss(a)
 !   use var_analysis
 !   
 !       implicit none
 !       integer::a
 !       
	!        if (a==2)then
 !    	        gpt(1)=-1.0/sqrt(3.0)
 !    	        gpt(2)=1.0/sqrt(3.0)
 !    	        wpt(1)=1.0
 !    	        wpt(2)=1.0
 !	        
 !	        elseif (a==3)THEN
 !    	        gpt(1)=-0.774596669241483
 !    	        gpt(2)=0.0
 !    	        gpt(3)=0.774596669241483
 !    	        wpt(1)=0.555555555555555
 !    	        wpt(2)=0.888888888888888
 !               wpt(3)=0.555555555555555
 !              
 !           elseif (a==4)THEN
 !    	        gpt(1)=-0.861136311594053
 !    	        gpt(2)=-0.339981043584856
 !    	        gpt(3)=0.339981043584856
 !               gpt(4)=0.861136311594053
 !    	        wpt(1)=0.347854845137454
 !    	        wpt(2)=0.652145154862546
 !    	        wpt(3)=0.652145154862546
 !               wpt(4)=0.347854845137454
 !       
 !           elseif (a==5)THEN
 !               gpt(1)=-0.906179845939
 !               gpt(2)=-0.538469310106
	!	        gpt(3)=0.0
	!	        gpt(4)=0.538469310106
	!	        gpt(5)=0.906179845939
 !    	        wpt(1)=0.236926885056
	!	        wpt(2)=0.478628670499
	!	        wpt(3)=0.568888888889
	!	        wpt(4)=0.478628670499
	!	        wpt(5)=0.236926885056
 !               
 !         elseif (a==6)THEN
 !               gpt(1)=-0.9324695142031520
 !               gpt(2)=-0.6612093864662645
	!	        gpt(3)=-0.2386191860831969
	!	        gpt(4)=0.2386191860831969
	!	        gpt(5)=0.6612093864662645
 !               gpt(6)=0.9324695142031520
 !               
 !    	        wpt(1)=0.1713244923791703
	!	        wpt(2)=0.3607615730481386
	!	        wpt(3)=0.4679139345726910
	!	        wpt(4)=0.4679139345726910
	!	        wpt(5)=0.3607615730481386
 !               wpt(6)=0.1713244923791703
 !               
 !           elseif (a==7)THEN
 !               gpt(1)=-0.9491079123427585
 !               gpt(2)=-0.7415311855993944
	!	        gpt(3)=-0.4058451513773971
	!	        gpt(4)=0.0
	!	        gpt(5)=0.4058451513773971
 !               gpt(6)=0.7415311855993944
	!	        gpt(7)=0.9491079123427585
 !               
 !    	        wpt(1)=0.1294849661688696
	!	        wpt(2)=0.2797053914892766
	!	        wpt(3)=0.3818300505051189
	!	        wpt(4)=0.4179591836734693
	!	        wpt(5)=0.3818300505051189
 !               wpt(6)=0.2797053914892766
	!	        wpt(7)=0.1294849661688696
 !               
 !           elseif (a==8)THEN
 !               gpt(1)=-0.9602898564975362
 !               gpt(2)=-0.7966664774136267
	!	        gpt(3)=-0.5255324099163289
	!	        gpt(4)=-0.1834346424956498
	!	        gpt(5)=0.1834346424956498
 !               gpt(6)=0.5255324099163289
	!	        gpt(7)=0.7966664774136267
 !               gpt(8)=0.9602898564975362
 !               
 !    	        wpt(1)=0.1012285362903762
	!	        wpt(2)=0.2223810344533744
	!	        wpt(3)=0.3137066458778872
	!	        wpt(4)=0.3626837833783619
	!	        wpt(5)=0.3626837833783619
 !               wpt(6)=0.3137066458778872
	!	        wpt(7)=0.2223810344533744
 !               wpt(8)=0.1012285362903762
 !               
 !           elseif (a==9)THEN
 !               gpt(1)=-0.9681602395076260
 !               gpt(2)=-0.8360311073266357
	!	        gpt(3)=-0.6133714327005903
	!	        gpt(4)=-0.3242534234038089
	!	        gpt(5)=0.0000000000000000
 !               gpt(6)=0.3242534234038089
	!	        gpt(7)=0.6133714327005903
 !               gpt(8)=0.8360311073266357
	!	        gpt(9)=0.9681602395076260
 !               
 !    	        wpt(1)=0.0812743883615744
	!	        wpt(2)=0.1806481606948574
	!	        wpt(3)=0.2606106964029354
	!	        wpt(4)=0.3123470770400028
	!	        wpt(5)=0.3302393550012597
 !               wpt(6)=0.3123470770400028
	!	        wpt(7)=0.2606106964029354
 !               wpt(8)=0.1806481606948574
	!	        wpt(9)=0.0812743883615744
 !               
 !           elseif (a==10)THEN
 !               gpt(1)=-0.9739065285171717
 !               gpt(2)=-0.8650633666889845
	!	        gpt(3)=-0.6794095682990244
	!	        gpt(4)=-0.4333953941292471
	!	        gpt(5)=-0.1488743389816312
 !               gpt(6)=0.1488743389816312
	!	        gpt(7)=0.4333953941292471
 !               gpt(8)=0.6794095682990244
	!	        gpt(9)=0.8650633666889845
 !               gpt(10)=0.9739065285171717
 !               
 !    	        wpt(1)=0.0666713443086881
	!	        wpt(2)=0.1494513491505805
	!	        wpt(3)=0.2190863625159820
	!	        wpt(4)=0.2692667193099963
	!	        wpt(5)=0.2955242247147528
 !               wpt(6)=0.2955242247147528
	!	        wpt(7)=0.2692667193099963
 !               wpt(8)=0.2190863625159820
	!	        wpt(9)=0.1494513491505805
 !               wpt(10)=0.0666713443086881
 !               
 !            elseif (a==11)THEN
 !               gpt(1)=-0.9782286581460569
 !               gpt(2)=-0.8870625997680952
	!	        gpt(3)=-0.7301520055740493
	!	        gpt(4)=-0.5190961292068118
	!	        gpt(5)=-0.2695431559523449
 !               gpt(6)=0.0000000000000000
	!	        gpt(7)=0.2695431559523449
 !               gpt(8)=0.5190961292068118
	!	        gpt(9)=0.7301520055740493
 !               gpt(10)=0.8870625997680952
 !               gpt(11)=0.9782286581460569
 !               
 !    	        wpt(1)=0.0556685671161736
	!	        wpt(2)=0.1255803694649046
	!	        wpt(3)=0.1862902109277342
	!	        wpt(4)=0.2331937645919904
	!	        wpt(5)=0.2628045445102466
 !               wpt(6)=0.2729250867779006
	!	        wpt(7)=0.2628045445102466
 !               wpt(8)=0.2331937645919904
	!	        wpt(9)=0.1862902109277342
 !               wpt(10)=0.1255803694649046
 !               wpt(11)=0.0556685671161736
 !               
 !            elseif (a==12)THEN
 !               gpt(1)=-0.9815606342467192
 !               gpt(2)=-0.9041172563704748
	!	        gpt(3)=-0.7699026741943046
	!	        gpt(4)=-0.5873179542866174
	!	        gpt(5)=-0.3678314989981801
 !               gpt(6)=-0.1252334085114689
	!	        gpt(7)=0.1252334085114689
 !               gpt(8)=0.3678314989981801
	!	        gpt(9)=0.5873179542866174
 !               gpt(10)=0.7699026741943046
 !               gpt(11)=0.9041172563704748
 !               gpt(12)=0.9815606342467192
 !               
 !    	        wpt(1)=0.0471753363865118
	!	        wpt(2)=0.1069393259953184
	!	        wpt(3)=0.1600783285433462
	!	        wpt(4)=0.2031674267230659
	!	        wpt(5)=0.2334925365383548
 !               wpt(6)=0.2491470458134027
	!	        wpt(7)=0.2491470458134027
 !               wpt(8)=0.2334925365383548
	!	        wpt(9)=0.2031674267230659
 !               wpt(10)=0.1600783285433462
 !               wpt(11)=0.1069393259953184
 !               wpt(12)=0.0471753363865118
 !               
 !           elseif (a==13)THEN
 !               gpt(1)=-0.9841830547185881d0
 !               gpt(2)=-0.9175983992229779d0
	!	        gpt(3)=-0.8015780907333099d0
	!	        gpt(4)=-0.6423493394403402d0
	!	        gpt(5)=-0.4484927510364468d0
 !               gpt(6)=-0.2304583159551347d0
	!	        gpt(7)=0.0000000000000000d0
 !               gpt(8)=0.2304583159551347d0
	!	        gpt(9)=0.4484927510364468d0
 !               gpt(10)=0.6423493394403402d0
 !               gpt(11)=0.8015780907333099d0
 !               gpt(12)=0.9175983992229779d0
 !               gpt(13)=0.9841830547185881d0
 !               
 !    	        wpt(1)=0.0404840047653158d0
	!	        wpt(2)=0.0921214998377284d0
	!	        wpt(3)=0.1388735102197872d0
	!	        wpt(4)=0.1781459807619457d0
	!	        wpt(5)=0.2078160475368885d0
 !               wpt(6)=0.2262831802628972d0
	!	        wpt(7)=0.23255155323087390d0
 !               wpt(8)=0.2262831802628972d0
	!	        wpt(9)=0.2078160475368885d0
 !               wpt(10)=0.1781459807619457d0
 !               wpt(11)=0.1388735102197872d0
 !               wpt(12)=0.0921214998377284d0
 !               wpt(13)=0.0404840047653158d0
 !           endif       
	!return
 !   end subroutine