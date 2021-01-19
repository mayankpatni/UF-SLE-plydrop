module struc_3D
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! for generating 3D brick element connectivity and coordinates !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine brick3D
    use var_inputdata
    
    implicit none
    integer ::div,cs,nbelem,bne3D_max,bnode3D_max,i,j,elbm,cnt,bne
    real*8  ::t,x_orig,z_orig,x_new,y_new,z_new,tmin,tmax
    
    bnode3D=anodes3D*(ne+1)+anode*ne
    bne3D=ane*ne

    allocate(bcord3D(bnode3D,3),belemcon3D(bne3D,20+1),belemcon3D_type(bne3D))
    bcord3D=0.0d0
    belemcon3D=0
    belemcon3D_type=0
    cnt=0
    tmin=y(1)                                               !yloc (or t) at first node 
    tmax=y(ne*(nne-1)+1)                                    !yloc (or t) at last node 
    do elbm=1,ne
        do i=1,2
            cnt=(anodes3D+anode)*(elbm-1)+anodes3D*(i-1)
            if (i.eq.1) then
                t=y((elbm-1)*(nne-1)+1)
                do j=1,anodes3D
                    x_orig=cscord3D(j,1)
                    z_orig=cscord3D(j,3)
                    call parametric(x_orig,t,tmin,tmax,z_orig,x_new,y_new,z_new)
                    bcord3D(cnt+j,1)=x_new
                    bcord3D(cnt+j,2)=y_new
                    bcord3D(cnt+j,3)=z_new
                enddo
            elseif (i.eq.2) then
                t=0.5d0*(y((elbm-1)*(nne-1)+1)+y(elbm*(nne-1)+1))
                do j=1,anode
                    x_orig=cscord3D(j,1)
                    z_orig=cscord3D(j,3)
                    call parametric(x_orig,t,tmin,tmax,z_orig,x_new,y_new,z_new)
                    bcord3D(cnt+j,1)=x_new
                    bcord3D(cnt+j,2)=y_new
                    bcord3D(cnt+j,3)=z_new
                enddo
            endif
        enddo
        if (elbm.eq.ne) then
            cnt=(anodes3D+anode)*elbm
            t=y(elbm*(nne-1)+1)
            do j=1,anodes3D
                x_orig=cscord3D(j,1)
                z_orig=cscord3D(j,3)
                call parametric(x_orig,t,tmin,tmax,z_orig,x_new,y_new,z_new)
                bcord3D(cnt+j,1)=x_new
                bcord3D(cnt+j,2)=y_new
                bcord3D(cnt+j,3)=z_new
            enddo
        endif
    enddo  
        
    do elbm=1,ne
        do i=1,ane
            bne=(elbm-1)*ane+i
            belemcon3D(bne,1)=bne
            if (cs_ele_type(i,2).eq.4) then            ! quad element  
                belemcon3D_type(bne)=20
                belemcon3D(bne,2:9)=elemcon3D(i,2:9)+(anodes3D+anode)*(elbm-1)                     
                belemcon3D(bne,10:13)=elemcon3D(i,2:5)+(anodes3D+anode)*(elbm-1)+anodes3D       
                belemcon3D(bne,14:21)=elemcon3D(i,2:9)+(anodes3D+anode)*elbm                   
            elseif (cs_ele_type(i,2).eq.3) then        ! tri element
                belemcon3D_type(bne)=15
                belemcon3D(bne,2:4)=elemcon3D(i,2:4)+(anodes3D+anode)*(elbm-1)    
                belemcon3D(bne,6:8)=elemcon3D(i,6:8)+(anodes3D+anode)*(elbm-1)
                belemcon3D(bne,10:12)=elemcon3D(i,2:4)+(anodes3D+anode)*(elbm-1)+anodes3D       
                belemcon3D(bne,14:16)=elemcon3D(i,2:4)+(anodes3D+anode)*elbm    
                belemcon3D(bne,18:20)=elemcon3D(i,6:8)+(anodes3D+anode)*elbm    
            endif
        enddo
    enddo
        
    write(11,*)bnode3D,bne3D
    do i=1,bnode3D
        write(11,*)bcord3D(i,1),bcord3D(i,2),bcord3D(i,3)
    enddo
         
    do i=1,bne3D
        write(11,90)(belemcon3D(i,j),j=1,21)
    enddo
 
    
90 format(28(I3,1X))    
    return
    end subroutine brick3D
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! parametric equation along the beam axis !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine parametric(x_orig,t,tmin,tmax,z_orig,x_new,y_new,z_new)
    use var_inputdata
    
    implicit none
    real*8::x_orig,t,z_orig,x_new,y_new,z_new,x1,x2,x3,x4,x5,x6,x7,x8,x9,x,m,c,tmax,tmin
    real*8::r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,x_1,m_1,x_2,m_2,x_3,m_3
    
    !!!!!!!!!! Prismatic !!!!!!!!!!!!!!!!!
    x_new=x_orig
    y_new=t
    z_new=z_orig
    
    !!!!!!!!!! tapered in z !!!!!!!!!!!!!!!!!
    !rmax=2.0d0
    !rmin=1.0d0
    !x_new=x_orig
    !y_new=t
    !z_new=(z_orig/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    
    !!!!!!!!! tapered in z multiple ply drops !!!!!!!!!!!!!!!!!
    !x1=32.0d0; x2=34.0d0; x3=36.0d0   
    !x4=38.0d0; x5=40.0d0; x6=42.0d0 
    !x7=44.0d0; x8=46.0d0; x9=48.0d0
    !
    !r1=0.25d0; r2=0.50d0; r3=0.75d0; r4=1.00d0; 
    !r5=1.25d0; r6=1.50d0; r7=1.75d0; r8=2.00d0; 
    !r9=2.25d0; r10=2.50d0; r11=2.75d0; r12=3.00d0; 
    !r13=3.25d0; r14=3.50d0; r15=3.75d0; r16=4.00d0
    !  
    !x_new=x_orig
    !y_new=t
    !
    !if (z_orig.le.r4) then
    !    z_new=z_orig
    !    
    !elseif (z_orig.gt.r4.and.z_orig.le.r5) then          ! element number 5 
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x3.le.1e-8) then     ! ply drop
    !        c=z_orig
    !        x=t-x2
    !        m=(z_orig-r4)/2.0d0
    !        z_new=c-m*x
    !    else
    !        z_new=r4
    !    endif
    !    
    !elseif (z_orig.gt.r5.and.z_orig.le.r6) then          ! element number 6
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x3.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x3.gt.1e-8.and.t-x6.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x6.gt.1e-8.and.t-x7.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=t-x6
    !        m=0.25d0/2.0d0
    !        m_1=(c-m*x-r4)/2.0d0
    !        z_new=c-m*x-m_1*x_1
    !    else
    !        z_new=r4
    !    endif
    !    
    !elseif (z_orig.gt.r6.and.z_orig.le.r7) then          ! element number 7
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x3.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x3.gt.1e-8.and.t-x6.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x6.gt.1e-8.and.t-x7.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=t-x6
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x7.gt.1e-8.and.t-x8.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x7-x6
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x8.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x7-x6
    !        x_2=t-x8
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        m_2=(c-m*x-m_1*x_1-r4)/2.0d0
    !        z_new=c-m*x-m_1*x_1-m_2*x_2
    !    else
    !        z_new=r4
    !    endif
    !elseif (z_orig.gt.r7.and.z_orig.le.r8) then          ! element number 8
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x3.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x3.gt.1e-8.and.t-x4.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x4.gt.1e-8.and.t-x5.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=t-x4
    !        m=0.25d0/2.0d0
    !        m_1=(c-m*x-r6)/2.0d0
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x5.gt.1e-8.and.t-x6.le.1e-8) then     
    !        z_new=r6
    !    elseif (t-x6.gt.1e-8.and.t-x7.le.1e-8) then     
    !         m=0.25d0/2.0d0
    !         x=t-x6
    !         z_new=r6-m*x
    !    elseif (t-x7.gt.1e-8.and.t-x8.le.1e-8) then     
    !         z_new=r5
    !    elseif (t-x8.gt.1e-8.and.t-x9.le.1e-8) then     
    !         m=0.25d0/2.0d0
    !         x=t-x8 
    !         z_new=r5-m*x
    !    else
    !        z_new=r4
    !    endif
    !elseif (z_orig.gt.r8.and.z_orig.le.r10) then          ! element number 9 and 10
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x3.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x3.gt.1e-8.and.t-x4.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x4.gt.1e-8.and.t-x5.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=t-x4
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x5.gt.1e-8.and.t-x6.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x5-x4
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x6.gt.1e-8.and.t-x7.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x5-x4
    !        x_2=t-x6
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        m_2=m
    !        z_new=c-m*x-m_1*x_1-m_2*x_2
    !    elseif (t-x7.gt.1e-8.and.t-x8.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x5-x4
    !        x_2=x7-x6
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        m_2=m
    !        z_new=c-m*x-m_1*x_1-m_2*x_2
    !    elseif (t-x8.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x5-x4
    !        x_2=x7-x6
    !        x_3=t-x8
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        m_2=m
    !        m_3=m
    !        z_new=c-m*x-m_1*x_1-m_2*x_2-m_3*x_3
    !    else
    !        c=z_orig
    !        x=x3-x2
    !        x_1=x5-x4
    !        x_2=x7-x6
    !        x_3=x9-x8
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        m_2=m
    !        m_3=m
    !        z_new=c-m*x-m_1*x_1-m_2*x_2-m_3*x_3
    !    endif
    !elseif (z_orig.gt.r10.and.z_orig.le.r11) then          ! element number 11
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x3.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x3.gt.1e-8.and.t-x4.le.1e-8) then     
    !        c=z_orig
    !        x=x3-x2
    !        x_1=t-x3
    !        m=0.25d0/2.0d0
    !        m_1=(c-m*x-r9)/2.0d0
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x4.gt.1e-8.and.t-x5.le.1e-8) then     
    !        c=z_orig
    !        x=t-x4
    !        m=0.25d0/2.0d0
    !        z_new=r9-m*x
    !    elseif (t-x5.gt.1e-8.and.t-x6.le.1e-8) then     
    !        z_new=r8
    !    elseif (t-x6.gt.1e-8.and.t-x7.le.1e-8) then     
    !        c=z_orig
    !        x=t-x6
    !        m=0.25d0/2.0d0
    !        z_new=r8-m*x
    !    elseif (t-x7.gt.1e-8.and.t-x8.le.1e-8) then     
    !        z_new=r7
    !    elseif (t-x8.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=t-x8
    !        m=0.25d0/2.0d0
    !        z_new=r7-m*x
    !    else
    !        z_new=r6
    !    endif
    !    
    !elseif (z_orig.gt.r11.and.z_orig.le.r12) then          ! element number 12
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x5.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x5.gt.1e-8.and.t-x6.le.1e-8) then     
    !        c=z_orig
    !        x=x5-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x6.gt.1e-8.and.t-x7.le.1e-8) then     
    !        c=z_orig
    !        x=x5-x2
    !        x_1=t-x6
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x7.gt.1e-8.and.t-x8.le.1e-8) then     
    !        c=z_orig
    !        x=x5-x2
    !        x_1=x7-x6
    !        x_2=t-x7
    !        m=0.25d0/2.0d0
    !        m_1=m
    !        m_2=(c-m*x-m_1*x_1-r7)/2.0d0
    !        z_new=c-m*x-m_1*x_1-m_2*x_2
    !    elseif (t-x8.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=t-x8
    !        m=0.25d0/2.0d0
    !        z_new=r7-m*x
    !    else
    !        z_new=r6
    !    endif
    !elseif (z_orig.gt.r12.and.z_orig.le.r13) then          ! element number 13
    !    if (t-x2.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x2.gt.1e-8.and.t-x5.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x5.gt.1e-8.and.t-x6.le.1e-8) then     
    !        c=z_orig
    !        x=x5-x2
    !        x_1=t-x5
    !        m=0.25d0/2.0d0
    !        m_1=(c-m*x-r9)/2.0d0
    !        z_new=c-m*x-m_1*x_1
    !    elseif (t-x6.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=t-x6
    !        m=0.25d0/2.0d0
    !        z_new=r9-m*x
    !    else
    !        z_new=r6
    !    endif
    !elseif (z_orig.gt.r13.and.z_orig.le.r14) then          ! element number 14
    !    if (t-x1.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x1.gt.1e-8.and.t-x2.le.1e-8) then     
    !        c=z_orig
    !        x=t-x1
    !        m=(c-r13)/2.0d0
    !        z_new=c-m*x
    !    elseif (t-x2.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=t-x2
    !        m=0.25d0/2.0d0
    !        z_new=r13-m*x
    !    else
    !        z_new=r6
    !    endif
    !elseif (z_orig.gt.r14.and.z_orig.le.r16) then          ! element number 15 and 16
    !    if (t-x1.le.1e-8) then
    !        z_new=z_orig
    !    elseif (t-x1.gt.1e-8.and.t-x9.le.1e-8) then     
    !        c=z_orig
    !        x=t-x1
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    else
    !        c=z_orig
    !        x=x9-x1
    !        m=0.25d0/2.0d0
    !        z_new=c-m*x
    !    endif
    !endif   
    
    !!!!!!!!!!! tapered in z ply drop !!!!!!!!!!!!!!!!!
    !tmin=y(4)                                              
    !tmax=y(7) 
    !rmax=25.0d0
    !rmin=0d0      
    !
    !if (t.lt.tmin) then
    !    x_new=x_orig
    !    y_new=t
    !    z_new=z_orig
    !elseif (t.ge.tmin.and.t.le.tmax) then
    !    x_new=x_orig
    !    y_new=t
    !    if (z_orig.ge.0.0d0) then
    !        if (z_orig.le.rmax) then
    !            z_new=(z_orig/rmax)*(rmax-(t-tmin)*((rmax-rmin)/(tmax-tmin)))
    !        else
    !            z_new=(rmax-(t-tmin)*((rmax-rmin)/(tmax-tmin)))+(z_orig-rmax)
    !        endif
    !    else
    !        if (abs(z_orig).le.rmax) then
    !            z_new=-1.0d0*(abs(z_orig)/rmax)*(rmax-(t-tmin)*((rmax-rmin)/(tmax-tmin)))
    !        else
    !            z_new=-1.0d0*((rmax-(t-tmin)*((rmax-rmin)/(tmax-tmin)))+(abs(z_orig)-rmax))
    !        endif
    !    endif
    !else
    !    x_new=x_orig
    !    y_new=t
    !    if (z_orig.ge.0.0d0) then
    !        if (z_orig.le.rmax) then
    !            z_new=(z_orig/rmax)*(rmax-(tmax-tmin)*((rmax-rmin)/(tmax-tmin)))
    !        else
    !            z_new=(rmax-(tmax-tmin)*((rmax-rmin)/(tmax-tmin)))+(z_orig-rmax)
    !        endif
    !    else
    !        if (abs(z_orig).le.rmax) then
    !            z_new=-1.0d0*(abs(z_orig)/rmax)*(rmax-(tmax-tmin)*((rmax-rmin)/(tmax-tmin)))
    !        else
    !            z_new=-1.0d0*((rmax-(tmax-tmin)*((rmax-rmin)/(tmax-tmin)))+(abs(z_orig)-rmax))
    !        endif
    !    endif
    !endif
    !
    !if (abs(z_new).lt.1e-8) then
    !    z_new=0.0d0
    !endif
    
    
    !!!!!!!!!! tapered in z (wedge) !!!!!!!!!!!!!!!!!
    !rmax=50.0d0
    !rmin=0.0d0
    !x_new=x_orig
    !y_new=t
    !z_new=(z_orig/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !if (z_orig.ge.0.0d0) then
    !    z_new=(z_orig/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !else
    !    z_new=-1.0d0*(abs(z_orig)/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !endif
    
    !!!!!!!!!! tapered in z for I-section !!!!!!!!!!!!!!!!!
    !rmax=450.0d0
    !rmin=50.0d0
    !x_new=x_orig
    !y_new=t
    !if (z_orig.ge.0.0d0) then
    !    if (z_orig.le.rmax) then
    !        z_new=(z_orig/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !    else
    !        z_new=(rmax-t*((rmax-rmin)/(tmax-tmin)))+(z_orig-rmax)
    !    endif
    !else
    !    if (abs(z_orig).le.rmax) then
    !        z_new=-1.0d0*(abs(z_orig)/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !    else
    !        z_new=-1.0d0*((rmax-t*((rmax-rmin)/(tmax-tmin)))+(abs(z_orig)-rmax))
    !    endif
    !endif
    !!!!!!!!!! tapered in z for tapered sandwich !!!!!!!!!!!!!!!!! 
    !rmax=150.0d0
    !rmin=9.45916d0                           !8deg (rmin=9.45916d0), 5deg (rmin=62.511d0), 3deg (rmin=97.59222d0)
    !x_new=x_orig
    !y_new=t
    !if (z_orig.ge.0.0d0) then
    !    if (z_orig.le.rmax) then
    !        z_new=(z_orig/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !    else
    !        z_new=(rmax-t*((rmax-rmin)/(tmax-tmin)))+(z_orig-rmax)
    !    endif
    !else
    !    if (abs(z_orig).le.rmax) then
    !        z_new=-1.0d0*(abs(z_orig)/rmax)*(rmax-t*((rmax-rmin)/(tmax-tmin)))
    !    else
    !        z_new=-1.0d0*((rmax-t*((rmax-rmin)/(tmax-tmin)))+(abs(z_orig)-rmax))
    !    endif
    !endif
    
    return
    end subroutine
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! brickshapefun for 3D brick shape functions and 3D Jacobian !!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine brickshapefun(bne,alpha,eta,beta)
    use var_inputdata
    use var_analysis
    use mod_math
    
    implicit none
    integer::bne
    real*8::alpha,beta,eta,beta_rt3,alpha_rt3,sq3
    real*8::dbsf(20,3),cxyz(20,3),dbsf_wedge(6,3),cxyz_wedge(6,3)
	real*8::l1l2,l1l2_a,l1l2_b,l2l3,l2l3_a,l2l3_b,l3l1,l3l1_a,l3l1_b
	real*8::l1,l2,l3,l1_a,l2_a,l3_a,l1_b,l2_b,l3_b
    
    jac3D=0.0d0
    jac3D_inv=0.0d0
    cxyz=0.0d0
    dbsf=0.0d0
	sq3=sqrt(3.0d0)
    beta_rt3=beta/sq3
    alpha_rt3=alpha/sq3
    
    if (belemcon3D_type(bne).eq.20) then
        cxyz(1:20,1)=bcord3D(belemcon3D(bne,2:21),1)
        cxyz(1:20,2)=bcord3D(belemcon3D(bne,2:21),2)
        cxyz(1:20,3)=bcord3D(belemcon3D(bne,2:21),3)
    
        dbsf(1,1)=0.125d0*(1.0d0-eta)*(1.0d0-beta)*((-1.0d0)*(-2.0d0-alpha-beta-eta)+(1.0d0-alpha)*(-1.0d0))
        dbsf(2,1)=0.125d0*(1.0d0-eta)*(1.0d0-beta)*((1.0d0)*(-2.0d0+alpha-beta-eta)+(1.0d0+alpha)*(1.0d0))
        dbsf(3,1)=0.125d0*(1.0d0-eta)*(1.0d0+beta)*((1.0d0)*(-2.0d0+alpha+beta-eta)+(1.0d0+alpha)*(1.0d0)) 
        dbsf(4,1)=0.125d0*(1.0d0-eta)*(1.0d0+beta)*((-1.0d0)*(-2.0d0-alpha+beta-eta)+(1.0d0-alpha)*(-1.0d0))
        dbsf(5,1)=0.25d0*(1.0d0-eta)*(1.0d0-beta)*(-2.0d0*alpha)
        dbsf(6,1)=0.25d0*(1.0d0-eta)*(1.0d0-beta**2)*(1.0d0)
        dbsf(7,1)=0.25d0*(1.0d0-eta)*(1.0d0+beta)*(-2.0d0*alpha)
        dbsf(8,1)=0.25d0*(1.0d0-eta)*(1.0d0-beta**2)*(-1.0d0)
        dbsf(9,1)=0.25d0*(1.0d0-eta**2)*(1.0d0-beta)*(-1.0d0)
        dbsf(10,1)=0.25d0*(1.0d0-eta**2)*(1.0d0-beta)*(1.0d0)
        dbsf(11,1)=0.25d0*(1.0d0-eta**2)*(1.0d0+beta)*(1.0d0)
        dbsf(12,1)=0.25d0*(1.0d0-eta**2)*(1.0d0+beta)*(-1.0d0)
        dbsf(13,1)=0.125d0*(1.0d0+eta)*(1.0d0-beta)*((-1.0d0)*(-2.0d0-alpha-beta+eta)+(1.0d0-alpha)*(-1.0d0))
        dbsf(14,1)=0.125d0*(1.0d0+eta)*(1.0d0-beta)*((1.0d0)*(-2.0d0+alpha-beta+eta)+(1.0d0+alpha)*(1.0d0))
        dbsf(15,1)=0.125d0*(1.0d0+eta)*(1.0d0+beta)*((1.0d0)*(-2.0d0+alpha+beta+eta)+(1.0d0+alpha)*(1.0d0)) 
        dbsf(16,1)=0.125d0*(1.0d0+eta)*(1.0d0+beta)*((-1.0d0)*(-2.0d0-alpha+beta+eta)+(1.0d0-alpha)*(-1.0d0))
        dbsf(17,1)=0.25d0*(1.0d0+eta)*(1.0d0-beta)*(-2.0d0*alpha)
        dbsf(18,1)=0.25d0*(1.0d0+eta)*(1.0d0-beta**2)*(1.0d0)
        dbsf(19,1)=0.25d0*(1.0d0+eta)*(1.0d0+beta)*(-2.0d0*alpha)
        dbsf(20,1)=0.25d0*(1.0d0+eta)*(1.0d0-beta**2)*(-1.0d0)
        
        dbsf(1,2)=0.125d0*(1.0d0-beta)*(1.0d0-alpha)*((-1.0d0)*(-2.0d0-alpha-beta-eta)+(1.0d0-eta)*(-1.0d0))
        dbsf(2,2)=0.125d0*(1.0d0-beta)*(1.0d0+alpha)*((-1.0d0)*(-2.0d0+alpha-beta-eta)+(1.0d0-eta)*(-1.0d0))
        dbsf(3,2)=0.125d0*(1.0d0+beta)*(1.0d0+alpha)*((-1.0d0)*(-2.0d0+alpha+beta-eta)+(1.0d0-eta)*(-1.0d0)) 
        dbsf(4,2)=0.125d0*(1.0d0+beta)*(1.0d0-alpha)*((-1.0d0)*(-2.0d0-alpha+beta-eta)+(1.0d0-eta)*(-1.0d0))
        dbsf(5,2)=0.25d0*(-1.0d0)*(1.0d0-beta)*(1.0d0-alpha**2)
        dbsf(6,2)=0.25d0*(-1.0d0)*(1.0d0-beta**2)*(1.0d0+alpha)
        dbsf(7,2)=0.25d0*(-1.0d0)*(1.0d0+beta)*(1.0d0-alpha**2)
        dbsf(8,2)=0.25d0*(-1.0d0)*(1.0d0-beta**2)*(1.0d0-alpha)
        dbsf(9,2)=0.25d0*(-2.0d0*eta)*(1.0d0-beta)*(1.0d0-alpha)
        dbsf(10,2)=0.25d0*(-2.0d0*eta)*(1.0d0-beta)*(1.0d0+alpha)
        dbsf(11,2)=0.25d0*(-2.0d0*eta)*(1.0d0+beta)*(1.0d0+alpha)
        dbsf(12,2)=0.25d0*(-2.0d0*eta)*(1.0d0+beta)*(1.0d0-alpha)
        dbsf(13,2)=0.125d0*(1.0d0-beta)*(1.0d0-alpha)*((1.0d0)*(-2.0d0-alpha-beta+eta)+(1.0d0+eta)*(1.0d0))
        dbsf(14,2)=0.125d0*(1.0d0-beta)*(1.0d0+alpha)*((1.0d0)*(-2.0d0+alpha-beta+eta)+(1.0d0+eta)*(1.0d0))
        dbsf(15,2)=0.125d0*(1.0d0+beta)*(1.0d0+alpha)*((1.0d0)*(-2.0d0+alpha+beta+eta)+(1.0d0+eta)*(1.0d0)) 
        dbsf(16,2)=0.125d0*(1.0d0+beta)*(1.0d0-alpha)*((1.0d0)*(-2.0d0-alpha+beta+eta)+(1.0d0+eta)*(1.0d0))
        dbsf(17,2)=0.25d0*(1.0d0)*(1.0d0-beta)*(1.0d0-alpha**2)
        dbsf(18,2)=0.25d0*(1.0d0)*(1.0d0-beta**2)*(1.0d0+alpha)
        dbsf(19,2)=0.25d0*(1.0d0)*(1.0d0+beta)*(1.0d0-alpha**2)
        dbsf(20,2)=0.25d0*(1.0d0)*(1.0d0-beta**2)*(1.0d0-alpha)
    
        dbsf(1,3)=0.125d0*(1.0d0-eta)*(1.0d0-alpha)*((-1.0d0)*(-2.0d0-alpha-beta-eta)+(1.0d0-beta)*(-1.0d0))
        dbsf(2,3)=0.125d0*(1.0d0-eta)*(1.0d0+alpha)*((-1.0d0)*(-2.0d0+alpha-beta-eta)+(1.0d0-beta)*(-1.0d0))
        dbsf(3,3)=0.125d0*(1.0d0-eta)*(1.0d0+alpha)*((1.0d0)*(-2.0d0+alpha+beta-eta)+(1.0d0+beta)*(1.0d0)) 
        dbsf(4,3)=0.125d0*(1.0d0-eta)*(1.0d0-alpha)*((1.0d0)*(-2.0d0-alpha+beta-eta)+(1.0d0+beta)*(1.0d0))
        dbsf(5,3)=0.25d0*(1.0d0-eta)*(-1.0d0)*(1.0d0-alpha**2)
        dbsf(6,3)=0.25d0*(1.0d0-eta)*(-2.0d0*beta)*(1.0d0+alpha)
        dbsf(7,3)=0.25d0*(1.0d0-eta)*(1.0d0)*(1.0d0-alpha**2)
        dbsf(8,3)=0.25d0*(1.0d0-eta)*(-2.0d0*beta)*(1.0d0-alpha)
        dbsf(9,3)=0.25d0*(1.0d0-eta**2)*(-1.0d0)*(1.0d0-alpha)
        dbsf(10,3)=0.25d0*(1.0d0-eta**2)*(-1.0d0)*(1.0d0+alpha)
        dbsf(11,3)=0.25d0*(1.0d0-eta**2)*(1.0d0)*(1.0d0+alpha)
        dbsf(12,3)=0.25d0*(1.0d0-eta**2)*(1.0d0)*(1.0d0-alpha)
        dbsf(13,3)=0.125d0*(1.0d0+eta)*(1.0d0-alpha)*((-1.0d0)*(-2.0d0-alpha-beta+eta)+(1.0d0-beta)*(-1.0d0))
        dbsf(14,3)=0.125d0*(1.0d0+eta)*(1.0d0+alpha)*((-1.0d0)*(-2.0d0+alpha-beta+eta)+(1.0d0-beta)*(-1.0d0))
        dbsf(15,3)=0.125d0*(1.0d0+eta)*(1.0d0+alpha)*((1.0d0)*(-2.0d0+alpha+beta+eta)+(1.0d0+beta)*(1.0d0)) 
        dbsf(16,3)=0.125d0*(1.0d0+eta)*(1.0d0-alpha)*((1.0d0)*(-2.0d0-alpha+beta+eta)+(1.0d0+beta)*(1.0d0))
        dbsf(17,3)=0.25d0*(1.0d0+eta)*(-1.0d0)*(1.0d0-alpha**2)
        dbsf(18,3)=0.25d0*(1.0d0+eta)*(-2.0d0*beta)*(1.0d0+alpha)
        dbsf(19,3)=0.25d0*(1.0d0+eta)*(1.0d0)*(1.0d0-alpha**2)
        dbsf(20,3)=0.25d0*(1.0d0+eta)*(-2.0d0*beta)*(1.0d0-alpha)
    
	elseif (belemcon3D_type(bne).eq.15) then
        cxyz(1:3,1)=bcord3D(belemcon3D(bne,2:4),1)
        cxyz(5:7,1)=bcord3D(belemcon3D(bne,6:8),1)
        cxyz(9:11,1)=bcord3D(belemcon3D(bne,10:12),1)
        cxyz(13:15,1)=bcord3D(belemcon3D(bne,14:16),1)
        cxyz(17:19,1)=bcord3D(belemcon3D(bne,18:20),1)
        cxyz(1:3,2)=bcord3D(belemcon3D(bne,2:4),2)
        cxyz(5:7,2)=bcord3D(belemcon3D(bne,6:8),2)
        cxyz(9:11,2)=bcord3D(belemcon3D(bne,10:12),2)
        cxyz(13:15,2)=bcord3D(belemcon3D(bne,14:16),2)
        cxyz(17:19,2)=bcord3D(belemcon3D(bne,18:20),2)
        cxyz(1:3,3)=bcord3D(belemcon3D(bne,2:4),3)
        cxyz(5:7,3)=bcord3D(belemcon3D(bne,6:8),3)
        cxyz(9:11,3)=bcord3D(belemcon3D(bne,10:12),3)
        cxyz(13:15,3)=bcord3D(belemcon3D(bne,14:16),3)
        cxyz(17:19,3)=bcord3D(belemcon3D(bne,18:20),3)
		
		
		l1=0.5d0*(1.0d0-alpha-beta/sq3)  
		l2=0.5d0*(1.0d0+alpha-beta/sq3) 
		l3=beta/sq3
		l1_a=-0.5d0
		l2_a=0.5d0
		l3_a=0.0d0
		l1_b=-0.5d0*(1.0d0/sq3)
		l2_b=-0.5d0*(1.0d0/sq3)
		l3_b=1.0d0/sq3
        l1l2=0.25d0*((1-beta/sq3)**2-alpha**2)
        l1l2_a=0.25d0*(-2.0d0*alpha)
        l1l2_b=0.25d0*(2.0d0*(1-beta/sq3)*(-1.0d0/sq3))
        l2l3=0.5d0*(beta/sq3+alpha*beta/sq3-(beta/sq3)**2)
        l2l3_a=0.5d0*(beta/sq3)
        l2l3_b=0.5d0*(1.0d0/sq3+alpha/sq3-2.0d0*beta/3.0d0)
        l3l1=0.5d0*(beta/sq3-alpha*beta/sq3-(beta/sq3)**2)
        l3l1_a=-0.5d0*(beta/sq3)
        l3l1_b=0.5d0*(1.0d0/sq3-alpha/sq3-2.0d0*beta/3.0d0)
		
		dbsf(1,1)=0.5d0*((1.0d0-eta)*(4.0d0*l1*l1_a-l1_a)-(1.0d0-eta**2)*l1_a)
		dbsf(1,2)=0.5d0*((-1.0d0)*(2.0d0*l1**2-l1)-(-2.0d0*eta)*l1)
		dbsf(1,3)=0.5d0*((1.0d0-eta)*(4.0d0*l1*l1_b-l1_b)-(1.0d0-eta**2)*l1_b)
		
		dbsf(2,1)=0.5d0*((1.0d0-eta)*(4.0d0*l2*l2_a-l2_a)-(1.0d0-eta**2)*l2_a)
		dbsf(2,2)=0.5d0*((-1.0d0)*(2.0d0*l2**2-l2)-(-2.0d0*eta)*l2)
		dbsf(2,3)=0.5d0*((1.0d0-eta)*(4.0d0*l2*l2_b-l2_b)-(1.0d0-eta**2)*l2_b)
		
		dbsf(3,1)=0.5d0*((1.0d0-eta)*(4.0d0*l3*l3_a-l3_a)-(1.0d0-eta**2)*l3_a)
		dbsf(3,2)=0.5d0*((-1.0d0)*(2.0d0*l3**2-l3)-(-2.0d0*eta)*l3)
		dbsf(3,3)=0.5d0*((1.0d0-eta)*(4.0d0*l3*l3_b-l3_b)-(1.0d0-eta**2)*l3_b)
		
		dbsf(5,1)=2.0d0*l1l2_a*(1.0d0-eta)
		dbsf(5,2)=2.0d0*l1l2*(-1.0d0)
		dbsf(5,3)=2.0d0*l1l2_b*(1.0d0-eta)
		
		dbsf(6,1)=2.0d0*l2l3_a*(1.0d0-eta)
		dbsf(6,2)=2.0d0*l2l3*(-1.0d0)
		dbsf(6,3)=2.0d0*l2l3_b*(1.0d0-eta)
		
		dbsf(7,1)=2.0d0*l3l1_a*(1.0d0-eta)
		dbsf(7,2)=2.0d0*l3l1*(-1.0d0)
		dbsf(7,3)=2.0d0*l3l1_b*(1.0d0-eta)
		
		dbsf(9,1)=l1_a*(1.0d0-eta**2)
		dbsf(9,2)=l1*(-2.0d0*eta)
		dbsf(9,3)=l1_b*(1.0d0-eta**2)
		
		dbsf(10,1)=l2_a*(1.0d0-eta**2)
		dbsf(10,2)=l2*(-2.0d0*eta)
		dbsf(10,3)=l2_b*(1.0d0-eta**2)
		
		dbsf(11,1)=l3_a*(1.0d0-eta**2)
		dbsf(11,2)=l3*(-2.0d0*eta)
		dbsf(11,3)=l3_b*(1.0d0-eta**2)
		
		dbsf(13,1)=0.5d0*((1.0d0+eta)*(4.0d0*l1*l1_a-l1_a)-(1.0d0-eta**2)*l1_a)
		dbsf(13,2)=0.5d0*((1.0d0)*(2.0d0*l1**2-l1)-(-2.0d0*eta)*l1)
		dbsf(13,3)=0.5d0*((1.0d0+eta)*(4.0d0*l1*l1_b-l1_b)-(1.0d0-eta**2)*l1_b)
		
		dbsf(14,1)=0.5d0*((1.0d0+eta)*(4.0d0*l2*l2_a-l2_a)-(1.0d0-eta**2)*l2_a)
		dbsf(14,2)=0.5d0*((1.0d0)*(2.0d0*l2**2-l2)-(-2.0d0*eta)*l2)
		dbsf(14,3)=0.5d0*((1.0d0+eta)*(4.0d0*l2*l2_b-l2_b)-(1.0d0-eta**2)*l2_b)
		
		dbsf(15,1)=0.5d0*((1.0d0+eta)*(4.0d0*l3*l3_a-l3_a)-(1.0d0-eta**2)*l3_a)
		dbsf(15,2)=0.5d0*((1.0d0)*(2.0d0*l3**2-l3)-(-2.0d0*eta)*l3)
		dbsf(15,3)=0.5d0*((1.0d0+eta)*(4.0d0*l3*l3_b-l3_b)-(1.0d0-eta**2)*l3_b)
		
		dbsf(17,1)=2.0d0*l1l2_a*(1.0d0+eta)
		dbsf(17,2)=2.0d0*l1l2*(1.0d0)
		dbsf(17,3)=2.0d0*l1l2_b*(1.0d0+eta)
		
		dbsf(18,1)=2.0d0*l2l3_a*(1.0d0+eta)
		dbsf(18,2)=2.0d0*l2l3*(1.0d0)
		dbsf(18,3)=2.0d0*l2l3_b*(1.0d0+eta)
		
		dbsf(19,1)=2.0d0*l3l1_a*(1.0d0+eta)
		dbsf(19,2)=2.0d0*l3l1*(1.0d0)
		dbsf(19,3)=2.0d0*l3l1_b*(1.0d0+eta)
        
        !dbsf(1,1)=0.25d0*((1.0d0-eta)*(2.0d0*alpha+2.0d0*beta_rt3-1.0d0)+(1.0d0-eta**2)*(1.0d0))
        !dbsf(2,1)=0.25d0*((1.0d0-eta)*(2.0d0*alpha-2.0d0*beta_rt3+1.0d0)-(1.0d0-eta**2)*(1.0d0))
        !dbsf(3,1)=0.0d0
        !dbsf(4,1)=0.0d0
        !dbsf(5,1)=0.5d0*(1.0d0-eta)*(-2.0d0*alpha)
        !dbsf(6,1)=(1.0d0-eta)*(beta_rt3)
        !dbsf(7,1)=(1.0d0-eta)*(-1.0d0*beta_rt3)
        !dbsf(8,1)=0.0d0
        !dbsf(9,1)=0.5d0*(1.0d0-eta**2)*(-1.0d0)
        !dbsf(10,1)=0.5d0*(1.0d0-eta**2)*(1.0d0)
        !dbsf(11,1)=0.0d0
        !dbsf(12,1)=0.0d0
        !dbsf(13,1)=0.25d0*((1.0d0+eta)*(2.0d0*alpha+2.0d0*beta_rt3-1.0d0)+(1.0d0-eta**2)*(1.0d0))
        !dbsf(14,1)=0.25d0*((1.0d0+eta)*(2.0d0*alpha-2.0d0*beta_rt3+1.0d0)-(1.0d0-eta**2)*(1.0d0))
        !dbsf(15,1)=0.0d0
        !dbsf(16,1)=0.0d0
        !dbsf(17,1)=0.5d0*(1.0d0+eta)*(-2.0d0*alpha)
        !dbsf(18,1)=(1.0d0+eta)*(beta_rt3)
        !dbsf(19,1)=(1.0d0+eta)*(-1.0d0*beta_rt3)
        !dbsf(20,1)=0.0d0
        !
        !dbsf(1,2)=0.25d0*((-1.0d0)*(alpha+beta_rt3)*(alpha+beta_rt3-1.0d0)+(-2.0d0*eta)*(alpha+beta_rt3-1.0d0))
        !dbsf(2,2)=0.25d0*((-1.0d0)*(alpha-beta_rt3)*(alpha-beta_rt3+1.0d0)-(-2.0d0*eta)*(alpha-beta_rt3+1.0d0))
        !dbsf(3,2)=0.5d0*((-1.0d0)*(2.0d0*beta_rt3**2-beta_rt3)-(-2.0d0*eta)*(beta_rt3))
        !dbsf(4,2)=0.0d0
        !dbsf(5,2)=0.5d0*(-1.0d0)*(1.0d0-alpha-beta_rt3)*(1.0d0+alpha-beta_rt3)
        !dbsf(6,2)=(-1.0d0)*(beta_rt3-beta_rt3**2+alpha*beta_rt3)
        !dbsf(7,2)=(-1.0d0)*(beta_rt3-beta_rt3**2-alpha*beta_rt3)
        !dbsf(8,2)=0.0d0
        !dbsf(9,2)=0.5d0*(-2.0d0*eta)*(1.0d0-alpha-beta_rt3)
        !dbsf(10,2)=0.5d0*(-2.0d0*eta)*(1.0d0+alpha-beta_rt3)
        !dbsf(11,2)=(-2.0d0*eta)*beta_rt3
        !dbsf(12,2)=0.0d0
        !dbsf(13,2)=0.25d0*((1.0d0)*(alpha+beta_rt3)*(alpha+beta_rt3-1.0d0)+(-2.0d0*eta)*(alpha+beta_rt3-1.0d0))
        !dbsf(14,2)=0.25d0*((1.0d0)*(alpha-beta_rt3)*(alpha-beta_rt3+1.0d0)-(-2.0d0*eta)*(alpha-beta_rt3+1.0d0))
        !dbsf(15,2)=0.5d0*((1.0d0)*(2.0d0*beta_rt3**2-beta_rt3)-(-2.0d0*eta)*(beta_rt3))
        !dbsf(16,2)=0.0d0
        !dbsf(17,2)=0.5d0*(1.0d0)*(1.0d0-alpha-beta_rt3)*(1.0d0+alpha-beta_rt3)
        !dbsf(18,2)=(1.0d0)*(beta_rt3-beta_rt3**2+alpha*beta_rt3)
        !dbsf(19,2)=(1.0d0)*(beta_rt3-beta_rt3**2-alpha*beta_rt3)
        !dbsf(20,2)=0.0d0
        !
        !dbsf(1,3)=0.25d0*((1.0d0-eta)*((2.0d0/3.0d0)*beta+2.0d0*alpha_rt3-1.0d0/sqrt(3.0d0))+(1.0d0-eta**2)*(1.0d0/sqrt(3.0d0)))
        !dbsf(2,3)=0.25d0*((1.0d0-eta)*((2.0d0/3.0d0)*beta-2.0d0*alpha_rt3-1.0d0/sqrt(3.0d0))-(1.0d0-eta**2)*(-1.0d0/sqrt(3.0d0)))
        !dbsf(3,3)=0.5d0*((1.0d0-eta)*((4.0d0/3.0d0)*beta-1.0d0/sqrt(3.0d0))-(1.0d0-eta**2)*(1.0d0/sqrt(3.0d0)))
        !dbsf(4,3)=0.0d0
        !dbsf(5,3)=0.5d0*(1.0d0-eta)*((2.0d0/3.0d0)*beta-2.0d0/sqrt(3.0d0))
        !dbsf(6,3)=(1.0d0-eta)*(1.0d0/sqrt(3.0d0)-(2.0d0/3.0d0)*beta+alpha_rt3)
        !dbsf(7,3)=(1.0d0-eta)*(1.0d0/sqrt(3.0d0)-(2.0d0/3.0d0)*beta-alpha_rt3)
        !dbsf(8,3)=0.0d0
        !dbsf(9,3)=0.5d0*(1.0d0-eta**2)*(-1.0d0/sqrt(3.0d0))
        !dbsf(10,3)=0.5d0*(1.0d0-eta**2)*(-1.0d0/sqrt(3.0d0))
        !dbsf(11,3)=(1.0d0-eta**2)/sqrt(3.0d0)
        !dbsf(12,3)=0.0d0
        !dbsf(13,3)=0.25d0*((1.0d0+eta)*((2.0d0/3.0d0)*beta+2.0d0*alpha_rt3-1.0d0/sqrt(3.0d0))+(1.0d0-eta**2)*(1.0d0/sqrt(3.0d0)))
        !dbsf(14,3)=0.25d0*((1.0d0+eta)*((2.0d0/3.0d0)*beta-2.0d0*alpha_rt3-1.0d0/sqrt(3.0d0))-(1.0d0-eta**2)*(-1.0d0/sqrt(3.0d0)))
        !dbsf(15,3)=0.5d0*((1.0d0+eta)*((4.0d0/3.0d0)*beta-1.0d0/sqrt(3.0d0))-(1.0d0-eta**2)*(1.0d0/sqrt(3.0d0)))
        !dbsf(16,3)=0.0d0
        !dbsf(17,3)=0.5d0*(1.0d0+eta)*((2.0d0/3.0d0)*beta-2.0d0/sqrt(3.0d0))
        !dbsf(18,3)=(1.0d0+eta)*(1.0d0/sqrt(3.0d0)-(2.0d0/3.0d0)*beta+alpha_rt3)
        !dbsf(19,3)=(1.0d0+eta)*(1.0d0/sqrt(3.0d0)-(2.0d0/3.0d0)*beta-alpha_rt3)
        !dbsf(20,3)=0.0d0
    endif
    
    jac3D=matmul(transpose(dbsf),cxyz)
    det_jac_3D=(jac3D(1,1)*(jac3D(2,2)*jac3D(3,3)-jac3D(2,3)*jac3D(3,2))-jac3D(1,2)*(jac3D(2,1)*jac3D(3,3)-jac3D(2,3)*jac3D(3,1))+jac3D(1,3)*(jac3D(2,1)*jac3D(3,2)-jac3D(2,2)*jac3D(3,1))) 

    call matinv3(jac3D,jac3D_inv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !cxyz_wedge(1:4,1)=bcord3D(belemcon3D(bne,2:5),1)
    !cxyz_wedge(1:4,2)=bcord3D(belemcon3D(bne,2:5),2)
    !cxyz_wedge(1:4,3)=bcord3D(belemcon3D(bne,2:5),3)
    !cxyz_wedge(5,1)=bcord3D(belemcon3D(bne,27),1)
    !cxyz_wedge(5,2)=bcord3D(belemcon3D(bne,27),2)
    !cxyz_wedge(5,3)=bcord3D(belemcon3D(bne,27),3)
    !cxyz_wedge(6,1)=bcord3D(belemcon3D(bne,25),1)
    !cxyz_wedge(6,2)=bcord3D(belemcon3D(bne,25),2)
    !cxyz_wedge(6,3)=bcord3D(belemcon3D(bne,25),3)
    
    !bsf_wedge(1,1)=0.125d0*(1.0d0-alpha)*(1.0d0-eta)*(1.0d0-beta)
    !dbsf_wedge(1,1)=-0.125d0*(1.0d0-eta)*(1.0d0-beta)
    !dbsf_wedge(1,2)=-0.125d0*(1.0d0-alpha)*(1.0d0-beta)
    !dbsf_wedge(1,3)=-0.125d0*(1.0d0-alpha)*(1.0d0-eta)
    
    !bsf_wedge(2,1)=0.125d0*(1.0d0+alpha)*(1.0d0-eta)*(1.0d0-beta)
    !dbsf_wedge(2,1)=0.125d0*(1.0d0-eta)*(1.0d0-beta)
    !dbsf_wedge(2,2)=-0.125d0*(1.0d0+alpha)*(1.0d0-beta)
    !dbsf_wedge(2,3)=-0.125d0*(1.0d0+alpha)*(1.0d0-eta)
    
    !bsf_wedge(3,1)=0.125d0*(1.0d0+alpha)*(1.0d0-eta)*(1.0d0+beta)
    !dbsf_wedge(3,1)=0.125d0*(1.0d0-eta)*(1.0d0+beta)
    !dbsf_wedge(3,2)=-0.125d0*(1.0d0+alpha)*(1.0d0+beta)
    !dbsf_wedge(3,3)=0.125d0*(1.0d0+alpha)*(1.0d0-eta)
    
    !bsf_wedge(4,1)=0.125d0*(1.0d0-alpha)*(1.0d0-eta)*(1.0d0+beta)
    !dbsf_wedge(4,1)=-0.125d0*(1.0d0-eta)*(1.0d0+beta)
    !dbsf_wedge(4,2)=-0.125d0*(1.0d0-alpha)*(1.0d0+beta)
    !dbsf_wedge(4,3)=0.125d0*(1.0d0-alpha)*(1.0d0-eta)
    
    !bsf_wedge(5,1)=0.25d0*(1.0d0-alpha)*(1.0d0+eta)
    !dbsf_wedge(5,1)=-0.25d0*(1.0d0+eta)
    !dbsf_wedge(5,2)=0.25d0*(1.0d0-alpha)
    !dbsf_wedge(5,3)=0.0d0
    
    !bsf_wedge(6,1)=0.25d0*(1.0d0+alpha)*(1.0d0+eta)
    !dbsf_wedge(6,1)=0.25d0*(1.0d0+eta)
    !dbsf_wedge(6,2)=0.25d0*(1.0d0+alpha)
    !dbsf_wedge(6,3)=0.0d0
    
   ! if (bne.eq.bne3D) then
       ! jac3D=matmul(transpose(dbsf_wedge),cxyz_wedge)
    !else
    !    jac3D=matmul(transpose(dbsf),cxyz)
   ! endif
    
    !det_jac_3D=(jac3D(1,1)*(jac3D(2,2)*jac3D(3,3)-jac3D(2,3)*jac3D(3,2))-jac3D(1,2)*(jac3D(2,1)*jac3D(3,3)-jac3D(2,3)*jac3D(3,1))+jac3D(1,3)*(jac3D(2,1)*jac3D(3,2)-jac3D(2,2)*jac3D(3,1))) 
    
    !call matinv3(jac3D,jac3D_inv)
    
    return
    end subroutine
    
    
end module
    