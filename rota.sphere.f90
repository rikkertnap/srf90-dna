! -----------------------------------------------------------|
! rotates a given chains conformation                        |  
! pre: xend = input chain                                    |
! post: xendr= rotated chain                                 |
!       return= is_larger_radius                             |
! -----------------------------------------------------------|

! old subroutine rota(xend,xendr,nseg,test,lseg)

subroutine rotation(xend,xendr,nseg,is_larger_radius)
  
    use precision_definition
    use mathconst
    use random
    use volume, only : delta,radius

    implicit none

    integer, intent(in) :: nseg
    real(dp) :: xend(3,nseg+5),xendr(3,nseg+5)
    logical, intent(out) :: is_larger_radius
        
    ! ..local arguments 

    real(dp) :: x(nseg+1),y(nseg+1),z(nseg+1)
    real(dp) :: theta,radius2
    integer :: k,i
    real(dp) :: fac,fac1,fac2,sbe,cbe,sal,cal,sga
    real(dp) :: alfa,gama,cga,a,b,c
     

    theta=pi/6.0d0
    radius2=radius**2
    
    fac=rands(seed)
    fac1=rands(seed)
    fac2=rands(seed)
    alfa=fac*2*pi
    cbe=fac1*2.0-1.0
    gama=fac2*2*pi
    
    sbe=(1-cbe**2)**0.5
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
      
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)
        xendr(1,i)=a*(-cbe*sal*sga+cal*cga) -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
        xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
        xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
    enddo
  
    !  further rotation and translation to grafted onto sphere 
     
    do k=2,nseg+1    
        x(k)=xendr(2,k)*cos(theta)+xendr(3,k)*sin(theta)+radius*sin(theta)
        y(k)=-xendr(2,k)*sin(theta)+xendr(3,k)*cos(theta)+radius*cos(theta)
        z(k)=xendr(1,k)
    enddo  
     
    is_larger_radius=.true.
    i=2 
    do while((i<=(nseg+1)).and.(is_larger_radius)) 
        if((x(i)**2+y(i)**2 + z(i)**2)<radius2) is_larger_radius=.false.
        i=i+1
    enddo

end subroutine rotation


    
