!     makes volume elements 
!     for spherical coordinates 

module volume 

    use precision_definition
    implicit none
  
  
    !     .. variables

    real(dp)  :: delta             ! delta  spacing of lattice site in r-direction
    integer :: nr                 ! nr number of lattice sites in r-direction  nr <= nsize
    real(dp)  :: radius             ! radius spherical surface	  
    real(dp), dimension(:), allocatable :: rc ! radial coordinate
    real(dp), dimension(:), allocatable :: G ! geometrical factor 
    real(dp), dimension(:), allocatable :: deltaG ! geometrical factor
    real(dp), dimension(:), allocatable :: Fplus ! factor in Poisson Eq  
    real(dp), dimension(:), allocatable :: Fmin

contains
  
  subroutine allocate_geometry(N)
    implicit none
    integer, intent(in) :: N
    allocate(rc(N))
    allocate(G(N))
    allocate(deltaG(N+1))
    allocate(Fplus(N))
    allocate(Fmin(N))
    
  end subroutine allocate_geometry
  
  
  
  subroutine  make_geometry()
    
    use globals
    
    implicit none
    
    integer ::  i
    real(dp)  ::  vol
    real(dp)  ::  vtest
    real(dp), parameter :: vtol=1.0d-5

  
    vol=0.0d0
    
    do i=1,nr
       rc(i)= (i-0.5d0) * delta + radius ! radial coordinate 
       G(i) =  (rc(i) /radius)**2 ! geometrical factor 
       deltaG(i)=G(i) +(delta*delta)/(12.0d0*radius*radius) ! delta G(i)= (1/delta) \int dr G(r) 
       Fplus(i)= 1.0d0+ delta/rc(i)
       Fmin(i) = 2.0d0-Fplus(i)      ! factors in Poisson Equation
       vol=vol+ deltaG(i)
    enddo
    
    vtest=(4.0d0/3.0d0)*pi*((nr*delta+radius)**3-radius**3)/(4.0d0*pi*(radius**2))
    vtest=vtest/delta
    if(dabs(vtest-vol)>=vtol) then 
       print*,"Warning vtest=",Vtest," not equal to vol=",vol
    end if
    
    return
      
    end subroutine make_geometry
    
  end module volume
  
