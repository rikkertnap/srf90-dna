module field
  
    !     .. variables
    
    use precision_definition
    implicit none
      
    real(dp), dimension(:), allocatable :: xpol  ! total volume fraction of polymer on sphere
    real(dp), dimension(:,:), allocatable :: rhopol! density  monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: xpi     ! volume fraction solvent with no chi-parameters = exp(-beta pi vsol) 
    real(dp), dimension(:), allocatable :: xsol    ! volume fraction solvent
    real(dp), dimension(:), allocatable :: psi     ! electrostatic potential 
    real(dp), dimension(:), allocatable :: xNa     ! volume fraction of positive Na+ ion
    real(dp), dimension(:), allocatable :: xK      ! volume fraction of positive K+ ion
    real(dp), dimension(:), allocatable :: xCa     ! volume fraction of positive Ca2+ ion
    real(dp), dimension(:), allocatable :: xNaCl   ! volume fraction of NaCl ion pair
    real(dp), dimension(:), allocatable :: xKCl    ! volume fraction of KCl  ion pair
    real(dp), dimension(:), allocatable :: xCl     ! volume fraction of negative ion
    real(dp), dimension(:), allocatable :: xHplus  ! volume fraction of Hplus
    real(dp), dimension(:), allocatable :: xOHmin  ! volume fraction of OHmin 
    real(dp), dimension(:), allocatable :: rhoq    ! total charge density in units of vsol
    real(dp), dimension(:), allocatable :: qpol    ! charge density of polymer
    real(dp), dimension(:,:), allocatable :: fdis  ! degree of dissociation 
   
    real(dp) :: q              ! normalization partion fnc polymer 
   
    contains

    subroutine allocate_field(N,nsegtypes)
        implicit none
    
        integer, intent(in) :: N, nsegtypes
        
        allocate(xpol(N))
        allocate(rhopol(N,nsegtypes))
        allocate(xpi(N))
        allocate(xsol(N))
        allocate(psi(N+1))
        allocate(xNa(N))
        allocate(xK(N))
        allocate(xCa(N))
        allocate(xNaCl(N)) 
        allocate(xKCl(N)) 
        allocate(xCl(N)) 
        allocate(xHplus(N))
        allocate(xOHmin(N))
        allocate(rhoq(N))
        allocate(qpol(N))
        allocate(fdis(N,nsegtypes))
        
    
    end subroutine allocate_field
  

    subroutine deallocate_field()
        implicit none


        deallocate(xpol)
        deallocate(rhopol)
        deallocate(xsol)
        deallocate(psi)
        deallocate(xNa)
        deallocate(xK)
        deallocate(xCa)
        deallocate(xNaCl) 
        deallocate(xKCl) 
        deallocate(xCl) 
        deallocate(xHplus)
        deallocate(xOHmin)
        deallocate(rhoq)
        deallocate(qpol)
        deallocate(fdis)


    end subroutine deallocate_field



    !     .. compute average height of tethered layer 
    !     .. first moment of density profile 

    subroutine average_height()

        use globals
        use parameters
        use volume

        implicit none

        integer :: i

        real(dp) :: zero,first

        first=0.0d0               ! first moment 
        zero=0.0d0                ! zero moment  
        do i=1,nr
            zero=zero+xpol(i)*deltaG(i)
            first=first+xpol(i)*rc(i)*deltaG(i)
        enddo

        if(zero>0.0) then 
            avheightpol=first/zero
        else
            avheightpol=0.0d0
        endif

        if(isNaN(avheightpol)) print*,"heightpol NaN"

    end subroutine average_height

    subroutine average_charge_polymer()

        use globals
        use volume
        use parameters

        implicit none

        integer :: i,t
        real(dp)  :: Asurf

        Asurf=4.0d0*pi*(radius**2)

        avqpol=0.0d0
        do t=1,nsegtypes
            do i=1,nr
                avqpol=avqpol+(zpol(t,1)*fdis(i,t)+zpol(t,2)*(1.0d0-fdis(i,t)))*rhopol(i,t)*deltaG(i)
            enddo
        enddo    
        avqpol=avqpol*Asurf*delta
        

    end subroutine average_charge_polymer

    ! .. post : return average charge of state of 
    !   of polymers

    subroutine average_fdis_polymer()

        use globals
        use volume
        use parameters
        use chains

        implicit none 

        integer :: i,s,t
        real(dp)  :: Asurf
        real(dp)  :: npol_t  ! .. number of monomers of type t

        Asurf=4.0d0*pi*(radius**2)

        do t=1,nsegtypes
            ! .. number monomors of type t  
            npol_t=0
            do s=1,nseg
                if(ismonomer_of_type(s,t).eqv..true.) then
                    npol_t=npol_t+1
                endif
            enddo
            if(npol_t.ne.0) then
                avfdis(t)=0.0d0
                do i=1,nr
                        avfdis(t)=avfdis(t)+fdis(i,t)*rhopol(i,t)*deltaG(i) 
                enddo
                avfdis(t)=avfdis(t)*delta/(sigma*delta*npol_t)
            else
                avfdis(t)=0.0d0
            endif
        enddo

    end subroutine average_fdis_polymer


    pure function MyisNaN(x) Result(l) 
        implicit none
        real(dp), intent(in) :: x
        logical :: l     
        if (x /= x) then
            l=.true.
        else
            l=.false.
        endif 
    end function MyisNaN 

    end module field

