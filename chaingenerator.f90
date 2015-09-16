! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f:                                            |       
! generator chains on spherical surface                         |
! --------------------------------------------------------------|


subroutine make_chains(chainmethod)

    implicit none

    character(len=8) :: chainmethod

    if(chainmethod.eq.'MC') then
        call make_chains_mc()
    elseif(chainmethod.eq.'FILE') then
        call make_chains_file()
    else
        print*,"chainmethod not equalt to MC or FILE"
        stop
    endif

end subroutine make_chains


!     purpose: init of cuantas polymer
!     configurations of polymer chain anchored onto a spherical surface 

subroutine make_chains_mc()
  
    use globals
    use chains
    use random
    use parameters, only : lseg
    use volume

    implicit none 

    !     .. variable and constant declaractions      

    integer :: j,s               ! dummy indices
    integer :: nchain            ! number of rotations
    integer :: maxnchains        ! number of rotations
    integer :: conf              ! counts number of conformations
    real(dp), dimension(:,:,:), allocatable :: chain !  chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp) :: x,y,z,r           ! coordinates
    integer :: ier

    !     .. executable statements 
    !     .. initializations of variables     
    !     .. init of  cuantas polymer configurations  
    !     .. anchoring of polymer chain onto spherical surface

    allocate(chain(3,nseg,200),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in chain : stat =",ier
        stop 
    endif        


    conf=0                    ! counter for conformations
    seed=435672               ! seed for random number generator 
    maxnchains=12

    do while (conf.le.cuantas)
        nchain= 0      ! init zero 
 
        call cadenas(chain,nseg,lseg,nchain,maxnchains) ! chain generator
        !     call cadenas(chain,nchain,maxnchains,nsegAB,lsegAB)	
     
        do j=1,nchain
        
            conf=conf +1
        
            do s=1,nseg        !     transforming form real- to lattice coordinates
                z = chain(1,s,j)
                x = chain(2,s,j)
                y = chain(3,s,j) 
                r = sqrt(x**2 + (y+radius)**2+ z**2) 
                indexchain(conf,s)= int((r-radius)/delta)+1  
            enddo
        
        enddo                  ! end j loop
     
    enddo                     ! end while loop
             
    !     .. end chains generation 
    print*,"Chains generated"      

    
    deallocate(chain)

end subroutine make_chains_mc


subroutine make_chains_file()
  
    !     .. variable and constant declaractions                                                           
    use globals
    use chains
    use random
    use parameters, only : lseg, vpol
    use volume
    use myutils, only : newunit

    implicit none

    integer :: s,rot,t           ! dummy indices 
    integer :: conf,conffile     ! counts number of conformations
    integer :: nsegfile          ! nseg in input file
    integer :: cuantasfile       ! cuantas in input file

    real(dp), dimension(:,:), allocatable :: chain ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp), dimension(:,:), allocatable :: chains_rot !(3,nsegAB+1) ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1

    real(dp)  :: x,y,z,r            ! coordinates
    real(dp)  :: x0,y0,z0           ! origin coordinates

    character(len=10) :: fname
    integer :: ios, un ,ier
    logical :: rottest
    integer :: nchain,rotmax,maxattempts

    !     .. executable statements                                                                          
    !     .. reading in chains  from file  

    write(fname,'(A10)')'chains.dat'
    open(unit=newunit(un),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        stop
    endif
  
    !     .. first read preamble of file
  
    read(un,*)nseg
    read(un,*)cuantasfile
    read(un,*)nsegtypes

    allocate(vpol(nsegtypes),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in vpol : stat =",ier
        stop 
    endif    

    allocate(type_of_monomer(nseg),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in type_of_monomer : stat =",ier
        stop 
    endif  

    do t=1,nsegtypes
        read(un,*)vpol(t)
    enddo    
    do s=1,nseg
        read(un,*)type_of_monomer(s)
    enddo
         
    conf=0                    ! counter for conformations                   
    conffile=0                ! counter for conformations in file
    seed=435672               ! seed for random number generator 
    rotmax=12                 ! maximum number of rotation conf 
    maxattempts=72            ! maximum number of attempts to rotate conf

!    if(nsegfile-1.ne.nseg) then
!        print*,"nseg chain file not equal internal nseg : stop program"
!       stop
!    endif

    allocate(chain(3,nseg),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in chain : stat =",ier
        stop 
    endif    



    do while ((conf.le.cuantas).and.(conffile.le.cuantasfile))
     
        conffile=conffile +1

        read(un,*)x0,y0,z0      ! .. origin

        chain(1,1) = 0.0
        chain(2,1) = 0.0
        chain(3,1) = 0.0

        do s=2,nsegfile        ! .. read form file
            read(un,*)x,y,z
            chain(1,s) = x-x0
            chain(2,s) = y-y0
            chain(3,s) = z-z0
        enddo

        do rot=1,rotmax        !  ..  transforming form real- to lattice coordinates                                  
            rottest=.FALSE.
            nchain=1
            
            do while ((rottest).and.(nchain.lt.maxattempts)) 
                call rotation(chain,chains_rot,nseg,rottest)
                nchain=nchain+1
            enddo
            
            if(rottest) then 
                conf=conf+1 
                do s=1,nseg
                    z = chains_rot(1,s+1)
                    x = chains_rot(2,s+1)
                    y = chains_rot(3,s+1)
                    r = sqrt(x**2 + (y+radius)**2+ z**2)
                    indexchain(conf,s)= int((r-radius)/delta)+1
                enddo
            endif
        enddo                  ! end rot loop  
     
    enddo                     ! end while loop                
  
    !     .. end chains generation 
  
    if(conf.le.cuantas) then
        print*,"Something went wrong"
        print*,"conf     = ",conf," cuantas      = ",cuantas
        print*,"conffile = ",conffile," cuantasfile = ",cuantasfile
        stop
    else
        print*,"Chains generated"
        cuantas=conffile
    endif
  
    deallocate(chain)

    close(un)
  
end subroutine make_chains_file




!subroutine make_sequence_chain(freq)
!  
!  use globals
!  use chains

!  implicit none
  
!  integer :: freq
!  integer :: s
  
!  do s=1,nseg
!     if(mod(s,freq).ne.0) then ! A segment 
!        isAmonomer(s)=.TRUE.
!     else
!        isAmonomer(s)=.FALSE.
!     endif
!  enddo
  
!end subroutine make_sequence_chain

! subroutine make_sequence_chain(freq,chaintype)
  
!   use globals
!   use chains

!   implicit none
  
!   integer :: freq
!   character(len=8)  :: chaintype
  
!   !     .. local variables
  
!   integer :: s
  
!   if(chaintype.eq.'altA') then
!      do s=1,nsegAB
!         if(mod(s,freq).ne.0) then ! A segment
!            isAmonomer(s)=.TRUE.
!         else
!            isAmonomer(s)=.FALSE.
!         endif
!      enddo
!   else if(chaintype.eq.'altB') then
!      do s=1,nsegAB
!         if(mod(s,freq).eq.0) then ! A segment
!            isAmonomer(s)=.TRUE.
!         else
!            isAmonomer(s)=.FALSE.
!         endif
!      enddo
!   else if(chaintype.eq.'diblock') then
!      do s=1,nsegAB
!         if(s.le.freq) then   ! A segment
!            isAmonomer(s)=.TRUE.
!         else
!            isAmonomer(s)=.FALSE.
!         endif
!      enddo
!   else
!      print*,"Wrong chaintype: aborting program"
!      stop
!   endif
  
! end subroutine make_sequence_chain



