! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f:                                            |       
! generator chains on spherical surface                         |
! --------------------------------------------------------------|


subroutine make_chains(chainmethod)
 
    use matrices

    implicit none

    character(len=8) :: chainmethod

    if(chainmethod.eq.'MC') then
        call init_matrices() ! make transition matrices
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
    use parameters, only : chainfname,spacer
    use volume, only : radius,delta
    use myutils, only : newunit

    implicit none

    integer :: s,rot,t           ! dummy indices 
    integer :: conf,conffile     ! counts number of conformations
    integer :: conffile_reject    ! counts number of conformations from file rejected 
    integer :: count_reject      ! flag for counting  rejected conffile_reject
    integer :: numconf           ! total number of conformations 
    integer :: nsegfile          ! nseg in input file
    integer :: cuantasfile       ! cuantas in input file
    integer :: nsegtypesfile     ! nsegtypes in input file
    
          
    real(dp), dimension(:,:), allocatable :: chain ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp), dimension(:,:), allocatable :: chains_rot !(3,nseg+5) ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1

    real(dp)  :: x,y,z,r            ! coordinates
    real(dp)  :: x0,y0,z0           ! origin coordinates

    character(len=10) :: fname
    integer :: ios, un ,ier
    logical :: rottest
    integer :: nchain,rotmax,maxattempts
    real(dp) :: reject_rate

    character(len=80) :: letter,comment,istr,str ! auxilary variable for file to be read in
    integer, dimension(:), allocatable :: type_of_monomer_file


    !     .. executable statements                                                                          
  
 
    !     .. reading in of chains  from file
    open(unit=newunit(un),file=chainfname,iostat=ios,status='old')
    if(ios/=0 ) then
        write(istr,'(I2)')ios
        str='Error opening file '//trim(adjustl(chainfname))//' : iostat = '//istr
        print*,str
        stop
    endif

    !     .. first read preamble of file
  
    read(un,*)nsegfile
    read(un,*)cuantasfile
    read(un,*)nsegtypesfile

    !     .. check preamble

    if(nsegfile/=nseg) then 
        print*,"Variable nseg read in from file does not match value in input.in"
        stop 
    endif
   
    if(nsegtypesfile/=nsegtypes) then 
        print*,"Variable nsegtypes read in from file does not match value in input.in"
        stop 
    endif

    allocate(type_of_monomer_file(nseg),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in type_of_monomer_file : stat =",ier
        stop      
    endif  

    do s=1,nsegfile
        read(un,*)type_of_monomer_file(s),type_of_monomer_char(s)
    enddo   

    allocate(chain(3,nseg+5),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in chain : stat =",ier
        stop 
    endif    

    allocate(chains_rot(3,nseg+5),stat=ier)
    if(ier.ne.0) then
        print*,"Allocation error in chain_rot : stat =",ier
        stop 
    endif 

         
    conf=0                    ! counter for accepted conformations                   
    conffile=0                ! counter for conformations read in file
    conffile_reject=0         ! number of conf read in rejected after rotmax*maxattempts trials  
    seed=435672               ! seed for random number generator 
    rotmax=5                  ! maximum number of rotation conf 
    maxattempts=72            ! maximum number of attempts to rotate conf
    y0=spacer                 ! spacer to add anchor point
    numconf=0                 ! total number of conformations

    do while ((conf<=cuantas).and.(conffile<=cuantasfile))
    
        ! .. origin

        chain(1,1) = 0.0
        chain(2,1) = 0.0
        chain(3,1) = 0.0

        
        !  xcoor = index 2 : ycoor = index 3 and zcoor = index 1
        !  coordinate transformation: 
        !  x'''=Rx''=RMx'=RM(x+T)  
        !  1) Transtation 2) mirror y=-(y-y0) and 3)  rotate over z-axis z-> -y and y-> z
 
        do s=1,nsegfile        
            read(un,*)x,y,z       !  read form file
            chain(2,s+1) = x            ! x    
            chain(3,s+1) = -z           ! -z 
            chain(1,s+1) = -(y-y0)      ! -(y-y0) 
        enddo
             
        conffile=conffile +1
        
        do rot=1,rotmax        !  ..  transforming form real- to lattice coordinates                                  
            rottest=.FALSE.
            nchain=1
            count_reject=0
            
            do while ((.not.rottest).and.(nchain<=maxattempts))  
                call rotation(chain,chains_rot,nseg,rottest)
                nchain=nchain+1
            enddo
        
            numconf=numconf+nchain ! number of attempted/rotated conformations

            if(rottest) then 
                conf=conf+1  
                do s=1,nseg
                    z = chains_rot(1,s+1)
                    x = chains_rot(2,s+1)
                    y = chains_rot(3,s+1)
                    r = sqrt(x**2 + (y+radius)**2+ z**2)
                    indexchain(conf,s)= int((r-radius)/delta)+1
                enddo
            else
               count_reject=count_reject+1
            endif
        enddo                  ! end rot loop  
      
        if(count_reject==rotmax) then
            conffile_reject=conffile_reject+1
        !    print*,"conformation ",conffile," rejected"
        endif

    enddo                     ! end while loop                
  
    !     .. end chains generation 
  
    reject_rate=conf*1.0d0/numconf

    if(conf.le.cuantas) then
        print*,"Something went wrong" 
        print*,"rotmax   = ",rotmax," attempts   = ",maxattempts
        print*,"conf     = ",conf," cuantas      = ",cuantas
        print*,"conffile = ",conffile," cuantasfile = ",cuantasfile
        print*,"rejection rate = ",reject_rate
        print*,"conffile reject= ",conffile_reject
        stop
    else
        print*,"Chains generated"
        print*,"rotmax   = ",rotmax," attempts   = ",maxattempts
        print*,"conf     = ",conf," cuantas      = ",cuantas
        print*,"conffile = ",conffile," cuantasfile = ",cuantasfile
        print*,"rejection rate = ",reject_rate
        print*,"conffile reject= ",conffile_reject
 
    endif
    
    deallocate(chain)
    deallocate(chains_rot)  
    deallocate(type_of_monomer_file)

    close(un)
  
end subroutine make_chains_file

! ismonomer_of_type is a table which row index is the segment  and column index correspond to the segment type  
! pre: type_of_monomer needs to be initialized
! post: table list of logicals indication is segement s is of type t

subroutine make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)
    
    implicit none

    logical, intent(out) :: ismonomer_of_type(:,:)
    integer, intent(in) :: type_of_monomer(:) 
    integer, intent(in) :: nseg
    integer, intent(in) :: nsegtypes

    ! local variable 
    integer :: s,t 

    do s=1,nseg 
        do t=1,nsegtypes   
            ismonomer_of_type(s,t)=.false.
        enddo
        ismonomer_of_type(s,type_of_monomer(s))=.true.    
    enddo

end subroutine make_type_table


! routine determines is segment type t is chargeable
! pre: zpol needs to initialized
! post: ismonomer_chargeable list of logicals

subroutine make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
    
    implicit none
 
    logical, intent(out) :: ismonomer_chargeable(:)
    integer, intent(in) :: zpol(:,:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t 

    do t=1,nsegtypes
        if(zpol(t,1)==0.and.zpol(t,2)==0) then 
            ismonomer_chargeable(t)=.false.
        else
            ismonomer_chargeable(t)=.true.
        endif
    enddo

end subroutine make_charge_table


logical function  is_polymer_neutral(ismonomer_chargeable, nsegtypes)
    
    implicit none
 
    logical, intent(in) :: ismonomer_chargeable(:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t,flag 

    flag=0
    do t=1,nsegtypes
        if(.not.ismonomer_chargeable(t)) flag=flag+1
    enddo

    if(flag==nsegtypes) then 
        is_polymer_neutral=.true. 
    else
        is_polymer_neutral=.false. 
    endif


end  function is_polymer_neutral
