module parameters

    use molecules
    use physconst
    use mathconst
    use random
    use volume
  
    implicit none
  

    !    .. list of parameters

    !    .. variables 

    type(moleclist) :: xbulk      ! bulk volume fraction
    type(moleclist) :: expmu
    type(moleclist) :: v          ! not yet used !!!
    type(moleclist) :: z 
    type(moleclist) :: R 
    
    !    .. volume
    real(dp), dimension(:), allocatable :: vpol  ! volume of polymer segment of given type, vpol in units of vsol
    
    real(dp) :: vsol               ! volume of solvent in nm^3     
    real(dp) :: vNa                ! volume positive ion in units of vsol
    real(dp) :: vK                 ! volume positive ion in units of vsol
    real(dp) :: vCl                ! volume negative ion in units of vsol   
    real(dp) :: vCa                ! volume positive divalent ion in units of vsol
    real(dp) :: vNaCl
    real(dp) :: vKCl

    !     .. radii
  
    real(dp) :: RNa
    real(dp) :: RK
    real(dp) :: RCl
    real(dp) :: RCa
    real(dp) :: RNaCl
    real(dp) :: RKCl
    
    ! ..charge 

    integer, dimension(:,:), allocatable :: zpol          ! valence charge polymer
   
    integer :: zNa               ! valence charge positive ion 
    integer :: zK                ! valence charge positive ion 
    integer :: zCa               ! valence charge divalent positive ion 
    integer :: zCl               ! valence charge negative ion 
  
    !  .. input filenames select if chainmethod==file
    integer, parameter :: lenfname=40
    character(len=lenfname) :: chainfname,vpolfname,pKafname,typesfname

    ! .. other physical quanties

    real(dp) :: lseg             ! segment length only relevant if chaimethod == mc
    real(dp) :: Temperature      ! temperature in K
    real(dp) :: dielectW         ! dielectric constant of water 
    real(dp) :: lb               ! Bjerrum length	   
    real(dp) :: constqW          ! constant in Poisson eq dielectric constant of water 

    real(dp) :: sigma            ! sigma polymer coated on planar surface
    real(dp) :: sigmamin         ! minimal surface coverage
    real(dp) :: sigmamax         ! maximum surface coverage 
    real(dp) :: sigmastepsize    ! stepsize for increase/decrease of surface coverage    
    real(dp) :: sigmadelta       ! tolerance: sigmastepsize will not be smaller then sigmadelta
    integer  :: num_sigma        ! number of different sigmas
    real(dp), dimension(:), allocatable :: sigma_array ! array of surface coverages
    
    real(dp) :: sigmaSurf        ! surface density of acid on surface in nm^2
    real(dp) :: sigmaqSurf       ! surface charge density on surface in nm^2
    real(dp) :: psiSurf          ! surface potential	 
  
  
    integer :: itmax             ! maximum number of iterations
    real(dp) :: error            ! error imposed accuaracy
    real(dp) :: fnorm            ! L2 norm of residual vector function fcn  
    integer :: infile            ! infile=1 read input files infile!=1 no input files 
    integer :: iter              ! counts number of iterations
  
    character(len=8) :: method   ! method="kinsol" 
    character(len=8) :: chainmethod     ! method of generating chains ="MC" or "FILE"   
    integer :: readinchains      ! number of used/readin chains
    real(dp) :: spacer           ! spacer to anchor DNA chain to surface NP
    
  
    real(dp) :: avheightpol        ! average height of layer
    real(dp) :: avqpol             ! charge polymer of layer 
    real(dp), dimension(:), allocatable  :: avfdis   ! average degree of dissociation 
   
  
    !     .. weak polyelectrolyte variables 
    !     .. equibrium constant
  
    real(dp), dimension(:), allocatable :: K0a              ! intrinsic equilibruim constant
    real(dp), dimension(:), allocatable :: Ka               ! experimemtal equilibruim constant 
    real(dp), dimension(:), allocatable :: pKa              ! experimental equilibruim constant pKa= -log[Ka]
   
    real(dp) :: pKw                 ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-] 
  
    real(dp) :: K0ionNa             ! intrinsic equilibruim constant
    real(dp) :: KionNa              ! experimemtal equilibruim constant 
    real(dp) :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	  
    real(dp) :: K0ionK              ! intrinsic equilibruim constant
    real(dp) :: KionK               ! experimemtal equilibruim constant 
    real(dp) :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 
  
    !     .. bulk volume fractions 
  
    real(dp) :: pibulk             ! -ln(xsolbulk)
    real(dp) :: xNaClsalt          ! volume fraction of salt in bulk
    real(dp) :: xKClsalt           ! volume fraction of salt in bulk
    real(dp) :: xCaCl2salt         ! volume fraction of divalent salt in bulk
  
    real(dp) :: cHplus             ! concentration of H+ in bulk in mol/liter
    real(dp) :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real(dp) :: cNaCl              ! concentration of salt in bulk in mol/liter
    real(dp) :: cKCl               ! concentration of salt in bulk in mol/liter
    real(dp) :: cCaCl2             ! concentration of divalent salt in bulk in mol/liter
    real(dp) :: pHbulk             ! pH of bulk pH = -log([H+])
    real(dp) :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
   
    integer  :: num_cNaCl          ! number of different NaCl concentrations
    real(dp), dimension(:), allocatable :: cNaCl_array ! array of salt concentrations
        
    contains

        !   determine total number of nonlinear equations

        subroutine set_size_neq()

            use globals

            implicit none

            if (sysflag=="elect") then
                neq = 2 * nsize
            else if (sysflag=="bulk water") then
                neq = 5 
            else 
                print*,"Wrong value sysflag:  ",sysflag
                stop
            endif
                
        
        end subroutine set_size_neq


        !     purpose: initialize all constants parameters 
        !     pre: first read_inputfile has to be called   

        
        subroutine init_constants()

            use globals
            use volume
            use random
            use physconst
            
            implicit none      
            
         
            !  .. initializations of variables
     
            pi=dacos(-1.0d0)          ! pi = arccos(-1)
            itmax=2000                ! maximum number of iterations
            nr=nsize                  ! size of lattice in r-direction 
            
            !     .. charges
            
            zNa   = 1                 ! valence positive charged ion
            zK    = 1                 ! valence positive charged ion
            zCa   = 2                 ! valence divalent positive charged ion
            zCl   =-1                 ! valence negative charged ion
            
              
            !     .. radii
            
            RNa = 0.102d0             ! radius of Na+ in nm
            RK  = 0.138d0             ! radius of K+ in nm
            RCl = 0.181d0             ! radius of Cl- in nm
            RCa = 0.106d0             ! radius of Ca2+ in nm
            RNaCl = 0.26d0            ! radius of ion pair: this value is strange and not used !!
            RKCl  = 0.26d0            ! radius of ion pair: not used
            
            !     .. volume

            vsol = 0.030d0              ! volume water solvent molecule in (nm)^3 
            vNa  = ((4.0d0/3.0d0)*pi*(RNa)**3)/vsol 
            vK   = ((4.0d0/3.0d0)*pi*(RK)**3)/vsol 
            vCl  = ((4.0d0/3.0d0)*pi*(RCl)**3)/vsol 
            vCa  = ((4.0d0/3.0d0)*pi*(RCa)**3)/vsol 
            vNaCl = (vNa+vCl)          ! contact ion pair
            vKCl = (vK+vCl)           ! contact ion pair
           
            !   physical variables

            pKw=14.0d0                ! water equilibruim constant
            Temperature=298.0d0                 ! temperature in Kelvin
            dielectW=78.54d0          ! dielectric constant water
            lb=(elemcharge**2)/(4.0d0*pi*dielectW*dielect0*kBoltzmann*Temperature) ! bjerrum length in water=solvent in m
            lb=lb/1.0d-9              ! bjerrum length in water in nm
            seed  = 435672            ! seed for random number generator
            constqW = delta*delta*(4.0d0*pi*lb)/vsol ! multiplicative constant Poisson Eq. 
            
            !  .. initializations of input dependent variables 
           
            !  .. sigma's scaled by delta
            sigma = sigma * (1.0d0/(delta))    ! sigma scaled by delta observe no vpol*vsol term  
            sigmamin = sigmamin * (1.0d0/(delta))  
            sigmamax = sigmamax * (1.0d0/(delta))
            sigmastepsize = sigmastepsize * (1.0d0/(delta))
            sigmadelta = sigmadelta * (1.0d0/(delta))
            sigmaSurf = sigmaSurf * (4.0d0*pi*lb)*delta ! dimensionless surface charge 
        
        end subroutine init_constants
 
        !     purpose: initialize all constant parameters realated to the chains
        !     pre: first read_inputfile has to be called   

        
        subroutine init_chain_parameters
            
            use globals, only : nsegtypes,nseg
            use chains, only :  type_of_monomer,ismonomer_chargeable


            implicit none  

            ! fortran90 routines need either an interface or a module !!!
            interface 
                subroutine make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
                    implicit none
                    logical, intent(out) :: ismonomer_chargeable(:)
                    integer, intent(in) :: zpol(:,:)
                    integer, intent(in) :: nsegtypes
                end subroutine make_charge_table
            end interface

            ! segment length in nm
            if(chainmethod=="MC") lseg=0.545d0
            
            !  allocate array depending on nsegtypes

            !  volume polymer segments, all volume scaled by vsol
            allocate(vpol(nsegtypes))
            
            !  equilibrium constants 
            allocate(pKa(nsegtypes)) 
            allocate(Ka(nsegtypes))   
            allocate(K0a(nsegtypes))  

            !  charge of segment of given type
            allocate(zpol(nsegtypes,2)) 
        
            !  .. assign value to vpol  
            call read_volume(vpol,vsol,vpolfname, nsegtypes)

            !  .. assign value to pKa and zpol
            call read_pKas_and_zpol(pKa,zpol,pKafname, nsegtypes)

            !  .. assign value to type_of_monomer
            call read_type_of_monomer(type_of_monomer,typesfname, nseg)
         
            !  .. assign value to ismonomer_chargeable 
            call make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
            
        end subroutine init_chain_parameters
   
   
        !     purpose: initialize expmu needed by fcn 
        !     pre: first read_inputfile has to be called

        subroutine init_expmu()
     
            use globals
            use physconst
            
            implicit none 
            
            !     .. local variable
            
            real(dp),  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
            real(dp),  dimension(:), allocatable :: xguess  
            
            integer :: t
            
            allocate(x(5))
            allocate(xguess(5))
            
            !     .. initializations of input dependent variables, electrostatic part 
            
            cHplus = (10d0)**(-pHbulk) ! concentration H+ in bulk
            pOHbulk = pKw -pHbulk       
            cOHmin  = (10d0)**(-pOHbulk) ! concentration OH- in bulk
            
            xbulk%Hplus = (cHplus*Na/(1.0d24))*(vsol) ! volume fraction H+ in bulk vH+=vsol
            xbulk%OHmin = (cOHmin*Na/(1.0d24))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
            
            xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
            
            if(pHbulk.le.7) then      ! pH<= 7
                xbulk%Na=xNaClsalt*vNa/(vNa+vCl)  
                xbulk%Cl=xNaClsalt*vCl/(vNa+vCl) +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
            else                      ! pH >7
                xbulk%Na=xNaClsalt*vNa/(vNa+vCl) +(xbulk%OHmin -xbulk%Hplus)*vNa ! NaCl+ NaOH  
                xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)  
            endif
            
            xKClsalt = (cKCl*Na/(1.0d24))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
            xbulk%K = xKClsalt*vK/(vK+vCl)  
            xbulk%Cl = xbulk%Cl+xKClsalt*vCl/(vK+vCl)  
            
            xCaCl2salt = (cCaCl2*Na/(1.0d24))*((vCa+2.0d0*vCl)*vsol) ! volume fraction CaCl2 in mol/l
            xbulk%Ca=xCaCl2salt*vCa/(vCa+2.0d0*vCl)
            xbulk%Cl=xbulk%Cl+ xCaCl2salt*2.0d0*vCl/(vCa+2.0d0*vCl)
            
            xbulk%NaCl=0.0d0    ! no ion pairing
            xbulk%KCl=0.0d0     ! no ion pairing
            
            xbulk%sol=1.0d0 -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca   
            
            !     .. if Kion neq 0 ion pairing !
            
            !     .. intrinstic equilibruim constant acid        
            !     Kion  = 0.246d0 ! unit 1/M= liter per mol !!!
            K0ionK  = KionK /(vsol*Na/1.0d24) ! intrinstic equilibruim constant 
            K0ionNa = KionNa/(vsol*Na/1.0d24) ! intrinstic equilibruim constant 
            
            if((KionNa.ne.0.0d0).or.(KionK.ne.0.0d0)) then  

                sysflag="bulk water"        ! set solver to fcnbulk
                call set_size_neq()              ! number of nonlinear equations
                
                x(1)=xbulk%Na
                x(2)=xbulk%Cl
                x(3)=xbulk%NaCl
                x(4)=xbulk%K
                x(5)=xbulk%KCl
                
                xguess(1)=x(1)
                xguess(2)=x(2)
                xguess(3)=x(3)
                xguess(4)=x(4)
                xguess(5)=x(5)
               
                call solver(x, xguess, error, fnorm) 
                
                !     .. return solution
                
                xbulk%Na  =x(1)
                xbulk%Cl  =x(2)
                xbulk%NaCl=x(3)
                xbulk%K   =x(4)
                xbulk%KCl =x(5)

                ! reset of flags
                iter=0
                sysflag="elect"             ! set solver to fcnelect
                call set_size_neq()         ! number of non-linear  equation        
                
                xbulk%sol=1.0d0-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca 
                
            endif
             
            !     .. intrinstic equilibruim constants      
            do t=1,nsegtypes
                Ka(t)  = 10d0**(-pKa(t)) ! experimental equilibruim constant acid 
                K0a(t) = (Ka(t)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant 
            enddo
           
              
              
            !     !! check the unit is is vNa or vNa*vsol 
              
            pibulk = -dlog(xbulk%sol)  ! pressure (pi) of bulk
            
            ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 
            expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
            expmu%K     = xbulk%K    /(xbulk%sol**vK)
            expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
            expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
            expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
            expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
            expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
            expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
              
            !     .. end init electrostatic part 

            deallocate(x)
            deallocate(xguess)
            
        end subroutine init_expmu


        !  .. assign vpol from values in file named filename
        !  .. values vpol are normalized by vsol
        !  .. checked if file exists- content and length not checked     

        subroutine read_volume(vpol,vsol,filename, ntypes)

            use  myutils

            implicit none 
            
            !     .. arguments 
            real(dp), intent(inout)  :: vpol(:) 
            real(dp), intent(in) :: vsol 
            character(40), intent(in) :: filename  
            integer,  intent(in) :: ntypes 

            !      .. local variables
            integer :: ios, un  ! un = unit number
            integer :: t 
            character(80) :: istr,str

            !     .. reading in of variables from file
            open(unit=newunit(un),file=filename,iostat=ios,status='old')
            if(ios/=0 ) then
                write(istr,'(I2)')ios
                str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
                print*,str
                stop
            endif
        
            do t=1,ntypes   
                read(un,*)vpol(t)
                vpol(t)=vpol(t)/vsol
            enddo    

            close(un)
       
        end subroutine read_volume

        
        !  .. assign  pKA and zpol from values in file named filename
       
        subroutine read_pKas_and_zpol(pKa,zpol,filename, ntypes)
       
            use  myutils
            implicit none 
            
            !     .. arguments 
            real(dp), intent(inout) :: pKa(:)
            integer, intent(inout) ::  zpol(:,:)
            character(40), intent(in) :: filename
            integer,  intent(in) :: ntypes

            !      .. local variables
            integer :: ios, un  ! un = unit number
            integer :: t 
            character(80) :: istr,str

            !     .. reading in of variables from file
            open(unit=newunit(un),file=filename,iostat=ios,status='old')
            if(ios/=0 ) then
                write(istr,'(I2)')ios
                str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
                print*,str
                stop
            endif
            
            do t=1,ntypes   
                read(un,*)pKa(t),zpol(t,1),zpol(t,2)
            enddo    
            
            close(un)


        end subroutine read_pKas_and_zpol 


        !  .. assign  pKA and zpol from values in file named filename
       
        subroutine read_type_of_monomer(type_of_monomer,filename, nseg)
       
            use  myutils
            implicit none 
            
            !     .. arguments 
            integer, intent(inout) :: type_of_monomer(:)
            character(40), intent(in) :: filename
            integer,  intent(in) :: nseg

            !      .. local variables
            integer :: ios, un  ! un = unit number
            integer :: s
            character(80) :: istr,str,letter

            !     .. reading in of variables from file
            open(unit=newunit(un),file=filename,iostat=ios,status='old')
            if(ios/=0 ) then
                write(istr,'(I2)')ios
                str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
                print*,str
                stop
            endif
            
            do s=1,nseg   
                read(un,*)type_of_monomer(s),letter
!                print*,type_of_monomer(s),letter
            enddo    
            
            close(un)


        end subroutine read_type_of_monomer



 end module parameters
