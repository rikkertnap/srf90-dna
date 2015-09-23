! --------------------------------------------------------------|
! myio.f90:                                                     |
! --------------------------------------------------------------|

module myio

    use precision_definition
    implicit none

    character(len=10) :: LogName        ! filename of status log file
    integer, parameter :: LogUnit=321   ! associated unit number
    real(dp), dimension(:), allocatable  :: cs          ! salt concentrations     
   
    contains    


    subroutine read_inputfile

        use globals  
        use parameters
        use volume  ! include in order to assign nradius
        use myutils, only : newunit
      
        implicit none
      
        character(len=8) :: fname
        integer :: ios, un_input  ! un = unit number
      
        !     .. reading in of variables from file
           
        write(fname,'(A8)')'input.in'
        open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
        if(ios > 0 ) then
            print*, 'Error opening file input.in : iostat =', ios
            stop
        endif
      
        read(un_input,*)method
        read(un_input,*)sysflag
        call check_value_sysflag()
        read(un_input,*)runflag
        call check_value_runflag()
        read(un_input,*)chainmethod
        if(chainmethod=="FILE") read(un_input,*)chainfname
        read(un_input,*)vpolfname
        read(un_input,*)pKafname
        read(un_input,*)typesfname
        if(chainmethod=="FILE") read(un_input,*)spacer
        if(chainmethod=="MC") read(un_input,*)lseg
        read(un_input,*)error             
        read(un_input,*)infile          ! guess  1==yes
        read(un_input,*)radius
        read(un_input,*)pHbulk
        read(un_input,*)KionNa
        read(un_input,*)KionK
        read(un_input,*)cKCl
        read(un_input,*)cCaCl2
        read(un_input,*)nsize
        read(un_input,*)nseg
        read(un_input,*)nsegtypes ! carefully need to be overwriten when chainmethod==FILE
        read(un_input,*)cuantas
        if(runflag=="range") then
            read(un_input,*)sigmamin
            read(un_input,*)sigmamax
            read(un_input,*)sigmastepsize
            read(un_input,*)sigmadelta
        endif 
        read(un_input,*)delta    ! delta 0.3/0.2 /0.4

        close(un_input)

        sigmaSurf=0.0d0   ! value not equal zero not yet properly workinf for free energy routine
        
        call check_input_value()
        
    end subroutine read_inputfile
 

    subroutine check_value_sysflag()

        use globals
        use parameters
            
        implicit none

        character(len=15) :: sysflagstr(2) 
        integer :: i
        logical :: flag

        ! permissible values of sysflag

        sysflagstr(1)="elect"
        sysflagstr(2)="electHP" ! not used yet 
        
        flag=.FALSE.

        do i=1,2
            if(sysflag==sysflagstr(i)) flag=.TRUE.
        enddo
        if (flag.eqv. .FALSE.) then 
            print*,"Error: value of sysflag is not permissible"
            print*,"sysflag = ",sysflag
            stop
        endif
        
    end subroutine
    
    subroutine check_value_runflag()

        use globals
        use parameters
            
        implicit none

        character(len=15) :: runflagstr(3) 
        integer :: i
        logical :: flag

        ! permissible values of runflag

        runflagstr(1)="default"
        runflagstr(2)="input"
        runflagstr(3)="range"  

        flag=.FALSE.

        do i=1,3
            if(runflag==runflagstr(i)) flag=.TRUE.
        enddo
        if (flag.eqv. .FALSE.) then 
            print*,"Error: value of runflag is not permissible"
            print*,"runflag = ",runflag
            stop
        endif

    end subroutine      


    ! check  input values on consistency with sysflag

    subroutine check_input_value()

        use globals
        use parameters

        implicit none


    end subroutine check_input_value

    ! purpose print input values for debugging only

    subroutine print_input_value()

        use globals
        use parameters

        implicit none

        print*,"nsegtypes=",nsegtypes    
        print*,"delta=",delta
        print*,"nseg=",nseg

    end subroutine print_input_value


    subroutine set_value_salt(runflag)

        use parameters, only : num_cNaCl,cNaCl_array
        use myutils, only : newunit

        implicit none

        character(len=8) :: runflag
        character(len=7) :: fname
        integer :: ios
        integer :: i    
        integer :: un_cs
        
        if(runflag=="default") then

            num_cNaCl=6  
            allocate(cNaCl_array(num_cNaCl))
            
            cNaCl_array(1)=0.25d0
            cNaCl_array(2)=0.20d0
            cNaCl_array(3)=0.15d0
            cNaCl_array(4)=0.10d0
            cNaCl_array(5)=0.05d0
            cNaCl_array(6)=0.01d0

        else if(runflag=="input".or.runflag=="range") then

            !     .. read salt concentrations from file
    
            write(fname,'(A7)')'salt.in'
            open(unit=newunit(un_cs),file=fname,iostat=ios,status='old')
            if(ios > 0 ) then
                print*, 'Error opening file salt.in : iostat =', ios
                stop
            endif

            read(un_cs,*)num_cNaCl ! read number of salt concentration form file
            allocate(cNaCl_array(num_cNaCl)) 
            
            do i=1,num_cNaCl     ! read value salt concentration
                read(un_cs,*)cNaCl_array(i)
            enddo    
            close(un_cs)
        endif 
        

    end subroutine  set_value_salt
    
    subroutine set_value_sigma(runflag)

        use parameters, only : num_sigma,sigma_array
        use myutils, only : newunit

        implicit none

        character(len=8) :: runflag

        character(len=8) :: fname
        integer :: ios, un_sg
        integer :: i 

        if(runflag=="default") then
            
            num_sigma=9
            allocate(sigma_array(num_sigma))

            sigma_array(1)=0.001d0
            sigma_array(2)=0.005d0
            sigma_array(3)=0.010d0
            sigma_array(4)=0.025d0
            sigma_array(5)=0.050d0
            sigma_array(6)=0.075d0
            sigma_array(7)=0.100d0
            sigma_array(8)=0.125d0
            sigma_array(9)=0.150d0

        
        else if(runflag=="input") then 
            ! read values of sigma from file

            write(fname,'(A8)')'sigma.in'
            open(unit=newunit(un_sg),file=fname,iostat=ios,status='old')
            if(ios > 0 ) then
                print*, 'Error opening file sigma.in: iostat =', ios
                stop
            endif

            read(un_sg,*)num_sigma
            allocate(sigma_array(num_sigma))

            do i=1,num_sigma
                read(un_sg,*)sigma_array(i)
            enddo    
            close(un_sg)
            
        endif    

    end subroutine  set_value_sigma


    subroutine output_elect()

        !     .. variables and constant declaractions
        use globals 
        use volume
        use parameters
        use field
        use energy 
        use chains
        use myutils
        !     use endpoint

        implicit none

        !     .. output file names       

        character(len=90)  :: sysfilename     
        character(len=90)  :: xsolfilename 
        character(len=90)  :: xpolfilename 
        character(len=100) :: xNafilename
        character(len=100) :: xKfilename
        character(len=90)  :: xCafilename
        character(len=110) :: xNaClfilename
        character(len=90)  :: xKClfilename
        character(len=90)  :: xClfilename
        character(len=90)  :: potentialfilename
        character(len=90)  :: chargefilename
        character(len=90)  :: xHplusfilename
        character(len=90)  :: xOHminfilename
        character(len=100) :: densfdisfilename
        character(len=110) :: densfracionpairfilename

        character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   
        integer :: un_sys,un_xpol,un_xsol,un_xNa,un_xCl,un_xNaCl,un_ionpair,un_xK,un_xKCl
        integer :: un_xOHmin,un_xHplus,un_fdis,un_xCa
        integer :: un_psi,un_charge  ! unit numbers

        !     .. local arguments

        integer :: i,n,t           ! dummy indexes
        character(len=100) :: fnamelabel
        character(len=10) :: rstr
        logical :: isneutral 
        
        !     .. executable statements 
        ! fortran90 routines need either an interface or a module !!!
        interface
            logical function  is_polymer_neutral(ismonomer_chargeable, nsegtypes)
    
                implicit none
 
                logical, intent(in) :: ismonomer_chargeable(:)
                integer, intent(in) :: nsegtypes
            end  function is_polymer_neutral
                
        end interface
        
        n=nr

        fmt2reals = "(2ES25.16E3)"  
        fmt3reals = "(3ES25.16E3)"  
        fmt4reals = "(4ES25.16E3)"  
        fmt5reals = "(5ES25.16E3)" 
        fmt6reals = "(6ES25.16E3)" 

        print*,"output-> sigma=",sigma
        !     .. make lable filenames 
        write(rstr,'(F5.3)')sigma*delta
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        if(cCaCl2.ne.0) then 
            write(rstr,'(F5.3)')cCaCl2
            fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr)) 
        endif    
        !     .. make filenames 
    
        sysfilename='system.'//trim(fnamelabel)//'.dat'
        xpolfilename='xpol.'//trim(fnamelabel)//'.dat'
        xsolfilename='xsol.'//trim(fnamelabel)//'.dat'
        xNafilename='xNaions.'//trim(fnamelabel)//'.dat'
        xKfilename='xKions.'//trim(fnamelabel)//'.dat'
        xCafilename='xCaions.'//trim(fnamelabel)//'.dat'
        xNaClfilename='xNaClionpair.'//trim(fnamelabel)//'.dat'
        xKClfilename='xKClionpair.'//trim(fnamelabel)//'.dat'
        xClfilename='xClions.'//trim(fnamelabel)//'.dat'
        potentialfilename='potential.'//trim(fnamelabel)//'.dat'
        chargefilename='charge.'//trim(fnamelabel)//'.dat'
        xHplusfilename='xHplus.'//trim(fnamelabel)//'.dat'
        xOHminfilename='xOHmin.'//trim(fnamelabel)//'.dat'
        densfdisfilename='densityfrac.'//trim(fnamelabel)//'.dat'
        densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)//'.dat'
          
        !     .. opening files        
               
        open(unit=newunit(un_sys),file=sysfilename)   
        open(unit=newunit(un_xpol),file=xpolfilename)
        open(unit=newunit(un_xsol),file=xsolfilename)
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_psi),file=potentialfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
        open(unit=newunit(un_fdis),file=densfdisfilename)

        if(cKCL/=0.0d0) open(unit=newunit(un_xK),file=xKfilename)
        if(cCaCl2/=0.0d0) open(unit=newunit(un_xCa),file=xCafilename)
        if(KionNa/=0.0d0) open(unit=newunit(un_xNaCl),file=xNaClfilename)
        if(KionK/=0.0d0) open(unit=newunit(un_xKCl),file=xKClfilename)
        if(KionNa/=0.0d0 .or.KionK/=0.0d0 ) open(unit=newunit(un_ionpair),file=densfracionpairfilename)
        
        
        write(un_psi,*)radius,psiSurf

        do i=1,n
         
            write(un_xpol,*)rc(i),xpol(i)
            write(un_xsol,*)rc(i),xsol(i)
            write(un_xNa,*)rc(i),xNa(i)
            write(un_xCl,*)rc(i),xCl(i)
            write(un_psi,*)rc(i),psi(i)
            write(un_charge,*)rc(i),rhoq(i)
            write(un_xHplus,*)rc(i),xHplus(i)
            write(un_xOHmin,*)rc(i),xOHmin(i)
            write(un_fdis,*)rc(i),(fdis(i,t),t=1,nsegtypes)

        enddo         
            
        if(cKCl/=0.0d0) then
            do i=1,n 
                write(un_xK,*)rc(i),xK(i)
            enddo
        endif    
        if(cCaCl2/=0.0d0) then
            do i=1,n
                write(un_xCA,*)rc(i),xCa(i)
            enddo
        endif
        if(KionNa/=0.0d0) then
            do i=1,n 
                write(un_xNaCl,*)rc(i),xNaCl(i)
            enddo
        endif
        if(KionK/=0.0d0.and.cKCl/=0.0d0) then    
            do i=1,n 
                write(un_xKCl,*)rc(i),xKCl(i)  
            enddo
        endif    
        if(KionNa/=0.0d0.and.KionK/=0.0d0) then  
            do i=1,n 
                write(un_ionpair,*)rc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            enddo
        endif   



        isneutral=is_polymer_neutral(ismonomer_chargeable, nsegtypes)


        !     .. system information 


        write(un_sys,*)'system      = spherical weakpolyelectrolyte brush'
        write(un_sys,*)'version     = ',VERSION
        write(un_sys,*)'free energy = ',FE
        write(un_sys,*)'energy bulk = ',FEbulk 
        write(un_sys,*)'deltafenergy = ',deltaFE
        if(isneutral) then 
            write(un_sys,*)'FEalt       = ',FEalt
            write(un_sys,*)'FEbulkalt   = ',FEbulkalt
            write(un_sys,*)'deltaFEalt  = ',deltaFEalt
        else
            write(un_sys,*)'FEalt       = not yet correct'
        endif     
        write(un_sys,*)'FEconf      = ',FEconf
        write(un_sys,*)'FEtrans%sol = ',FEtrans%sol  
        write(un_sys,*)'fnorm       = ',fnorm
        write(un_sys,*)'q residual  = ',qres
        write(un_sys,*)'error       = ',error
        write(un_sys,*)'sumphi      = ',sumphi 
        write(un_sys,*)'check phi   = ',sumphi-sigma*delta*nseg
        write(un_sys,*)'FEq         = ',FEq 
        write(un_sys,*)'FEpi        = ',FEpi
        write(un_sys,*)'FErho       = ',FErho
        write(un_sys,*)'FEel        = ',FEel
        write(un_sys,*)'FEelsurf    = ',FEelsurf 
        write(un_sys,*)'FEbind      = ',FEbind
        write(un_sys,*)'q           = ',q
        write(un_sys,*)'mu          = ',-dlog(q)
        write(un_sys,*)'nseg        = ',nseg
        write(un_sys,*)'nsegtypes   = ',nsegtypes
        write(un_sys,*)'nr          = ',nr
        write(un_sys,*)'delta       = ',delta 
        write(un_sys,*)'radius      = ',radius
        write(un_sys,*)'vsol        = ',vsol
        write(un_sys,*)'vNa         = ',vNa*vsol
        write(un_sys,*)'vCl         = ',vCl*vsol
        write(un_sys,*)'vCa         = ',vCa*vsol
        write(un_sys,*)'vK          = ',vK*vsol
        write(un_sys,*)'vNaCl       = ',vNaCl*vsol
        write(un_sys,*)'vKCl        = ',vKCl*vsol
        write(un_sys,*)'cNaCl       = ',cNaCl
        write(un_sys,*)'cKCl        = ',cKCl
        write(un_sys,*)'cCaCl2      = ',cCaCl2
        write(un_sys,*)'pHbulk      = ',pHbulk
        write(un_sys,*)'KionNa      = ',KionNa
        write(un_sys,*)'KionK       = ',KionK
        write(un_sys,*)'K0ionNa     = ',K0ionNa
        write(un_sys,*)'K0ionK      = ',K0ionK
        write(un_sys,*)'xNabulk     = ',xbulk%Na
        write(un_sys,*)'xClbulk     = ',xbulk%Cl
        write(un_sys,*)'xKbulk      = ',xbulk%K
        write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
        write(un_sys,*)'xKClbulk    = ',xbulk%KCl
        write(un_sys,*)'xCabulk     = ',xbulk%Ca
        write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
        write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
        write(un_sys,*)'sigma       = ',sigma*delta
        write(un_sys,*)'psiSurf     = ',psiSurf
        write(un_sys,*)'avheightpol = ',avheightpol
        write(un_sys,*)'avqpol      = ',avqpol
        write(un_sys,*)'dielectW    = ',dielectW
        write(un_sys,*)'lb          = ',lb
        write(un_sys,*)'T           = ',Temperature        
        write(un_sys,*)'zNa         = ',zNa
        write(un_sys,*)'zCa         = ',zCa
        write(un_sys,*)'zK          = ',zK
        write(un_sys,*)'zCl         = ',zCl
        write(un_sys,*)'nsize       = ',nsize  
        write(un_sys,*)'cuantas     = ',cuantas
        write(un_sys,*)'nseg        = ',nseg
        write(un_sys,*)'nsegtypes   = ',nsegtypes
        write(un_sys,*)'iterations  = ',iter
        write(un_sys,*)'chainmethod = ',chainmethod 
        if(chainmethod.eq."FILE") then
            write(un_sys,*)'readinchains = ',readinchains
            write(un_sys,*)'spacer  = ',spacer
        endif
        write(un_sys,*)'sysflag     = ',sysflag
        write(un_sys,*)'runflag     = ',runflag

        close(un_sys)
        close(un_xpol)
        close(un_xsol)
        close(un_xNa)
        close(un_xCl)
        close(un_psi)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
        close(un_ionpair)
        close(un_fdis)


        close(un_xCa)
        close(un_xK)
        close(un_xNaCl)

    end subroutine output_elect



    subroutine output()

        use globals, only : sysflag
        implicit none

        if(sysflag=="elect") call output_elect()

    end subroutine output

     
end module myio
