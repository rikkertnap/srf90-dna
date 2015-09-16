! --------------------------------------------------------------|
! Solves the SCMFT eqs for WEAK polyelectrolytes polymers       |
! coated onto a spherical surface,                              | 
! input: see myio.f90                                           | 
! --------------------------------------------------------------|
      
program brushweakpolyelectrolyte 
    
    !     .. variable and constant declaractions 
    use globals       ! parameters definitions 
    use mathconst
    use volume
    use random
    use field
    use parameters
    use matrices
    use energy
    use chains
    use listfcn 
    use fcnpointer
    use initxvector
    use myio
    
    implicit none  
    
    real(dp),  dimension(:), allocatable :: x         ! iteration vector 
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec     
  
    integer :: i,c,sg,s         ! dummy indices      
    logical :: use_xstored      
  
    !     .. executable statements 
    !     .. init 
  
    call read_inputfile()
    call init_constants()
    call init_matrices()        ! init matrices for chain generation
    call allocate_chains(cuantas,nseg)
    call make_chains(chainmethod) ! generate polymer configurations 
!    call make_sequence_chain(period,chaintype)
    call allocate_geometry(nsize)
    call make_geometry()        ! generate volume elements lattice 
    call allocate_field(nsize,nsegtypes) 
    call set_size_neq()         ! number of non-linear equation neq
  
    allocate(x(neq))
    allocate(xguess(neq))
    allocate(xstored(neq))
    allocate(fvec(neq))
  
    !     .. computation starts
    
    call set_value_salt(runflag)
    call set_value_sigma(runflag)
      
     
    if(runflag=='range') then 
        
        use_xstored=.false.

        do c=1,num_cNaCl                ! loop over salt concentration 
            cNaCl=cNaCl_array(c)
            call init_expmu()
            sigma=sigmamin
            do while (sigmamin<=sigma.and.sigma<=sigmamax.and.sigmastepsize>=sigmadelta) 
                call make_guess(x,xguess,sigma,sigmamin,use_xstored,xstored)
                call solver(x, xguess, error, fnorm) 
                if(isNaN(fnorm)) then 
                    print*,"no solution: backstep" 
                    sigmastepsize=sigmastepsize/2.0d0 ! smaller 
                    sigma=sigma-sigmastepsize ! step back
                    do i=1,neq
                        x(i)=xguess(i)
                    enddo       
                else 

                    do i=1,nr                
                        xsol(i)=x(i)        ! solvent volume fraction   
                        psi(i) =x(i+nr)     ! potential           
                    enddo
               
                    call fcnenergy()        ! free energy
                    call average_height()      
                    call average_charge_polymer()
                    call average_fdis_polymer()
                    call output()           ! writing of output
                    sigma=sigma+sigmastepsize

                endif

                iter  = 0                ! reset of iteration counter 
               
                if(sigma==sigmamin) then
                    do i=1,neq
                        xstored(i)=x(i)
                    enddo
                endif
                if(sigma==sigmamax) then 
                    use_xstored=.true.
                else
                    use_xstored=.false.
                endif     

            enddo

        enddo
    else

        use_xstored=.false.
        do c=1,num_cNaCl                ! loop over salt concentration 
            cNaCl=cNaCl_array(c)
            call init_expmu()
            do  sg=1,num_sigma               ! loop over sigma values
                sigma=sigma_array(sg)
                call make_guess(x,xguess,sigma,sigma_array(1),use_xstored,xstored)
                call solver(x, xguess, error, fnorm) 
                do i=1,nr                
                    xsol(i)=x(i)        ! solvent volume fraction   
                    psi(i) =x(i+nr)     ! potential           
                enddo
           
                call fcnenergy()        ! free energy
                call average_height()      
                call average_charge_polymer()
                call output()           ! writing of output
              
                iter  = 0                ! reset of iteration counter 
               
                if(sg==1) then
                    do i=1,neq
                        xstored(i)=x(i)
                    enddo
                endif
                if(sg==num_sigma) then 
                    use_xstored=.true.
                else
                    use_xstored=.false.
                endif     
            enddo
        enddo

    endif    
    
    deallocate(x)
    deallocate(xguess)
    deallocate(xstored)
    deallocate(fvec)
    if(runflag/="range") then
        deallocate(cs)  
        deallocate(sigma_array)     
    endif   
    call deallocate_field()

end program brushweakpolyelectrolyte
