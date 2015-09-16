module initxvector

    implicit none
  
    
    contains

    !     purpose: initalize x and xguess

    subroutine init_guess_elect(x, xguess)
          
        use globals
        use parameters
        use field
        use volume
        use myutils, only : newunit
        
        implicit none
      
        real(dp) :: x(neq)       ! iteration vector 
        real(dp) :: xguess(neq)  ! guess vector
      
        !     ..local variables 
        integer :: i
        character(len=8) :: fname(2)
        integer :: ios,un_file(2)
      
      
        ! .. init guess all xbulk     

        do i=1,nr
            x(i)=xbulk%sol
            x(i+nr)=0.000d0
        enddo
      
        if (infile.eq.1) then   ! infile is read in from file/stdio  
        
            write(fname(1),'(A7)')'xsol.in'
            write(fname(2),'(A6)')'psi.in'
                  
            do i=1,2 ! loop files
                open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
                if(ios >0 ) then    
                    print*, 'file unit num  =',un_file(i),' file name =',fname(i)
                    print*, 'Error opening file : iostat =', ios
                    stop
                endif
            enddo
         
            read(un_file(2),*),psisurf          ! electrostatix potential
            do i=1,nr
                read(un_file(1),*),xsol(i)      ! solvent
                read(un_file(2),*),psi(i)       ! degree of complexation A
                x(i)      = xsol(i)             ! placing xsol  in vector x
                x(i+nr)   = psi(i)               ! placing psi  in vector x
            enddo
        
            do i=1,2
                close(un_file(i))
            enddo

        endif
        !     .. end init from file 
      
        do i=1,neq
            xguess(i)=x(i)
        enddo

    end subroutine init_guess_elect


    subroutine init_guess_neutral(x, xguess)
          
        use globals
        use parameters
        use field
        use volume
        use myutils, only : newunit

        implicit none
      
        real(dp) :: x(neq)       ! volume fraction solvent iteration vector 
        real(dp) :: xguess(neq)  ! guess fraction  solvent 
      
        !     ..local variables 
        integer :: n, i
        character(len=10) :: fname
        integer :: ios,un_file
      
        !     .. init guess all xbulk      

        do i=1,nr
            x(i)=xbulk%sol
        enddo
      
      
        if (infile.eq.1) then   ! infile is read in from file/stdio  
            write(fname,'(A7)')'xsol.in'
         
            open(unit=newunit(un_file),file=fname,iostat=ios,status='old')
            if(ios >0 ) then
                print*, 'file number =',un_file,' file name =',fname
                print*, 'Error opening file : iostat =', ios
                stop
            endif
            
            do i=1,nr
                read(un_file,*),xsol(i)       ! solvent
                x(i)      = xsol(i)       ! placing xsol in vector x
            enddo
            
            close(un_file)
            
        endif
        !     .. end init from file 
      
        do i=1,neq
            xguess(i)=x(i)
        enddo
    end subroutine init_guess_neutral


    !     .. copy solution of previous solution ( distance ) to create new guess
    !     .. data x=(pi,psi) and pi and psi order and split into  blocks
    !     .. of size (nptso,nptsi,nptsb,nptss) =( outside, inside, boundary, on sphere )

    subroutine make_guess(x, xguess,val,valfirst,flagstored,xstored)
      
        use globals
        use parameters
        use volume

        implicit none

        real(dp), intent(in) :: x(neq)         ! iteration vector 
        real(dp), intent(out) :: xguess(neq)    ! guess volume fraction solvent and potential 
        real(dp), intent(in) :: val,valfirst         ! current and first salt concentration
        logical, optional, intent(in) :: flagstored
        real(dp), optional, intent(out) :: xstored(neq)

        !     ..local variables 
        integer :: i
      
        if(present(flagstored)) then
            if(present(xstored)) then
                if(flagstored) then   
#ifdef DEBUG 
                        print*,"flagstored==true"   
#endif
                    ! copy solution xstored in xguess and x
                    do i=1,neq  
                        xguess(i)=xstored(i)
                    enddo
                else if(val==valfirst) then       ! first guess
#ifdef DEBUG 
                        print*,"hello val==valfirst"
#endif    
                    if(sysflag.eq."elect") call init_guess_elect(x,xguess)
                    if(sysflag.eq."neutral") call init_guess_neutral(x,xguess)
                else  
                    ! print*,"hello val/=valfirst"   
                    do i=1,neq
                        xguess(i)=x(i)      ! volume fraction solvent 
                   enddo
                endif
            else
                print*,"Error: argument xstored not present, while flagstored present"
                stop 
            endif 
        else if(val==valfirst) then       ! first guess
#ifdef DEBUG 
                print*,"hello val==valfirst"
#endif
            if(sysflag.eq."elect") call init_guess_elect(x,xguess)
            if(sysflag.eq."neutral") call init_guess_neutral(x,xguess)
        else  
            ! print*,"hello val/=valfirst"    
            do i=1,neq
                xguess(i)=x(i)      ! volume fraction solvent 
            enddo

        endif

    end subroutine make_guess

    
    


end module initxvector
