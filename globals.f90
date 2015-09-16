!     .. module file of global variables 

module globals
  
    use precision_definition
    use mathconst
    use mathconst

    implicit none

    !     .. variables

    integer  :: nsize         ! size lattice, numer of layers
    integer(8) :: neq         ! number of non-linear equations
    integer  :: nseg          ! length of polymer
    integer  :: nsegtypes     ! number of segment types 
    integer  :: cuantas       ! number of polymer configurations
    character(len=15) :: sysflag       ! sysflag selects fcn    
    character(len=2) :: bcflag(2)      ! bcflag selects bc of surface

end module globals

