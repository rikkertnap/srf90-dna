!     .. module file of chains variables

module chains
  
    !     .. variables

    implicit none

    integer(2), dimension(:,:), allocatable :: indexchain ! indexchain(conf,s) lattice coordinate of segment s of conformation conf
    integer, dimension(:), allocatable :: type_of_monomer               ! type of monomer represented as a number
    character(len=2), dimension(:), allocatable :: type_of_monomer_char ! type of monomer represented as a letter
    logical, dimension(:,:), allocatable :: ismonomer_of_type   ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable :: ismonomer_chargeable  ! ismonomer_chargeabl(s)=true if segment number "s" is acid or base  

    contains
  
    subroutine allocate_chains(cuantas,nseg,nsegtypes)
        implicit none

        integer, intent(in) :: cuantas,nseg,nsegtypes
        !     .. extra 100 in index because of  nchain rotations

        allocate(indexchain(cuantas+100,nseg))
        allocate(type_of_monomer(nseg)) 
        allocate(type_of_monomer_char(nseg))
        allocate(ismonomer_of_type(nseg,nsegtypes)) 
        allocate(ismonomer_chargeable(nsegtypes))
     
    end subroutine allocate_chains

end module chains
