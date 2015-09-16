!     .. module file of chains variables
module chains
  
    !     .. variables

    use globals
    implicit none

    integer(2), dimension(:,:), allocatable :: indexchain ! indexchain(conf,s) lattice coordinate of segment s of conformation conf
    integer, dimension(:), allocatable :: type_of_monomer 
    logical, dimension(:,:), allocatable :: ismonomer_of_type
    logical, dimension(:), allocatable :: ismonomer_chargeable

    contains
  
    subroutine allocate_chains(cuantas,nseg)
        implicit none

        integer, intent(in) :: cuantas,nseg

        !     .. extra 100 in index because of  nchain rotations

        allocate(indexchain(cuantas+100,nseg))
        allocate(type_of_monomer(nseg))
        allocate(ismonomer_of_type(nseg,nsegtypes))
        allocate(ismonomer_chargeable(nsegtypes))
     
    end subroutine allocate_chains

end module chains
