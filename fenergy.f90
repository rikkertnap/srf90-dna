! --------------------------------------------------------------|
! april  2015                                                   | 
! fenergy.f:                                                    |
! calculates the free energy and bulk free energy               |
! --------------------------------------------------------------|


!     .. module file of energy variables

module energy 

    use precision_definition
    use molecules

    implicit none
  
    !     .. variables
  
    real(dp) :: FE                  ! free energy
    real(dp) :: FEbulk              ! free energybulk
    real(dp) :: deltaFE             ! free energy difference
  
    !     .. auxiliary variable used in free energy computation  
    real(dp) :: FEq                 ! partition function poly A and B 
    real(dp) :: FEpi                ! sum over pi
    real(dp) :: FErho               ! sum over densities
    real(dp) :: FEel                ! electrostatics energ
    real(dp) :: FEelsurf            ! electrostatics energy from  surface
    real(dp) :: FEchem              ! chemical energy
    real(dp) :: FEbind              ! complexation contribution
  
    real(dp) :: sumphi              ! check integrale over phiA
    real(dp) :: qres                ! charge charge
    
    real(dp) :: FEconf
    real(dp) :: FEtranssol

    real(dp) :: FEalt               ! free energy
    real(dp) :: FEbulkalt           ! free energybulk
    real(dp) :: deltaFEalt          ! free energy difference delteFE=FE-FEbulk

    type(moleclist) :: FEtrans,FEchempot,FEtransbulk,FEchempotbulk
    type(moleclist) :: deltaFEtrans,deltaFEchempot


    contains

    subroutine fcnenergy()
 
        use globals
        implicit none
    
        call fcnenergy_elect()
        call fcnenergy_elect_alternative()
   
    end subroutine fcnenergy

   
    subroutine fcnenergy_elect()
    
        use globals
        use volume
        use parameters
        use field

        implicit none

        !  .. local arguments 
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf              ! total charge on surface 
        real(dp), parameter :: sigmaTOL = 0.00000001d0             ! tolerance of surface coverage, below no polymers 
        real(dp) :: Asurf              ! area surface sphere
        integer :: i,t                   ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 


        !  .. computation of free energy 

        Asurf = 4.0d0*pi*radius*radius
    
        FEpi  = 0.0d0
        FErho = 0.0d0
        FEel  = 0.0d0
        FEelsurf = 0.0d0
        sumphi = 0.0d0
        
        FEq = 0.0d0
        FEbind = 0.0d0
        FEchem = 0.0d0
        qres = 0.0d0

        do i=1,nr
            FEpi = FEpi  + deltaG(i)*dlog(xsol(i))
            FErho = FErho - deltaG(i)*(xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel = FEel  - deltaG(i) * rhoq(i) * psi(i)/2.0d0        
           
            qres = qres + deltaG(i) * rhoq(i)
            do t=1,nsegtypes
                sumphi = sumphi + deltaG(i) * rhopol(i,t)
            enddo   
        enddo
    
        FEel  = (delta/vsol)*FEel
        FEpi  = (delta/vsol)*FEpi
        FErho = (delta/vsol)*FErho
        qres = (delta/vsol)*qres
        sumphi = delta*sumphi

        FEq = -delta*sigma*dlog(q)
    
        if(sigma <= sigmaTOL) then 
            FEq = 0.0d0 
        endif
    
        !     .. surface charge constribution  
    
        sigmaq0= sigmaqSurf/(delta*4.0d0*pi*lb) ! dimensional charge density  
        psi0 = psi(1)+sigmaqSurf/2.0d0 ! surface potential units !!
    
        FEelsurf = FEelsurf + sigmaq0 * psi0 
        FEelsurf = FEelsurf/2.0d0
    
    !    FEchem = (sigmaSurf/(4.0d0*pi*lb*delta))*dlog(1.0-fdisSurf) -2.0d0*FEelsurf
        FEchem  =0.0d0

        !     .. total free energy per area of surface 

        FE = FEq + FEpi + FErho + FEel + FEelsurf
        
        print*,"FE = " ,FE
    
        qsurf  = sigmaqSurf/(4.0d0*pi*lb*delta) ! total surface charge
  
        print*,"qsurf=",qsurf,"qres=",qres    

        qres = qres + qsurf  ! total residual charge 

        volumelat=(4.0d0/3.0d0)*pi*((nr*delta+radius)**3-(radius**3))   !volume lattice divide by area NP
        volumelat=volumelat/(4.0d0*pi*(radius**2))
        FEbulk   = dlog(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ & 
            xbulk%Na/vNa +xbulk%Ca/vCa +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_elect

    

    

    ! .. computes conformational entropy contrbution to free energy 


    subroutine FEconf_elect()

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none
        
        !     .. declare local variables

        real(dp) :: exppi(nsize,nsegtypes)   ! auxilairy variable for computing P(\alpha) 
        integer :: i,k,c,s,t         ! dummy indices
        real(dp) :: pro
        
        real(dp), parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 

        !     .. executable statements 
        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,nr
                    exppi(i,t)  = (xsol(i)**vpol(t))*exp(-zpol(t,2)*psi(i) )/fdis(i,t)   ! auxilary variable palpha
                enddo  
            else
                do i=1,nr
                    exppi(i,t)  = xsol(i)**vpol(t)
                enddo  
            endif   
        enddo   
     
        !   .. computation conformational entropy polymer 

        FEconf=0.0d0
        
        do c=1,cuantas               ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nseg              ! loop over segments 
                k=indexchain(c,s)
                t=type_of_monomer(s)
                pro = pro*exppi(k,t)
            enddo
            FEconf=FEconf+pro*dlog(pro)
        enddo

        ! normalize
        FEconf=(FEconf/q-dlog(q))*(sigma*delta)   
            
        
    end subroutine FEconf_elect


    ! .. older routine FEtrans_entroy more general

    real(dp) function FEtrans_fcn(xvol,vol)

        use globals
        use field
        use parameters
        implicit none

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: vol    

        integer :: i

        FEtrans_fcn=0.0d0

        do i=1,nr
            FEtrans_fcn=FEtrans_fcn + deltaG(i)*xvol(i)*(dlog(xvol(i))-1.0d0)
        enddo
        FEtrans_fcn= delta *FEtrans_fcn /vol
        return 

    end function FEtrans_fcn   
 
    real(dp) function FEtrans_entropy(xvol,xvolbulk,vol,flag)
    
        use globals
        use parameters
        implicit none

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: xvolbulk 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        integer :: i


        if(xvolbulk==0.0d0) then 
            FEtrans_entropy=0.0d0
        else
            FEtrans_entropy=0.0d0
            if(present(flag)) then
            ! water special case because vsol treated diffetent then vi  
                do i=1,nr
                    FEtrans_entropy=FEtrans_entropy + deltaG(i)*xvol(i)*(dlog(xvol(i))-1.0d0)
                enddo 
                FEtrans_entropy = delta*FEtrans_entropy/vol
            else 
                do i=1,nr
                    FEtrans_entropy = FEtrans_entropy +deltaG(i)*xvol(i)*(dlog(xvol(i)/vol)-1.0d0)
                enddo
                FEtrans_entropy = delta*FEtrans_entropy/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy   
    
 

    real(dp) function FEchem_pot(xvol,expchempot,vol,flag)
    
        use globals
        use field
        use parameters
        implicit none

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: expchempot 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        ! .. local 
        integer :: i
        real(dp) :: chempot ! chemical potential difference 
        real(dp) :: sumdens 

        if(expchempot==0.0d0) then 
            FEchem_pot=0.0d0
        else     
            sumdens=0.0d0
            if(present(flag)) then  ! water special case because vsol treated different then vi
                chempot = -dlog(expchempot)    
                do i=1,nr
                    sumdens=sumdens +deltaG(i)*xvol(i)
                enddo
                FEchem_pot=delta*chempot*sumdens/vol            
            else
                chempot = -dlog(expchempot/vol)    
                do i=1,nr
                    sumdens=sumdens +deltaG(i)*xvol(i)
                enddo
                FEchem_pot=delta*chempot*sumdens/(vol*vsol)               
            endif
        endif

    end function FEchem_pot   


    real(dp) function FEtrans_entropy_bulk(xvolbulk,vol,flag)
    
        use globals
        use parameters
        implicit none

        real(dp), intent(in) :: xvolbulk 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        if(xvolbulk==0.0d0) then 
            FEtrans_entropy_bulk=0.0d0
        else
            FEtrans_entropy_bulk=0.0d0
            if(present(flag)) then
                ! water special case because vsol treated different then vol  
                FEtrans_entropy_bulk=xvolbulk*(dlog(xvolbulk)-1.0d0)/vol
            else 
                FEtrans_entropy_bulk=xvolbulk*(dlog(xvolbulk/vol)-1.0d0)/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy_bulk   

    real(dp) function FEchem_pot_bulk(xvolbulk,expchempot,vol,flag)
    
        use globals
        use field
        use parameters
        implicit none

        real(dp), intent(in) :: xvolbulk
        real(dp), intent(in) :: expchempot 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

    
        if(expchempot==0.0d0) then 
            FEchem_pot_bulk=0.0d0
        else     
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                FEchem_pot_bulk=-dlog(expchempot)*xvolbulk/vol            
            else
                FEchem_pot_bulk=-dlog(expchempot/vol)*xvolbulk/(vol*vsol)               
            endif
        endif

    end function FEchem_pot_bulk   


    !  chem ical contrubution of polyelectrolyte chain not yet implmentated
    
    subroutine fcnenergy_elect_alternative()
    
        use globals
        use volume
        use parameters
        use field

        implicit none

        !  .. local arguments 
          
        real(dp) :: volumelat          ! volume lattice 
        
        !  .. computation of free energy 
    
        !  .. alternative computation free energy

        call FEconf_elect()

        ! .. translational entropy 

        FEtrans%sol   = FEtrans_entropy(xsol,xbulk%sol,vsol,"w")   
        FEtrans%Na    = FEtrans_entropy(xNa,xbulk%Na,vNa)
        FEtrans%Cl    = FEtrans_entropy(xCl,xbulk%Cl,vCl)
        FEtrans%Ca    = FEtrans_entropy(xCa,xbulk%Ca,vCa)
        FEtrans%K     = FEtrans_entropy(xK,xbulk%K,vK)
        FEtrans%KCl   = FEtrans_entropy(xKCl,xbulk%KCl,vKCl)
        FEtrans%NaCl  = FEtrans_entropy(xNaCl,xbulk%NaCl,vNaCl)
        FEtrans%Hplus = FEtrans_entropy(xHplus,xbulk%Hplus,vsol,"w")
        FEtrans%OHmin = FEtrans_entropy(xOHmin,xbulk%OHmin,vsol,"w")    


        ! .. chemical potential + standard chemical potential 

        FEchempot%sol   = 0.0d0 ! by construction  
        FEchempot%Na    = FEchem_pot(xNa,expmu%Na,vNa)
        FEchempot%Cl    = FEchem_pot(xCl,expmu%Cl,vCl)
        FEchempot%Ca    = FEchem_pot(xCa,expmu%Ca,vCa)
        FEchempot%K     = FEchem_pot(xK,expmu%K,vK) 
        FEchempot%KCl   = FEchem_pot(xKCl,expmu%KCl,vKCl)
        FEchempot%NaCl  = FEchem_pot(xNaCl,expmu%NaCl,vNaCl)
        FEchempot%Hplus = FEchem_pot(xHplus,expmu%Hplus,vsol,"w")
        FEchempot%OHmin = FEchem_pot(xOHmin,expmu%OHmin,vsol,"w")


        ! .. summing all contrubutions
        
        FEalt = FEtrans%sol +FEtrans%Na+ FEtrans%Cl +FEtrans%NaCl+FEtrans%Ca 
        FEalt = FEalt+FEtrans%OHmin +FEtrans%Hplus +FEtrans%K +FEtrans%KCl
        FEalt = FEalt+FEchempot%sol +FEchempot%Na+ FEchempot%Cl +FEchempot%NaCl+FEchempot%Ca 
        FEalt = FEalt+FEchempot%OHmin +FEchempot%Hplus+ FEchempot%K +FEchempot%K+FEchempot%KCl
        
        ! be vary carefull FE = -1/2 \int dr G(r) rho_q(r) psi(r) and FEVdW= - \int dr G(r) chi_IJ ...

        FEalt = FEalt - FEel + FEelSurf + FEconf  !-FEVdW

        ! .. delta translational entropy

        FEtransbulk%sol   = FEtrans_entropy_bulk(xbulk%sol,vsol,"w")   
        FEtransbulk%Na    = FEtrans_entropy_bulk(xbulk%Na,vNa)
        FEtransbulk%Cl    = FEtrans_entropy_bulk(xbulk%Cl,vCl)
        FEtransbulk%Ca    = FEtrans_entropy_bulk(xbulk%Ca,vCa)
        FEtransbulk%K     = FEtrans_entropy_bulk(xbulk%K,vK)
        FEtransbulk%KCl   = FEtrans_entropy_bulk(xbulk%KCl,vKCl)
        FEtransbulk%NaCl  = FEtrans_entropy_bulk(xbulk%NaCl,vNaCl)
        FEtransbulk%Hplus = FEtrans_entropy_bulk(xbulk%Hplus,vsol,"w")
        FEtransbulk%OHmin = FEtrans_entropy_bulk(xbulk%OHmin,vsol,"w")    

        ! .. delta chemical potential + standard chemical potential 

        FEchempotbulk%sol   = 0.0d0 ! by construction  
        FEchempotbulk%Na    = FEchem_pot_bulk(xbulk%Na,expmu%Na,vNa)
        FEchempotbulk%Cl    = FEchem_pot_bulk(xbulk%Cl,expmu%Cl,vCl)
        FEchempotbulk%Ca    = FEchem_pot_bulk(xbulk%Ca,expmu%Ca,vCa)
        FEchempotbulk%K     = FEchem_pot_bulk(xbulk%K,expmu%K,vK) 
        FEchempotbulk%KCl   = FEchem_pot_bulk(xbulk%KCl,expmu%KCl,vKCl)
        FEchempotbulk%NaCl  = FEchem_pot_bulk(xbulk%NaCl,expmu%NaCl,vNaCl)
        FEchempotbulk%Hplus = FEchem_pot_bulk(xbulk%Hplus,expmu%Hplus,vsol,"w")
        FEchempotbulk%OHmin = FEchem_pot_bulk(xbulk%OHmin,expmu%OHmin,vsol,"w")

        
        ! .. bulk free energy
 
        FEbulkalt = FEtransbulk%sol +FEtransbulk%Na+ FEtransbulk%Cl +FEtransbulk%NaCl+FEtransbulk%Ca 
        FEbulkalt = FEbulkalt+FEtransbulk%OHmin +FEtransbulk%Hplus +FEtransbulk%K +FEtransbulk%KCl
        FEbulkalt = FEbulkalt+FEchempotbulk%sol +FEchempotbulk%Na+FEchempotbulk%Cl +FEchempotbulk%NaCl+FEchempotbulk%Ca 
        FEbulkalt = FEbulkalt+FEchempotbulk%OHmin +FEchempotbulk%Hplus -FEchempotbulk%K -FEchempotbulk%KCl
        
        !   .. volume lattice divide by area NP
        volumelat = (4.0d0/3.0d0)*pi*((delta*nr+radius)**3-radius**3)   
        volumelat = volumelat/(4.0d0*pi*(radius**2))
        FEbulkalt = volumelat*FEbulkalt

        ! .. delta

        deltaFEtrans%sol   = FEtrans%sol  - FEtransbulk%sol * volumelat
        deltaFEtrans%Na    = FEtrans%Na   - FEtransbulk%Na * volumelat
        deltaFEtrans%Cl    = FEtrans%Cl   - FEtransbulk%Cl * volumelat
        deltaFEtrans%Ca    = FEtrans%Ca   - FEtransbulk%Ca * volumelat
        deltaFEtrans%K     = FEtrans%K    - FEtransbulk%K * volumelat
        deltaFEtrans%KCl   = FEtrans%KCl  - FEtransbulk%KCl * volumelat
        deltaFEtrans%NaCl  = FEtrans%NaCl - FEtransbulk%NaCl * volumelat
        deltaFEtrans%Hplus = FEtrans%Hplus- FEtransbulk%Hplus * volumelat
        deltaFEtrans%OHmin = FEtrans%OHmin- FEtransbulk%OHmin * volumelat
         
        deltaFEchempot%sol   = FEchempot%sol  - FEchempotbulk%sol * volumelat
        deltaFEchempot%Na    = FEchempot%Na   - FEchempotbulk%Na * volumelat
        deltaFEchempot%Cl    = FEchempot%Cl   - FEchempotbulk%Cl * volumelat
        deltaFEchempot%Ca    = FEchempot%Ca   - FEchempotbulk%Ca * volumelat
        deltaFEchempot%K     = FEchempot%K    - FEchempotbulk%K * volumelat
        deltaFEchempot%KCl   = FEchempot%KCl  - FEchempotbulk%KCl * volumelat
        deltaFEchempot%NaCl  = FEchempot%NaCl - FEchempotbulk%NaCl * volumelat
        deltaFEchempot%Hplus = FEchempot%Hplus- FEchempotbulk%Hplus * volumelat
        deltaFEchempot%OHmin = FEchempot%OHmin- FEchempotbulk%OHmin * volumelat


        ! .. differences

        deltaFEalt = FEalt - FEbulkalt

    !    print*,"FEbulkalt=",FEbulkalt

    !        print*,"delta FEchemsurfalt(RIGHT)=",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)



    end subroutine fcnenergy_elect_alternative
    

end module energy
