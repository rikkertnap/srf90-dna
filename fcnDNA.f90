! --------------------------------------------------------------|
! fcnCa.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes on a coated onto a spherical surface             |
! --------------------------------------------------------------|


module fcnpointer

    implicit none

    abstract interface
        subroutine fcn(x,f,n)
            use precision_definition
            implicit none
            real(dp), dimension(n), intent(in) :: x
            real(dp), dimension(n), intent(out) :: f
            integer*8, intent(in) :: n
        end subroutine fcn
    end interface

    procedure(fcn), pointer :: fcnptr => null()

end module fcnpointer


module listfcn

    implicit none

    contains

    ! pre:  vector x=xsol+psi                   
    ! post: vector f=xpol+Sum_i x_i -1,poisson equation              

    subroutine fcnelect(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        interface 
            function L2norm(f,n)
                use precision_definition
                implicit none
                real(dp) :: L2norm
                integer, intent(in) :: n
                real(dp), intent(in) :: f(n)
            end function L2norm 
        end interface

        !     .. declare local variables

        real(dp) :: exppi(nsize,nsegtypes)   ! auxilairy variable for computing P(\alpha) 
        real(dp) :: pro,rhopol0
        integer :: n                 ! half of n   
        integer :: i,k,c,s,t         ! dummy indices
        real(dp) :: norm


        real(dp), parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nr                      ! size vector neq=2*nr x=(pi,psi)

        do i=1,n                  ! init x 
            xsol(i)= x(i)         ! solvent volume fraction 
            psi(i) = x(i+n)       ! potential
        enddo

        psiSurf = psi(1)          ! bc surface potentail

        do i=1,n                  ! init volume fractions 
    
            xNa(i)    = expmu%Na   *(xsol(i)**vNa  )*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xCl(i)    = expmu%Cl   *(xsol(i)**vCl  )*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i)       )*dexp(-psi(i))     ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i)       )*dexp(+psi(i))     ! OH- volume fraction
            xK(i)     = expmu%K    *(xsol(i)**vK   )*dexp(-psi(i)*zK ) ! ion plus volume fraction
            xCa(i)    = expmu%Ca   *(xsol(i)**vCa  )*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i)  = expmu%NaCl *(xsol(i)**vNaCl)                   ! ion pair volume fraction
            xKCl(i)   = expmu%KCl  *(xsol(i)**vKCl )                   ! ion pair volume fraction

        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,n
                    rhopol(i,t) = 0.0d0                                                ! init polymer density          
                    fdis(i,t)  = 1.0d0/(1.0d0+xHplus(i)/(K0a(t)*xsol(i)))      
                    exppi(i,t)  = (xsol(i)**vpol(t))*exp(zpol(t,1)*psi(i) )/fdis(i,t)   ! auxilary variable palpha
                enddo  
            else
                do i=1,n
                    rhopol(i,t) = 0.0d0
                    fdis(i,t)  = 1.0d0
                    exppi(i,t)  = xsol(i)**vpol(t)
                enddo  
            endif   
        enddo   
     
        !   .. computation polymer volume fraction 

        q = 0.0d0                 ! init q

        do c=1,cuantas               ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nseg              ! loop over segments 
                k=indexchain(c,s)
                t=type_of_monomer(s)
                pro = pro*exppi(k,t)
            enddo

            q= q+pro
            do s=1,nseg
                k=indexchain(c,s)
                t=type_of_monomer(s)
                rhopol(k,t)=rhopol(k,t)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopol0=sigma/q

        do i=1,n
            xpol(i) =0.0d0
            rhoq(i) = 0.0d0
            do t=1,nsegtypes
                rhopol(i,t)= rhopol0*rhopol(i,t)/deltaG(i) ! density of polymer type t
                xpol(i)=xpol(i)+rhopol(i,t)*vpol(t)*vsol 
                rhoq(i)=(zpol(t,1)*fdis(i,t)+zpol(t,2)*(1.0d0-fdis(i,t)))*rhopol(i,t)*vsol ! charge from polymer 
            enddo

            f(i)=xpol(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0d0
       
            rhoq(i)= rhoq(i)+zNa*xNa(i)/vNa +zCa*xCa(i)/vCa +zK*xK(i)/vK &
                +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)
                      
            !   ..  total charge density in units of vsol
        enddo  
        !  .. end computation polymer density and charge density  

        ! .. electrostatics 

        sigmaqSurf=0.0d0 ! surface charge 
        psi(n+1)=0.0d0   ! bulk potential

        !    .. Poisson Eq 

        f(n+1)= -0.5d0*(Fplus(1)*(psi(2)-psi(1)) + Fmin(1)*sigmaqSurf +rhoq(1)*constqW)      !     boundary

        do i=2,n
            f(n+i)= -0.5d0*(Fplus(i)*psi(i+1)-2.0d0*psi(i) + Fmin(i)*psi(i-1) +rhoq(i)*constqW)
        enddo


        norm=l2norm(f,2*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnelect




    !     .. function solves for bulk volume fraction 

    subroutine fcnbulk(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use physconst

        implicit none

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out):: f(neq)

        interface 
            function L2norm(f,n)
                use precision_definition
                real(dp) :: L2norm
                integer, intent(in) :: n
                real(dp), intent(in) :: f(n)
            end function L2norm
        end interface

        !     .. local variables

        real(dp) :: phiNaCl,phiNa,phiCl,phiK,phiKCl
        real(dp) :: deltavolNaCl,deltavolKCl,norm

        !     .. executable statements 


        deltavolNaCl=(vNaCl-vNa-vCl)
        deltavolKCl=(vKCl-vNa-vCl)

        phiNa   =x(1)
        phiCl   =x(2)
        phiNaCl =x(3)
        phiK    =x(4)
        phiKCl  =x(5)

        !     .. condensation equilibrium eq NaCl

        f(1)= phiNaCl-(phiNa*phiCl*K0ionNa*vsol*vNaCl/(vNa*vCl*vsol))*      & 
         ((1.0d0 -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolNaCl)

        !     .. charge neutrality
        f(2)=phiNa/vNa-phiCl/vCl+2.0d0*xbulk%Ca/vCa+xbulk%Hplus-xbulk%OHmin+phiK/vK

        !     .. conservation of number NaCl

        f(3)=phiNa/vNa+phiCl/vCl+2.0d0*phiNaCl/vNaCl + phiKCL/vKCL &    
             -2.0d0*vsol*(Na/1.0d24)*cNaCl-dabs(xbulk%Hplus-xbulk%OHmin) &
            -2.0d0*xbulk%Ca/vCa-vsol*(Na/1.0d24)*cKCl

        !     .. condensation equilibrium eq KCl

        f(4)= phiKCl-(phiK*phiCl*K0ionK*vsol*vKCl/(vK*vCl*vsol))* &
             ((1.0d0 -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolKCl)

        !     .. conservation of number K
        f(5)= phiK/vK+phiKCl/vKCl-vsol*(Na/1.0d24)*cKCl

        norm=l2norm(f,5)
        iter=iter+1
     
        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulk


    end module listfcn
