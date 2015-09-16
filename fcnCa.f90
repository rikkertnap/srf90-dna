!---------------------------------------------------------------|
! fcnCa.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes on a coated onto a spherical surface             |
!---------------------------------------------------------------|


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

    ! pre:  vector x=xsol+psi+rhopolA+rhoolB+xpolC                   
    ! post: vector f=xpolAB+xpolC+xsol+Sum_i x_i -1,poisson equation,rhoA,rhoB,xpolC              

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

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize),xpolCin(nsize)
        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,k,c,s         ! dummy indices
        real(dp) :: norm


        real(dp), parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nr                      ! size vector neq=5*nr x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                  ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
            xpolCin(i)=x(i+4*n)    ! volume fraction C-polymer
       enddo

        psiSurf = psi(1)          ! surface potentail

        do i=1,n                  ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            xpolC(i)   = 0.0d0     ! C polymer volume fraction 
            rhopolA(i) = 0.0d0     ! A polymer density 
            rhopolB(i) = 0.0d0
            rhopolC(i) = 0.0d0
            xpi(i) = xsol(i)/dexp(-chiCS*xpolCin(i)-(chiAS*rhopolAin(i)*vpolA(3)+&
                    chiBS*rhopolBin(i)*vpolB(3))*vsol)  

            xNa(i)   = expmu%Na*(xpi(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xpi(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xpi(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xpi(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xpi(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xpi(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xpi(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xpi(i))*dexp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xpi(i)**deltavA(1)))     ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xpi(i)**deltavA(2)))  ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xpi(i)**deltavA(3)))  ! ACa+/A-
       
            sumxA=xA(1)+xA(2)+xA(3)
            constA=(2.0d0*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xpi(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0d0/(1.0d0+sumxA)
                fdisA(5,i)=0.0d0
            else
                fdisA(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constA/((sumxA+1.0d0)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0d0)/(2.0d0*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xpi(i) **deltavB(1)))    ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xpi(i)**deltavB(2)))  ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xpi(i)**deltavB(3)))  ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0d0*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xpi(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0d0/(1.0d0+sumxB)
                fdisB(5,i)=0.0d0
            else
                fdisB(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constB/((sumxB+1.0d0)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0d0)/(2.0d0*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       

    !        exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
    !        exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

    !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
    !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
     
            ! Na condensed ANa reference state  
            exppiA(i) = (xpi(i)**vpolA(3))*dexp(-chiAS*xsol(i)*vpolA(3)-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i) = (xpi(i)**vpolB(3))*dexp(-chiBS*xsol(i)*vpolB(3)-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
    
            exppiC(i)=(xpi(i)**vpolC)*dexp(-chiCS*xsol(i)*vpolC)! auxiliary variable
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0d0                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0d0                 ! init q   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0d0               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            do k=1,4               ! polymer volume fraction
                xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                    +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo    
            xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0d0)
            xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0d0)
       
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0d0
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
            !   ..  total charge density in units of vsol
        enddo  !  .. end computation polymer density and charge density  

        ! .. electrostatics 

        sigmaqSurf=0.0d0 ! charge regulating surface charge 
        psi(n+1)=0.0d0   ! bulk potential

        !    .. Poisson Eq 

        f(n+1)= -0.5d0*(Fplus(1)*(psi(2)-psi(1)) + Fmin(1)*sigmaqSurf +rhoq(1)*constqW)      !     boundary

        do i=2,n
            f(n+i)= -0.5d0*(Fplus(i)*psi(i+1)-2.0d0*psi(i) + Fmin(i)*psi(i-1) +rhoq(i)*constqW)
        enddo

        do i=1,n
            f(2*n+i)=rhopolA(i)-rhopolAin(i)
            f(3*n+i)=rhopolB(i)-rhopolBin(i)
            f(4*n+i)=xpolC(i)-xpolCin(i)
        enddo

        norm=l2norm(f,5*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnelect



    subroutine fcnneutralAB(x,f,nn)

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
                real(dp) :: L2norm
                integer, intent(in) :: n
                real(dp), intent(in) :: f(n)
            end function L2norm 
        end interface

        !     .. declare local variables


        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: xpolBin(nsize),xpolAin(nsize)
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,k,c,s         ! dummy indices
        real(dp) :: norm
        real(dp), parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 

        !     .. executable statements 

        n=nr                      ! size vector neq=3*nr x(pi,rhopolA,rhopolB)

        do i=1,n                  ! init x 
            xsol(i)= x(i)         ! solvent volume fraction 
            xpolBin(i)=x(i+n)     ! volume fraction B-polymer 
            xpolAin(i)=x(i+2*n)     ! volume fraction A-polymer 
        enddo

        do i=1,n                   ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            xpolC(i)   = 0.0d0     ! C polymer volume fraction 
            rhopolA(i) = 0.0d0     ! A polymer density 
            rhopolB(i) = 0.0d0
            rhopolC(i) = 0.0d0
        
            xpi(i) = xsol(i)/dexp(-chiAS*xpolAin(i)-chiBS*xpolBin(i))
            exppiC(i) = (xpi(i)**vpolC)
            exppiB(i) = (xpi(i)**vpolB(3))*dexp(-chiBS*xsol(i)*vpolB(3)) 
            exppiA(i) = (xpi(i)**vpolA(3))*dexp(-chiAS*xsol(i)*vpolA(3))
             
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0d0                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0d0                 ! init q   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0d0               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            ! polymer volume fraction A and B all in state 3 ANa and BNa
            xpolAB(i)=(rhopolA(i)*vpolA(3)+rhopolB(i)*vpolB(3))*vsol
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)-1.0d0
            f(n+i)=rhopolB(i)*vpolB(3)*vsol-xpolBin(i)
            f(2*n+i)=rhopolA(i)*vpolA(3)*vsol-xpolAin(i)

        enddo  !  .. end computation polymer density 

        norm=l2norm(f,3*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnneutralAB



    subroutine fcnneutralB(x,f,nn)

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
        integer*8, intent(in) :: nn

        interface 
            function L2norm(f,n)
                use precision_definition
                real(dp) :: L2norm
                integer, intent(in) :: n
                real(dp), intent(in) :: f(n)
            end function L2norm 
        end interface

        !     .. declare local variables

        real(dp) :: exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: xpolBin(nsize)
        real(dp) :: xpolCin(nsize)
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,k,c,s         ! dummy indices
        real(dp) :: norm
       
        ! real(dp) :: cn                ! auxilary variable for Poisson Eq

        real(dp), parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 

        !     .. executable statements 

        n=nr                      ! size vector neq=3*nr x(pi,rhopolA,rhopolB)

        do i=1,n                  ! init x 
            xsol(i)= x(i)         ! solvent volume fraction 
            xpolBin(i)=x(i+n)     ! volume fraction B-polymer 
            xpolCin(i)=x(i+2*n)     ! volume fraction B-polymer 
        enddo

        do i=1,n                   ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            xpolC(i)   = 0.0d0     ! C polymer volume fraction 
            rhopolB(i) = 0.0d0
            rhopolC(i) = 0.0d0
        
            xpi(i) = xsol(i)/dexp(-chiBS*xpolBin(i)-chiCS*xpolCin(i))
            exppiC(i)=(xpi(i)**vpolC)*dexp(-chiCS*xsol(i)*vpolC)
            exppiB(i) = (xpi(i)**vpolB(3))*dexp(-chiBS*xsol(i)*vpolB(3)) !*dexp(-zpol*psi(i))/fdisB(i) 
            !     .. VdW interaction   
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0d0                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                pro = pro*exppiB(k)
             enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                rhopolB(k)=rhopolB(k)+pro
            enddo
        enddo

        qC = 0.0d0                 ! init q   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0d0               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= 0.0d0
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            ! polymer volume fraction B only  all in state 3  BNa
            xpolAB(i)=rhopolB(i)*vpolB(3)*vsol
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)-1.0d0
            f(n+i)=rhopolB(i)*vpolB(3)*vsol-xpolBin(i)
            f(2*n+i)=xpolC(i)-xpolCin(i)
         enddo  !  .. end computation polymer density 

        norm=l2norm(f,3*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnneutralB








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
