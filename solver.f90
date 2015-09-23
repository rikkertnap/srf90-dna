subroutine solver(x, xguess, accuracy, residual)
  
    use globals
    use parameters
    use listfcn
    use fcnpointer

    
    implicit none
  
    !     .. variables and constant declaractions

    !     .. array arguments
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)  
    !     .. scalar arguments
    real(dp) :: accuracy
    real(dp) :: residual
 
    call set_size_neq()
 
    if(sysflag=="elect") then
        fcnptr => fcnelect
    elseif(sysflag=="bulk water") then 
        fcnptr => fcnbulk
!     elseif(sysflag=="neutralAB") then 
!         fcnptr => fcnneutralAB
!     elseif(sysflag=="neutral") then
!         fcnptr => fcnneutralB
    else
        print*,"Wrong value sysflag : ", sysflag
        stop
    endif

 
    if(method.eq."kinsol") then
     
        call kinsol_gmres_solver(x, xguess, neq, accuracy, residual)
    
     !  elseif (method.eq."zspow") then 
     !     call zspow_solver(x, xguess, neq, accuracy, residual)
     
    else  
        print*,"Solver method incorrect"
        stop
    endif
  
end subroutine solver
