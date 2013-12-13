    subroutine jacobi(Aii,Cx,Cy,Nx,i1,iN,B,U,Uold)
    implicit none

    double precision  :: Aii,Cx,Cy, invAii,err,Tol
    integer     :: max_iter,i1,iN,Nx,M
    double precision, dimension(i1:iN), intent(inout) :: U,Uold
    double precision, dimension(i1:iN), intent(in) :: B

    Tol = 1.d-12!1d-18
    invAii = 1.d0/Aii
!~  j = 0
    max_iter = 100000
    M = (iN-i1+1)/Nx
    err = 100.d0

!~ Do while( (err > Tol) .AND. (j < max_iter) )
!~         Uold(i1:iN)=U(i1:iN)
        call MATVEC(0.d+00,Cx,Cy,Nx,M,i1,Uold,U)
        U(i1:iN)=(B(i1:iN)-U(i1:iN))*invAii

        !pour le calcul de l'erreur
!~         Uold(i1:iN)=U(i1:iN)-Uold(i1:iN)
!~         err=sqrt(dot_product(Uold(i1:iN),Uold(i1:iN)))

!~         j=j+1
!~ enddo
!~  write(*,*) 'fin jacobi, convergence en ', j, 'iterations', err

    end subroutine jacobi
