    subroutine jacobi(Aii,Cx,Cy,Nx,i1,iN,B,U,Uold)
    implicit none

    double precision  :: Aii,Cx,Cy, invAii
    integer     :: i1,iN,Nx,M
    double precision, dimension(i1:iN), intent(inout) :: U,Uold
    double precision, dimension(i1:iN), intent(in) :: B

    invAii = 1.d0/Aii
    M = (iN-i1+1)/Nx

	Uold(i1:iN)=U(i1:iN)
	call MATVEC(0.d+00,Cx,Cy,Nx,M,i1,Uold,U)
	U(i1:iN)=(B(i1:iN)-U(i1:iN))*invAii

    end subroutine jacobi
