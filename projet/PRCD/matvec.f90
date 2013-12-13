subroutine matvec(Aii,Cx,Cy,Nx,m,i1,U,AU)

! Produit matrice vecteur dans le cas ou A pentadiagonale de
! la forme:
!
! A = B C             matrice pentadiagonale (m*Nx,m*Nx)
!     C B C
!       C B C
!         . . .
!          . . .
!            C B
! avec
! C = Cy Id            matrice diagonale (Nx,Nx)
!
! B = Aii Cx           matrice tridiagonale (Nx,Nx)
!     Cx  Aii Cx
!          .   .   .
!              Cx Aii

    implicit none
    integer, intent(in) :: Nx,m,i1
    double precision, intent(in)  :: Aii,Cx,Cy
    double precision, dimension(i1:Nx*m), intent(in)  :: U
    double precision, dimension(i1:Nx*m), intent(out) :: AU
    integer :: i,j,k

!Premier bloc
AU(i1) = Aii*U(i1) + Cx*U(i1+1) + Cy*U(i1+Nx)
i = i1
DO j = 1, Nx-2
    i = i1+j
    AU(i) = Aii*U(i) + Cx*U(i-1) + Cx*U(i+1) + Cy*U(i+Nx)
ENDDO
AU(i1+(Nx-1)) = Aii*U(i1+(Nx-1)) + Cx*U(i1+(Nx-1)-1) + Cy*U(i1+(Nx-1)+Nx)

i = i1 + (Nx-1)

!bloc general, il y a m-2 blocs generaux
DO k = 1, (m-2)
    !Premiere ligne
    i = i+1
    AU(i) = Aii*U(i) + Cx*U(i+1) + Cy*U(i-Nx) + Cy*U(i+Nx)
    !ligne generale
    DO j = 1, Nx-2
        i = i+1
        AU(i) = Aii*U(i) + Cx*U(i-1) + Cx*U(i+1) + Cy*U(i+Nx) + Cy*U(i-Nx)
    ENDDO
    !Derniere ligne
    i = i + 1
    AU(i) = Aii*U(i) + Cx*U(i-1) + Cy*U(i-Nx) + Cy*U(i+Nx)
ENDDO
i = i+1
!Dernier bloc
AU(i) = Aii*U(i) + Cx*U(i+1) + Cy*U(i-Nx)
DO j = 1, Nx-2
    i = i+1
    AU(i) = Aii*U(i) + Cx*U(i+1) + Cx*U(i-1) + Cy*U(i-Nx)
ENDDO
i = i+1
AU(i) = Aii*U(i) + Cx*U(i-1) + Cy*U(i-Nx)

end subroutine matvec
