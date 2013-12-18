subroutine second_membre(B,U,Utop,Ubottom,Me,Np,i1,iN,Nx,dx,dy,Cx,Cy)
    USE utils
    implicit none
    integer, intent(in) :: i1,iN,Nx,Me,Np
    double precision, intent(in)  :: dx,dy,Cx,Cy
    !double precision :: left,right,bottom,top
        double precision  :: posx,posy
    double precision, dimension(i1:iN), intent(out):: B
    double precision, dimension(i1:iN), intent(in):: U
    double precision, dimension(1:Nx) :: Utop,Ubottom

    integer :: i,j,k,l,M

    B=0.d0
    M = (iN-i1+1)/Nx !Taille du sous domaine en # de lignes
    DO i = i1,iN
        CALL nloc(j,k,i,Nx)
        posx = (k)*dx
        posy = (j)*dy
        B(i) = f(posx,posy,0.d0)
    ENDDO
!   IF( me == 0 )THEN
        l = 1
        DO i = i1, i1+(Nx-1)
            CALL nloc(j,k,i,Nx)
            posx = (k)*dx
            posy = (j)*dy
            Ubottom(l) = bottom(k*dx,0.0d0,0.d0)
            l = l+1
        ENDDO
!   ENDIF
    l = 1
!   IF( me == Np-1 )THEN
        DO i = iN-Nx+1,iN
            CALL nloc(j,k,i,Nx)
            posx = (k)*dx
            posy = (j)*dy
            Utop(l) = top(posx,1.d0,0.0d0)
            l = l+1
        ENDDO
!   ENDIF
!Ligne du bas
l =1
do j= i1,i1+(Nx-1)
    B(j) = B(j) - Ubottom(l)*Cy
    l = l+1
enddo
CALL nloc(i,j,i1,Nx)
B(i1) = B(i1) -left(0.d0,i*dy,0.d0)*Cx
B(i1+Nx-1) = B(i1+Nx-1) -right(1.d0,i*dy,0.d0)*Cx

!Ligne du milieu
j = i1+Nx
do i = 2, M-1
    call nloc(k,l,j,Nx)
    B(j) = B(j) -left(0.d0,k*dy,0.d0)*Cx
    B(j+Nx-1) = B(j+Nx-1) -right(1.d0,k*dy,0.d0)*Cx
    j = i1 + (i)*Nx

enddo
!ligne du haut
call nloc(k,l,iN,Nx)
B(iN-Nx+1) = B(iN-Nx+1) -left(0.d0,k*dy,0.d0)*Cx
B(iN) = B(iN) -right(1.d0,k*dy,0.d0)*Cx
l  = 1
do j = iN-Nx+1,iN
    B(j) = B(j) - Utop(l)*Cy
        l = l+1
enddo

end subroutine second_membre
