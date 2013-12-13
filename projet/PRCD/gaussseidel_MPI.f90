   SUBROUTINE gauss(Aii,Cx,Cy,Nx,i1,iN,B,U, Uold)
    implicit none
    integer :: j, i, k, min_boucle, max_boucle
    integer :: i1,iN, max_iter, Nx
    real*8  :: Aii,Cx,Cy,invAii,err,Tol, M
    real*8, DIMENSION(i1:iN)    :: U,B, Uold

    j = 0
    max_iter = 100000
    err = 100.d0
    invAii = 1.d0/Aii
    Tol = 1.d-12
    M = (iN-i1+1)/Nx

!~ http://fr.wikipedia.org/wiki/M%C3%A9thode_de_Gauss-Seidel

!~ boucle methode gauss-eidel
Do while( (err > Tol) .AND. (j < max_iter) )
    Uold(i1:iN)=U(i1:iN)
    DO i=i1, iN
        U(i)=B(i)
        min_boucle = MIN(MAX(i1, i-1),MAX(i1, i-2))
        max_boucle = MAX(MIN(iN, i+1) ,MIN(iN, i+2))
        DO k=min_boucle,max_boucle
            if( ((k-1)== i) .OR. ((k+1) ==i)) then
                U(i) = U(i) - Cx*U(k)
            elseif((k-2)==i .OR. (k+2 ==i))then
                U(i) = U(i) - Cy*U(k)
            endif
        ENDDO
        U(i)=U(i)*invAii
    ENDDO
    j=j+1
    !pour le calcul de l'erreur
    Uold(i1:iN)=U(i1:iN)-Uold(i1:iN)
    err=sqrt(dot_product(Uold(i1:iN),Uold(i1:iN)))

enddo
    write(*,*) 'fin gauss, convergence en ', j, 'iterations', err

   RETURN
   END SUBROUTINE gauss
