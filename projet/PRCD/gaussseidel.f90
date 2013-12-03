   SUBROUTINE gauss(Aii,Cx,Cy,Nx,i1,iN,B,U,l)
      implicit none
	integer, intent(out)	:: l
	integer			:: i1,iN,i,Nx,Nl,M
	real*8	:: residu,alpha,beta,drl,dwl, epsilon
	real*8	:: Aii,Cx,Cy
	real*8, DIMENSION(i1:iN)	:: U,B,r,kappa,d,W
      
!~     formule :(L + D)x^(k+1) = b-Ux^k par blocs
!~     D partie diagonale, L partie triangulaire basse, U partie triangulaire haute
!~ 	   U = A-L-D
	j = 0
	max_iter = 10000
	M = (iN-i1+1)/Nx
	err = 100.d0


      M = (iN-i1+1)/Nx !nombre de blocs

!~ http://fr.wikipedia.org/wiki/M%C3%A9thode_de_Gauss-Seidel

!~ initialisation du vecteur x (vecteur par bloc)
!~ x0, x1 (besoin d'une copie) 


!~ boucle methode gausseidel
Do while( (err > Tol) .AND. (j < max_iter) )
	DO l=1, M
	!~ formule :(L + D)x^(k+1) = b-Ux^k par blocs
    !~ A11 X1^(k+1)1 = b - U X^k solve(AX=B)
	!~ ...	

	ENDDO
enddo

   RETURN
   END SUBROUTINE res
