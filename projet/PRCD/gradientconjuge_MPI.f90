   SUBROUTINE grad(Aii,Cx,Cy,Nx,i1,iN,B,U,l)
      implicit none
    integer, intent(out)    :: l
    integer         :: i1,iN,i,Nx,Nl,M
    double precision  :: residu,alpha,beta,drl,dwl, epsilon
    double precision  :: Aii,Cx,Cy
    double precision, DIMENSION(i1:iN)    :: U,B,r,kappa,d,W


      epsilon = 1.d-12
      Nl =100000
      M = (iN-i1+1)/Nx

      !initialisation Gradient conjugue
      kappa(i1:iN) = U(i1:iN)
      call MATVEC(Aii,Cx,Cy,Nx,M,i1,kappa,r)
      r(i1:iN)     = r(i1:iN) - B(i1:iN)
      residu = SUM(r(i1:iN)*r(i1:iN))

      d(i1:iN)=r(i1:iN)

      ! boucle du Gradient conjugue
      DO l=1, Nl
         IF( SQRT(residu) .lt. epsilon ) EXIT

         call MATVEC(Aii,Cx,Cy,Nx,M,i1,d,W)
         drl = 0.0d0
         dwl = 0.0d0
         DO i = i1, iN
            drl = drl + d(i)*r(i)
            dwl = dwl + d(i)*w(i)
         ENDDO

         alpha = drl/dwl
         kappa(i1:iN) = kappa(i1:iN) - alpha*d(i1:iN)
         r(i1:iN) = r(i1:iN) - alpha*W(i1:iN)
         beta = SUM(r(i1:iN)*r(i1:iN))/residu
         d(i1:iN) = r(i1:iN) + beta*d(i1:iN)
         residu    = SUM(r(i1:iN)*r(i1:iN))
         !Fin Gradient conjugue
       ENDDO

       U(i1:iN) = kappa(i1:iN)
!~        print*, '        Gradient Conjugue',l,residu

   RETURN
   END SUBROUTINE grad
