MODULE utils
   implicit none

CONTAINS
   SUBROUTINE charge(Me,Np,N,i1,iN)
      implicit none
      integer :: i1,iN,Me,Np,N
      integer :: r,q
       
       q = N/Np
       r = N - q*Np
       
       IF( Me < r ) THEN
          i1 = Me*( q+1 ) + 1
	  iN = ( Me + 1 )*( q + 1 )
       ELSE
          i1 = 1 + r + Me*q
	  iN = i1 + q - 1
       END IF

   END SUBROUTINE charge

   SUBROUTINE nloc(i,j,n,Nx)
      integer :: n,Nx
      integer :: i,j
      integer :: q,r

      q = n/Nx
      r = n - q*Nx
      if ( r == 0 )then
           i = q
           j = Nx
      else
           i = 1+q
           j = r
      endif


   END SUBROUTINE nloc

   SUBROUTINE Rename(Me,name)
      implicit none
      integer :: Me
      character*13  ::name
      character*3 :: tn
      integer :: i1,i2,i3

      i1 = Me/100
      i2 =( Me - 100*i1)/10
      i3 = Me - 100*i1 -10*i2
      tn = char(i1+48)//char(i2+48)//char(i3+48)
      name='sol'//tn//'.dat'

   END SUBROUTINE Rename



   SUBROUTINE wrisol(U,Nx,Ny,dx,dy,Me,i1,iN)
      implicit none
      integer :: Me,Nx,Ny,i1,iN
      real*8  :: dx,dy
      real*8, dimension(i1:iN) :: U
      integer :: unit,i,j,k
      real*8  :: posx,posy
      CHARACTER(LEN=13) :: name

      CALL RENAME(Me,name) 

      unit = 7  
      open(unit,file=name,status='unknown',form='formatted')
      DO i = i1, iN
           CALL nloc(j,k,i,Nx)
           posx = (k)*dx
           posy = (j)*dy
           write(unit,*) posx,posy,U(i)
       ENDDO
       close(unit)

      return
   END SUBROUTINE wrisol

   FUNCTION right(x,y,t)
	implicit none
	real*8, intent(in)	:: x,y,t
	real*8	:: right

!~ 	right = 0.d0
	right = sin(x)+cos(y)
   END FUNCTION right

   FUNCTION left(x,y,t)
	implicit none
	real*8, intent(in)	:: x,y,t
	real*8	:: left

!~ 	left = 0.d0
	left = sin(x)+cos(y)
   END FUNCTION left

   FUNCTION bottom(x,y,t)
	implicit none
	real*8, intent(in)	:: x,y,t
	real*8	:: bottom

!~ 	bottom = 0.d0
	bottom = sin(x)+cos(y)
   END FUNCTION bottom

   FUNCTION top(x,y,t)
	implicit none
	real*8, intent(in)	:: x,y,t
	real*8	:: top

!~ 	top = 0.0d0
	top = sin(x)+cos(y)
   END FUNCTION top

   FUNCTION f(x,y,t)
	implicit none
	real*8	:: x,y,t
	real*8	:: f

!~ 	f = 2.d0*(y - y*y + x - x*x) 
	f = sin(x)+cos(y)
   END FUNCTION f

end MODULE utils
