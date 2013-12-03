program schwarz_additif
!
! Resolution du placien par l'algorithme de schwarz-additif
!   ∆U = f
! 
        USE UTILS
	implicit none
!	include 'mpif.h'

	! Declaration des variables pour schwarz
	real*8, parameter  :: eps = 1.d-12
	real*8	:: ERR,ERR_glob
	! Decalaration des variables MPI et charges	
        integer	:: itop,ibottom	 !indice de localisation des debuts de msg a envoyer
	! Declaration des variables numeriques
	integer	:: l,Max_l
	integer :: R
	real*8	:: t_debut,t_fin,t_max
	! Declaration des variables du pb a resoudre
	integer	:: Nx,Ny,N
	real*8	:: dx,dy,Lx,Ly,posx,posy
	real*8	:: Aii,Cx,Cy,D
	real*8, dimension(:), allocatable :: U, Uold, RHS, Uschwarz
	real*8, dimension(:), allocatable :: Utop, Ubottom
	integer	:: i,j,k

        Ny = 10!*Np
        Nx = 10 
!       R  = 1 !Recouvrement en Nbr de lignes (1 est le minimum)
        Lx = 1.0d0
        Ly = 1.0d0 
	D  = 1.0d0 !coefficient de diffusion de l equation
	!Lecture des variables dans un fichier
	open(11, file='param', status ='old')
	read(11,*) Ny !nb lignes
	read(11,*) Nx !nb colonnes
	read(11,*) R !recouvrementS
	read(11,*) Lx
	read(11,*) Ly
	read(11,*) D
	close(11)
	
	
	
        dx = Lx/(1+Nx)
        dy = Ly/(1+Ny)
	Aii = 2.d0*D/(dx*dx)+2.d0*D/(dy*dy) ! Terme diagonal de la matrice
	Cx  = -1.d0*D/(dx*dx)	! Terme extra-diagonale proche
	Cy  = -1.d0*D/(dy*dy)	! Terme extra-diagonale eloigne
        N  = Nx*Ny

        ibottom = 1
        itop = N
        
        ! Allocation dynamique pour chaque proc en tenant compte du recouvrement
        ALLOCATE(U(1:N));ALLOCATE(RHS(1:N));ALLOCATE(Uold(1:N));
	! Allocation dynamique pour la communication
	ALLOCATE(Utop(1:Nx));ALLOCATE(Ubottom(1:Nx))

	! Initialisation du probleme
        U(1:N)    = 0.0d0
	Utop	    = 0.0d0
	Ubottom	    = 0.0d0

!	schw_iter = 0
	ERR_glob = 1.d3
		! Remplissage du second membre incluant les conditions de Bords Dirichlet
		call second_membre(RHS,U,Utop,Ubottom,77,77,1,N,Nx,dx,dy,Cx,Cy)
		! Resolution du systeme lineaire
		call jacobi(Aii,Cx,Cy,Nx,1,N,RHS,U,Uold,l)
!~ 		call   grad(Aii,Cx,Cy,Nx,1,N,RHS,U,l)

	call wrisol( U,Nx,Ny,dx,dy,77,1,N )
	DEALLOCATE(U,RHS,Uold)
	DEALLOCATE(Utop,Ubottom)
   end program schwarz_additif
