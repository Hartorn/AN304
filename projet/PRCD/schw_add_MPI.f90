program schwarz_additif
!
! Resolution du placien par l'algorithme de schwarz-additif
!   ∆U = f
! 
        USE UTILS
	implicit none
	include 'mpif.h'

	! Declaration des variables pour schwarz
	real*8, parameter  :: eps = 1.d-12
	real*8	:: error_calcul ,ERR_glob
	! Decalaration des variables MPI et charges	
        integer	:: itop,ibottom	 !indice de localisation des debuts de msg a envoyer
	! Declaration des variables numeriques
	integer	:: l,Max_l
	integer :: R
	real*8	:: t_debut,t_fin,t_max
	! Declaration des variables du pb a resoudre
	integer	:: Nx,Ny,N, myrank, wsize, IERROR, Nb_ligne_proc, max_iter, reste
	real*8	:: dx,dy,Lx,Ly,posx,posy
	real*8	:: Aii,Cx,Cy,D
	real*8, dimension(:), allocatable :: U, Uold, RHS, Uschwarz
	real*8, dimension(:), allocatable :: Utop, Ubottom
	integer	:: i,j,k
	character*10 :: mode
	integer :: send_bottom, recv_bottom, send_top, recv_top

	call MPI_Init(IERROR)
	call MPI_Comm_rank(MPI_COMM_WORLD, myrank, IERROR)
	call MPI_Comm_size(MPI_COMM_WORLD, wsize, IERROR)
	 
	max_iter = 100000 
	j=0
	if(iargc()/=0) then
		call getarg(1, mode)
		mode = TRIM(mode)
	else
		mode = "0"
	endif

        Ny = 10!*Np
        Nx = 10 
        R  = 1 !Recouvrement en Nbr de lignes (1 est le minimum)
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
	
	!calcul des zones pour chaque processus
	Nb_ligne_proc = Ny/wsize
	reste = mod(Ny,wsize)
	
    dx = Lx/(1+Nx)
    dy = Ly/(1+Ny)
	Aii = 2.d0*D/(dx*dx)+2.d0*D/(dy*dy) ! Terme diagonal de la matrice
	Cx  = -1.d0*D/(dx*dx)	! Terme extra-diagonale proche
	Cy  = -1.d0*D/(dy*dy)	! Terme extra-diagonale eloigne
    N  = Nx*Ny
    
    if(myrank < reste) then
    ibottom = myrank*(Nb_ligne_proc+1)*Nx+1
    else
    ibottom = (myrank+reste)*(Nb_ligne_proc)*Nx+1
    endif
    itop = ibottom + (Nb_ligne_proc-1)*Nx
!~     if(myrank < reste) then
!~    	itop = Nx*(Nb_ligne_proc)+1
!~    	N = Nx*(Nb_ligne_proc+1+R)+1
!~     else
!~    	itop = Nx*(Nb_ligne_proc-1)+1
!~ 	N = Nx*(Nb_ligne_proc+R)+1
!~     endif
!~ 	ibottom = 1
	
!~ 	if(reste < myrank) then
!~ 	itop = itop+Nx
!~ 	N = Nx*(Nb_ligne_proc+1+R)
!~ 	else
!~ 	N = Nx*(Nb_ligne_proc+R)
!~ 	endif
	
	 write(*,*) 'ibottom,', ibottom,'itop', itop+Nx-1

    ! Allocation dynamique pour chaque proc en tenant compte du recouvrement
    ALLOCATE(U(1:N));ALLOCATE(RHS(1:N));ALLOCATE(Uold(1:N));
	! Allocation dynamique pour la communication
	ALLOCATE(Utop(1:Nx));ALLOCATE(Ubottom(1:Nx))

	! Initialisation du probleme
    U(1:N)    = 0.0d0
	Utop(1:Nx)  = 0.0d0
	Ubottom(1:Nx)  = 0.0d0
	error_calcul=1

!	schw_iter = 0
	ERR_glob = 1.d3
	! Remplissage du second membre incluant les conditions de Bords Dirichlet
	call second_membre(RHS,U,Utop,Ubottom,77,77,1,N,Nx,dx,dy,Cx,Cy)
	!Preparation de requete persistante
	if (myrank == 0)then
		call MPI_Send_init(U(itop:itop+Nx-1), Nx, MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, send_top, IERROR)
		call MPI_Recv_init(Utop(1:Nx), Nx,MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, recv_top, IERROR)
	elseif (myrank == (wsize - 1)) then
		call MPI_Send_init(U(ibottom:ibottom+Nx-1), Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, send_bottom, IERROR)
		call MPI_Recv_init(Ubottom(1:Nx),Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, recv_bottom, IERROR)
	else
		call MPI_Send_init(U(itop:itop+Nx-1), Nx, MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, send_top, IERROR)
		call MPI_Recv_init(Utop(1:Nx),Nx, MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, recv_top, IERROR)
		call MPI_Send_init(U(ibottom:ibottom +Nx-1), Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, send_bottom, IERROR)
		call MPI_Recv_init(Ubottom(1:Nx), Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, recv_bottom, IERROR)
	endif
	
	! Initialisation des bords
	! Communication entre les vecteurs
		if (myrank == 0)then
		call MPI_Start(send_top, IERROR)
		call MPI_Start(recv_top, IERROR)
		elseif (myrank == (wsize - 1)) then
		call MPI_Start(send_bottom, IERROR)
		call MPI_Start(recv_bottom, IERROR)
		else
		call MPI_Start(send_top, IERROR)
		call MPI_Start(recv_top, IERROR)
		call MPI_Start(send_bottom, IERROR)
		call MPI_Start(recv_bottom, IERROR)
		endif
	
	if (myrank == 0)then
		call MPI_Wait(recv_top, IERROR)
		elseif (myrank == (wsize - 1)) then
		call MPI_Wait(recv_bottom, IERROR)
		else
		call MPI_Wait(recv_top, IERROR)
		call MPI_Wait(recv_bottom, IERROR)
		endif
        error_calcul=1	
	!write(*,*) 'debut boucle err:', error_calcul,' j:', j

	! Boucle tant que non convergence
	Do while ((error_calcul > eps) .AND. (j < max_iter))
	write(*,*) 'debut boucle err:', error_calcul,' j:', j
		! Communication entre les vecteurs
		if (myrank == 0)then
		call MPI_Start(send_top, IERROR)
		call MPI_Start(recv_top, IERROR)
		elseif (myrank == (wsize - 1)) then
		call MPI_Start(send_bottom, IERROR)
		call MPI_Start(recv_bottom, IERROR)
		else
		call MPI_Start(send_top, IERROR)
		call MPI_Start(recv_top, IERROR)
		call MPI_Start(send_bottom, IERROR)
		call MPI_Start(recv_bottom, IERROR)
		endif
		
!~         Uold(1:N)=U(1:N)
	Uold(ibottom:itop+Nx-1)=U(ibottom:itop+Nx-1)

		! Resolution du systeme lineaire (une iteration)
			if(mode == "0")then
				call grad(Aii,Cx,Cy,Nx,1,N,RHS,U,l)
			else if(mode == "1") then
				call gauss(Aii,Cx,Cy,Nx,ibottom,itop+Nx-1,RHS,U, Uold)
			else
				call jacobi(Aii,Cx,Cy,Nx,ibottom,itop+Nx-1,RHS,U,Uold,l)
			endif
			
		! Attente de la fin des envois
		if (myrank == 0)then
			call MPI_Wait(send_top, IERROR)
		elseif (myrank == (wsize - 1)) then
			call MPI_Wait(send_bottom, IERROR)
		else
			call MPI_Wait(send_top, IERROR)
			call MPI_Wait(send_bottom, IERROR)
		endif
		
		!calcul de l erreur
		Uold(ibottom:itop+Nx-1)=U(ibottom:itop+Nx-1)-Uold(ibottom:itop+Nx-1)
        error_calcul = sqrt(dot_product(Uold(ibottom:itop+Nx-1),Uold(ibottom:itop+Nx-1)))
        write(*,*)'avant allreduce err :', error_calcul		
		
		! Test de verification de convergence (MPI_all_reduce)
		call MPI_Allreduce(MPI_IN_PLACE, error_calcul,1,  MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
		write(*,*)'apres allreduce err :', error_calcul	
		j=j+1
		
		! Attente de la fin des receptions
		if (myrank == 0)then
			call MPI_Wait(recv_top, IERROR)
		elseif (myrank == (wsize - 1)) then
			call MPI_Wait(recv_bottom, IERROR)
		else
			call MPI_Wait(recv_top, IERROR)
			call MPI_Wait(recv_bottom, IERROR)
		endif
		write(*,*)'fin boucle :', error_calcul, 'j:', j	
	ENDDO
		write(*,*) 'fin , convergence en ', j, 'iterations', error_calcul
!~ 		call jacobi(Aii,Cx,Cy,Nx,1,N,RHS,U,Uold,l)

	call wrisol( U,Nx,Ny,dx,dy,77,1,N )
	
	if (myrank == 0)then
		call MPI_Request_free(send_top, IERROR)
		call MPI_Request_free(recv_top, IERROR)
	elseif (myrank == (wsize - 1)) then
		call MPI_Request_free(send_bottom, IERROR)
		call MPI_Request_free(recv_bottom, IERROR)
	else
		call MPI_Request_free(send_bottom, IERROR)
		call MPI_Request_free(recv_bottom, IERROR)
		call MPI_Request_free(send_top, IERROR)
		call MPI_Request_free(recv_top, IERROR)
	endif

	call MPI_FINALIZE(IERROR)

	DEALLOCATE(U,RHS,Uold)
	DEALLOCATE(Utop,Ubottom)
   end program schwarz_additif
