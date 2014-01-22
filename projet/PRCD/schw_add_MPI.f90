program schwarz_additif
!
! Resolution du placien par l'algorithme de schwarz-additif
!   ∆U = f
!
        USE UTILS
    implicit none
    include 'mpif.h'

    ! Declaration des variables pour schwarz
    double precision, parameter  :: eps = 1.d-12
    double precision  :: error_calcul ,ERR_glob

    ! Decalaration des variables MPI et charges
    integer :: itop,ibottom  !indice de localisation des debuts de msg a envoyer

    ! Declaration des variables numeriques
    integer :: l,Max_l, Nb_elem, aux, aux2, aux3
    integer :: R, modif_d, modif_f, idebut, ifin
    double precision  :: t_debut,t_fin,t_max

    ! Declaration des variables du pb a resoudre
    integer :: Nx,Ny,N, myrank, wsize, IERROR, Nb_ligne_proc, max_iter, reste
    double precision  :: dx,dy,Lx,Ly,posx,posy
    double precision  :: Aii,Cx,Cy,D
    double precision, dimension(:), allocatable :: U, Uold,Uold2, RHS, Uschwarz
    double precision, dimension(:), allocatable :: Utop, Ubottom
    integer :: i,j,k
    character(len=10) :: mode
    integer :: send_bottom, recv_bottom, send_top, recv_top

    call MPI_Init(IERROR)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, IERROR)
    call MPI_Comm_size(MPI_COMM_WORLD, wsize, IERROR)

    max_iter = 100000
    j=0

    if(COMMAND_ARGUMENT_COUNT()/=0) then
        call GET_COMMAND_ARGUMENT(1, mode)
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

!~  write(*,*) 'nb lignes ', Nb_ligne_proc, ' reste ', reste

    dx = Lx/(1+Nx)
    dy = Ly/(1+Ny)
    Aii = 2.d0*D/(dx*dx)+2.d0*D/(dy*dy) ! Terme diagonal de la matrice
    Cx  = -1.d0*D/(dx*dx)   ! Terme extra-diagonale proche
    Cy  = -1.d0*D/(dy*dy)   ! Terme extra-diagonale eloigne
    N  = Nx*Ny
    Nb_elem = Nb_ligne_proc*Nx

    if(myrank < reste) then
        ibottom = myrank*(Nb_ligne_proc+1)*Nx+1
        itop = ibottom + Nb_ligne_proc*Nx
        Nb_elem = Nb_elem+Nx
    elseif (reste /= 0) then
        ibottom = (myrank+reste-1)*Nb_ligne_proc*Nx+1
        itop = ibottom + (Nb_ligne_proc-1)*Nx
    else
        ibottom = (myrank)*Nb_ligne_proc*Nx+1
        itop = ibottom + (Nb_ligne_proc-1)*Nx
    endif

!~      write(*,*) 'ibottom,', ibottom,'itop', itop+Nx-1

    ! Allocation dynamique pour chaque proc en tenant compte du recouvrement
    ALLOCATE(U(1:N));ALLOCATE(RHS(1:N));ALLOCATE(Uold(1:N));ALLOCATE(Uold2(1:N));

    ! Allocation dynamique pour la communication
    ALLOCATE(Utop(1:Nx));ALLOCATE(Ubottom(1:Nx))

    ! Initialisation du probleme
    U(1:N)    = 0.0d0
    Utop(1:Nx)  = 0.0d0
    Ubottom(1:Nx)  = 0.0d0
    error_calcul=1

    if(myrank ==0)then
        idebut = 1
        ifin = itop + R*Nx + Nx -1
    else if (myrank == wsize -1) then
        idebut = ibottom -R*Nx
        ifin = itop + Nx - 1
    else
        idebut = ibottom - R*Nx
        ifin = itop + R*Nx + Nx -1
    endif

!   schw_iter = 0
    ERR_glob = 1.d3
    !Preparation de requete persistante
    if (myrank == 0)then
!~         write(*,*) 'proc0 envoi top', itop-Nx, 'fin', itop-1, 'nb' ,Nx
        call MPI_Send_init(U(itop-Nx:itop-1), Nx, MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, send_top, IERROR)
        call MPI_Recv_init(Utop(1:Nx), Nx,MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, recv_top, IERROR)
    elseif (myrank == (wsize - 1)) then
!~         write(*,*) 'proc n-1 envoi bottom', ibottom+Nx, 'fin', ibottom+2*Nx-1, 'nb' ,Nx
        call MPI_Send_init(U(ibottom+Nx:ibottom+2*Nx-1), Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, send_bottom, IERROR)
        call MPI_Recv_init(Ubottom(1:Nx),Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, recv_bottom, IERROR)
    else
!~         write(*,*) 'proc', myrank, 'envoi top', itop-Nx, 'fin', itop-1, 'nb' ,Nx
!~         write(*,*) 'proc',myrank,'envoi bottom', ibottom+Nx, 'fin', ibottom+2*Nx-1, 'nb' ,Nx
        call MPI_Send_init(U(itop-Nx:itop), Nx, MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, send_top, IERROR)
        call MPI_Recv_init(Utop(1:Nx),Nx, MPI_DOUBLE_PRECISION, myrank+1, 1, MPI_COMM_WORLD, recv_top, IERROR)
        call MPI_Send_init(U(ibottom+Nx:ibottom+2*Nx-1), Nx, MPI_DOUBLE_PRECISION, myrank-1, 1, MPI_COMM_WORLD, send_bottom, IERROR)
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
        call MPI_Wait(recv_top, MPI_STATUS_IGNORE,IERROR)
        elseif (myrank == (wsize - 1)) then
        call MPI_Wait(recv_bottom, MPI_STATUS_IGNORE,IERROR)
    else
        call MPI_Wait(recv_top, MPI_STATUS_IGNORE,IERROR)
        call MPI_Wait(recv_bottom, MPI_STATUS_IGNORE,IERROR)
    endif
        write(*,*) 'proc', myrank, idebut, ifin

    ! Remplissage du second membre incluant les conditions de Bords Dirichlet
    call second_membre(RHS(idebut:ifin),U(idebut:ifin),Utop,Ubottom,myrank,wsize,idebut,ifin,Nx,dx,dy,Cx,Cy)

    ! Boucle tant que non convergence
    Do while ((error_calcul > eps) .AND. (j < max_iter))

        Uold2(idebut:ifin)=U(idebut:ifin)

        ! Resolution du systeme lineaire
        if(mode == "0")then
            call grad(Aii,Cx,Cy,Nx,idebut,ifin,RHS(idebut:ifin),U(idebut:ifin),l)
        else if(mode == "1") then
            call gauss(Aii,Cx,Cy,Nx,idebut,ifin,RHS,U, Uold)
        else
            call jacobi(Aii,Cx,Cy,Nx,idebut,ifin,RHS(idebut:ifin),U(idebut:ifin),Uold(idebut:ifin),k)
        endif

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

        !calcul de l erreur
        Uold2(idebut:ifin)=U(idebut:ifin)-Uold2(idebut:ifin)
        error_calcul = sqrt(dot_product(Uold2(idebut:ifin),Uold2(idebut:ifin)))

        ! Test de verification de convergence (MPI_all_reduce)
        call MPI_Allreduce(MPI_IN_PLACE, error_calcul,1,  MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)

        j=j+1

        ! Attente de la fin des receptions
        if (myrank == 0)then
            call MPI_Wait(recv_top, MPI_STATUS_IGNORE, IERROR)
        elseif (myrank == (wsize - 1)) then
            call MPI_Wait(recv_bottom, MPI_STATUS_IGNORE, IERROR)
        else
            call MPI_Wait(recv_top, MPI_STATUS_IGNORE, IERROR)
            call MPI_Wait(recv_bottom, MPI_STATUS_IGNORE, IERROR)
        endif

        !Modification des bords
        call second_membre(RHS(idebut:ifin),U(idebut:ifin),Utop,Ubottom,myrank,wsize,idebut,ifin,Nx,dx,dy,Cx,Cy)

        ! Attente de la fin des envois
        if (myrank == 0)then
            call MPI_Wait(send_top, MPI_STATUS_IGNORE, IERROR)
        elseif (myrank == (wsize - 1)) then
            call MPI_Wait(send_bottom, MPI_STATUS_IGNORE, IERROR)
        else
            call MPI_Wait(send_top, MPI_STATUS_IGNORE, IERROR)
            call MPI_Wait(send_bottom, MPI_STATUS_IGNORE, IERROR)
        endif
    ENDDO

    if(myrank==0)then
        write(*,*) 'fin , convergence en ', j, 'iterations', error_calcul
    endif

    ! Ecriture de chaque morceau dans son fichier (sol0??.dat)
    call wrisol(U(idebut:ifin),Nx,Ny,dx,dy,myrank,idebut,ifin )

	! Libération des requêtes préparées
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
	! Libération des matrices et vecteurs
    DEALLOCATE(U,RHS,Uold, Uold2)
    DEALLOCATE(Utop,Ubottom)
   end program schwarz_additif
