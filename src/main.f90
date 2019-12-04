Program Parallel

  Use function
  Use matrice
  Implicit None

  ! Definition des variables
  Real(PR), dimension(:), Allocatable :: Fx,D1,D2_m,D2_p,D3_m,D3_p
  character(len=50) :: file_name, Me_string
  Real(PR):: t1,t2,temps

  ! Mise en place de l'environnement parallele
  Call MPI_INIT(statinfo)
  Call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

  t=0._PR

  ! Choix condition initiale

  call Init()

  if (me == 0) then
     print *,
     print *, "Choisissez le jeu de conditions initiales (1,2 ou 3)"
     print *, "1 : f = 2(y - y^2 + x - x^2) ; g = 0 : h = 0"
     print *, "2 : f = sin(x) + cos(y) ; g = sin(x) + cos(y) : h = sin(x) + cos(y)"
     print *, "3 : f = exp(-(x - Lx/2)^2).exp(-(y - Ly/2)^2).cos(t.pi/2) ; g = 0 : h = 1"

  read *, cond_init

     if (cond_init > 3 .or. cond_init < 1) then
        print*, "Attention ! Le set (1) a ete pris par default"
     end if
  end if

  ! Partage du choix des conditions initiales aux autres processeurs
  call MPI_BCAST(cond_init,1,MPI_INTEGER,0,MPI_COMM_WORLD,statinfo)

!!$  ! On mesure le temps de calcul par proc, premier compteur t1
!!$  call CPU_TIME(t1)

  ! Répartition des procs
  call charge2(Ny,Np,me,i1,iN,q)


  Allocate(Fx(i1:iN),D1(i1:iN),D2_m(i1:iN),D2_p(i1:iN),D3_m(i1:iN),D3_p(i1:iN))
  Fx = 0._PR
  Call MatriceDF(D1,D2_m,D2_p,D3_m,D3_p,D,sx,sy,i1,iN)

  k = 0
  !Boucle en temps
  Do while (t<tf)
     Call sec_membre(Nx,Ny,dx,dy,dt,Lx,Ly,D,Fx,t,i1,iN)
     Call grad_conj(D1,D2_m,D2_p,D3_m,D3_p,U,Fx+U0,i1,iN)
     U0 = U
     t = t +dt
     k = k+1
     if (mod(k,100) == 0) then
        print*, "nb_iter : ", k
     end if
  end Do


  ! Ecriture de la solution sur des fichiers .dat pour chaque proc
  Write( Me_string, '(i10)' )  Me
  file_name = 'sol00' // trim(adjustl(Me_string)) // '.dat'
  Open(10+Me, File=trim(file_name))
  Do k=i1,iN
     Write(10+Me,*) Reste(k,Nx)*dx, ((k-1)/Nx+1)*dy, U(k), k
  End Do
  close(10+Me)

!!$  ! Deuxième compteur en temps
!!$  call CPU_TIME(t2)
!!$
!!$  ! Temps de calcul final
!!$  temps=t2-t1
!!$  print*,'Le temps de calcul du proc',me,' pour résoudre le problème est :',temps,'s'

  call MPI_FINALIZE(statinfo)


Contains




End Program Parallel
