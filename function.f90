Module Function

  Implicit None
  Include 'mpif.h'

  Integer,Dimension(MPI_STATUS_SIZE) :: Status
  Integer,Parameter :: PR=4,  tag=100
  Real(PR),Parameter :: pi=4*Atan(1.)
  Real, Dimension(:), Allocatable :: U, U0
  Integer :: cond_init
  Real(Kind=PR) :: sx,sy,dx,dy,dt,Lx,Ly,D,eps1,t,tf
  Integer :: i,me,Np,statinfo,i1,iN,j1,jN,Ny,Nx,n,q,k,kmax


Contains

  Subroutine Init()     ! Lecture du fichier permettant d'initialiser le syst�me
    Integer :: i1,iN,q

    Open(1,file = "data.txt")
    Read(1,*)Nx
    Read(1,*)Ny
    Read(1,*)Lx
    Read(1,*)Ly
    Read(1,*)D
    Read(1,*)dt
    Read(1,*)eps1
    Read(1,*)kmax
    Read(1,*)tf

    dx = Lx/(Nx+1)
    dy = Ly/(Ny+1)

    sx = D*dt/(dx*dx)
    sy = D*dt/(dy*dy)

    Call charge2(Ny,Np,me,i1,iN,q)


    Allocate(U(i1:iN))
    Allocate(U0(i1:iN))
    U0 = 1._PR
    U = 0._PR

    Close(1)

  End Subroutine Init




  Function f(x,y,t)Result(f_xy)         !!! f(x,y,t)
    Real,Intent(in)::x,y,t
    Real :: f_xy

    Select Case (cond_init)
    Case(1)
       f_xy = 2*(y - y*y + x - x*x)
    Case(2)
       f_xy = Sin(x) + Cos(y)
    Case(3)
       f_xy = Exp(-((x - Lx/2)**2))*Exp(-((y - Ly/2)**2))*Cos(t*pi/2)
    Case default
       f_xy = 2*(y - y*y + x - x*x)
    End Select

  End Function f


  Function g(x,y,t)Result(g_xy)         !!! g(x,y,t)
    Real,Intent(in)::x,y,t
    Real :: g_xy

    Select Case (cond_init)
    Case(1)
       g_xy = 0.
    Case(2)
       g_xy = Sin(x) + Cos(y)
    Case(3)
       g_xy = 0.
    Case default
       g_xy = 0.
    End Select

  End Function g


  Function h(x,y,t)Result(h_xy)         !!! h(x,y,t)
    Real,Intent(in)::x,y,t
    Real :: h_xy

    Select Case (cond_init)
    Case(1)
       h_xy = 0.
    Case(2)
       h_xy = Sin(x) + Cos(y)
    Case(3)
       h_xy = 1.
    Case default
       h_xy = 0.
    End Select

  End Function h




  Subroutine ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,x,y,i1,iN)     ! Produit Mat/Vect en creux y=Ax, A=(D1,D2,D3)
    Integer,Intent(In) :: i1, iN
    Real(PR),Dimension(i1:iN),Intent(In) :: D1
    Real(PR),Dimension(i1:iN),Intent(In) :: D2_p,D2_m
    Real(PR),Dimension(i1:iN),Intent(In) :: D3_p,D3_m
    Real(PR),Dimension(i1:iN),Intent(In) :: x
    Real(PR),Dimension(i1:iN+Nx) :: x1              ! Stock pour le premier proc
    Real(PR),Dimension(i1-Nx:iN+Nx) :: x2           ! Stock pour les autres procs
    Real(PR),Dimension(i1-Nx:iN) :: x3              ! Stock pour le dernier proc
    Real(PR),Dimension(i1:iN), Intent(Out) :: y
    Integer :: i


    If (i1==1) Then          ! Premier proc
       do i = i1,iN
          x1(i) = x(i)
       end do

       Call charge2(Ny,Np,me+1,j1,jN,q)
       Call MPI_RECV(x1(j1:j1+Nx-1),Nx,MPI_REAL,1,tag,MPI_COMM_WORLD,status,statinfo)
       Call MPI_SEND(x1(iN-Nx+1:iN),Nx,MPI_REAL,1,tag,MPI_COMM_WORLD,statinfo)

       y(i1) = D1(i1)*x1(i1)+D2_p(i1)*x1(i1+1)+D3_p(i1)*x1(i1+Nx)     ! 1ere ligne (1er �lement)
       Do i=2,iN                                                      ! 1ere ligne
          y(i) = D2_m(i)*x1(i-1)+D1(i)*x1(i)+D2_p(i)*x1(i+1)+D3_p(i)*x1(i+Nx)
       End Do
       Do i = Nx+1,iN
          y(i) = D3_m(i)*x1(i-Nx)+D2_m(i)*x1(i-1)+D1(i)*x1(i)+D2_p(i)*x1(i+1)+D3_p(i)*x1(i+Nx)
       End Do


    Else If (iN==Nx*Ny) Then        ! Dernier proc
       do i = i1,iN
          x3(i) = x(i)
       end do

       Call MPI_SEND(x3(i1:i1+Nx-1),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
       Call charge2(Ny,Np,me-1,j1,jN,q)
       Call MPI_RECV(x3(jN-Nx+1:jN),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,status,statinfo)

       Do i = i1,iN-Nx
          y(i) = D3_m(i)*x3(i-Nx)+D2_m(i)*x3(i-1)+D1(i)*x3(i)+D2_p(i)*x3(i+1)+D3_p(i)*x3(i+Nx)
       End Do
       Do i= iN - Nx + 1,iN - 1
          y(i) = D3_m(i)*x3(i-Nx)+D2_m(i)*x3(i-1)+D1(i)*x3(i)+D2_p(i)*x3(i+1)
       End Do
       y(Nx*Ny) = D3_m(Nx*Ny)*x3(Nx*Ny-Nx)+D2_m(Nx*Ny)*x3(Nx*Ny-1)+D1(Nx*Ny)*x3(Nx*Ny)


    Else                   ! autres procs
       do i = i1,iN
          x2(i) = x(i)
       end do

       Call MPI_SEND(x2(i1:i1+Nx-1),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
       Call charge2(Ny,Np,me-1,j1,jN,q)
       Call MPI_RECV(x2(jN-Nx+1:jN),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,status,statinfo)
       Call charge2(Ny,Np,me+1,j1,jN,q)
       Call MPI_RECV(x2(j1:j1+Nx-1),Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,status,statinfo)
       Call MPI_SEND(x2(iN-Nx+1:iN),Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)

       Do i=i1,iN
          y(i) = D3_m(i)*x2(i-Nx)+D2_m(i)*x2(i-1)+D1(i)*x2(i)+D2_p(i)*x2(i+1)+D3_p(i)*x2(i+Nx)
       End Do

    end If

  End Subroutine ProdMatVect





  Subroutine Charge2(n,Np,me,i1,iN,q2)     ! Charge adapt�

    Integer, Intent(In) :: n,Np,me
    Integer, Intent(Out) :: i1,iN,q2
    Integer :: q,r

    q = n/Np
    r = Mod(n,Np)

    If (r==0) Then
       i1 = 1 + q*me
       iN = q*(me+1)

    Else
       If (me<r) Then
          i1=me*(q+1)+1
          iN=(me+1)*(q+1)
       Else
          i1=1+ me*q +r
          iN=i1+q-1
       End If
    End If

    ! q est le nombre de lignes connues par le proc
    q2 = iN-i1 + 1

    ! On travaille sur n = nx*ny éléments
    i1 = (i1-1)*Nx +1
    iN = iN*Nx

  End Subroutine Charge2

  subroutine charge(n,Np,me,i1,iN)
    integer,intent(in)::n,Np,me
    integer,intent(out)::i1,iN
    integer::q,r
    q = n/Np
    r = mod(n,Np)
    if (me<r) then
       i1 = me*(q+1) + 1
       iN = (me+1)*(q+1)
    else
       i1 = 1 + r + me*q
       iN = i1 + q - 1
    end if

  end subroutine charge



  Subroutine grad_conj(D1,D2_m,D2_p,D3_m,D3_p,x,b,i1,iN)   ! Gradient conjugué creux pour Ax=b, A=(D1,D2,D3)
    Real(PR),Dimension(i1:iN),Intent(In) :: D1,D2_m,D2_p,D3_m,D3_p
    Real(PR),Dimension(i1:iN),Intent(Out) :: x
    Real(PR),Dimension(i1:iN),Intent(In) :: b
    Real(PR),Dimension(i1:iN) :: r,r2,z,p
    Real(PR),Dimension(i1:iN) :: y
    Real(PR) :: alpha,beta,gamma
    Integer,Intent(In) :: i1,iN
    Integer :: k

    x(i1:iN) = 293._PR
    Call ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,x,y,i1,iN)
    r = b-y
    p = r
    beta = Sqrt(prodscal(r,r,i1,iN))
    k=0
    Do While ((beta>eps1) .And. (k<=kmax))
       Call ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,p,z,i1,iN)
       alpha = prodscal(r,r,i1,iN)/(prodscal(z,p,i1,iN))
       x = x + alpha*p
       r2 = r
       r = r - alpha*z
       gamma = prodscal(r,r,i1,iN)/prodscal(r2,r2,i1,iN)
       p = r + gamma*p
       beta = Sqrt(prodscal(r2,r2,i1,iN))
       k = k+1


    End Do

  End Subroutine grad_conj



  Function prodscal(X,Y,i1,iN)Result(Somme)                 !!! Produit scalaire
    Integer, Intent(In) :: i1,iN
    Real(Kind=PR), Dimension(i1:iN), Intent(In) :: X,Y
    Real(Kind=PR) :: Z,Somme
    Integer :: k,statinfo

    Somme=0._PR
    Z=0._PR


    Do k=i1,iN
       Z = Z + X(k)*Y(k)
    End Do

    Call MPI_ALLREDUCE(Z,Somme,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)

  End Function prodscal

  Function Reste(k,Nx)Result(r)
    Integer,intent(in) :: k,Nx
    integer :: r
    if (modulo(k,Nx) == 0)then
       r = Nx
    else
       r = modulo(k,Nx)
    end if
  end Function Reste




End Module Function
