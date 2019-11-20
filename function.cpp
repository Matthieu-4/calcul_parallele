#include <mpi>
#include <math>
#include <iostream>
#include <fstream>

int Status[MPI_STATUS_SIZE] = {};
int PR = 4;
int tag = 100;
double *U, *U0;
int cond_init;
double sx,sy,dx,dy,dt,Lx,Ly,D,eps1,t,tf;
int i,me,Np,statinfo,i1,iN,j1,jN,Ny,Nx,n,q,k,kmax;

void Init(void){
  int i1, iN, q;
  int* tab[9] = {&Nx, &Ny, &Lx, &Ly, &D, &dt, &eps1, &kmax, &tf};
  string line;
  ifstream myfile ("data.txt");
  int i = 0;
  if (myfile.is_open())
  {
    for(i = 0; i < 9; i++){
      *tab[i] = std::atof(getline(myfile,line));
    }
    myfile.close();
  } else {
    cout << "Error when openning 'data.txt' file\n";
    return;
  }
  file.close();

  dx = Lx/(Nx+1);
  dy = Ly/(Ny+1);
  sx = D*dt/(dx*dx);
  sy = D*dt/(dy*dy);

  charge2(Ny,Np,me,&i1,&iN,&q);

  U = calloc((iN-i1) * sizeof(double));
  U0 = calloc((iN-i1) * sizeof(double));
  for(i = 0; i < iN-i1; i++){
    U0[i] = 1.;
  }
}



double f(x, y, t){
  if(cond_init == 3){
    return std::exp(-((x - Lx/2) * (x - Lx/2))) * std::exp(-((y - Ly/2) * (y - Ly/2))) * std::cos(t*pi/2);
  }else if(cond_init == 2){
    return std::sin(x) + std::cos(y);
  }else{
    return 2 * (y - y * y + x - x * x);
  }
}


double g(x, y, t){
  (void) t;
  if(cond_init == 2){
    return std::sin(x) + std::cos(y);
  }else{
    return 2 * (y - y * y + x - x * x);
  }
}

double h(x, y, t){
  if(cond_init == 3){
    return 1.;
  }else if(cond_init == 2){
    return std::sin(x) + std::cos(y);
  }else{
    return 2 * (y - y * y + x - x * x);
  }
}



void ProdMatVect(double D1[],
                 double D2_m[],
                 double D2_p[],
                 double D3_m[],
                 double D3_p[],
                 double x[],
                 double y[],
                 int i1,
                 int iN){


  double* x1 = (double*)malloc(sizeof(double) * (iN + Nx - i1));
  double* x2 = (double*)malloc(sizeof(double) * (iN + Nx - i1));
  double* x3 = (double*)malloc(sizeof(double) * (iN + Nx - i1));

  int i = 0;
  if(i1 == 1){
    for(i = 0; i < iN - i1; i++){
      x1[i] = x[i];
    }
    charge2(Ny, Np, me+1, &j1, &jN, &q);
    //MPI_RECV(x1(j1:j1+Nx-1),Nx,MPI_REAL,1,tag,MPI_COMM_WORLD,status,statinfo)
    //MPI_SEND(x1(iN-Nx+1:iN),Nx,MPI_REAL,1,tag,MPI_COMM_WORLD,statinfo)

    y[0] = D1[0] * x1[0] + D2_p[0] * x1[1] + D3_p[0] * x1[Nx]
    for(i = 1; i < iN; i++){
      y[i] = D2_m[i]*x1[i-1]+D1[i]*x1[i]+D2_p[i]*x1[i+1]+D3_p[i]*x1[i+Nx]
    }
    for(i = Nx; i < iN; i++){
      y[i] = D3_m[i] * x1[i-Nx] + D2_m[i] * x1[i-1] + D1[i] * x1[i] + D2_p[i] * x1[i+1] + D3_p[i] * x1[i+Nx]
    }
  } else if (iN == Nx*Ny){
    double *x1 = (double*)malloc(sizeof(double) * (iN + Nx - i1));
    for(i = 0; i < iN - i1; i++){
      x3[i] = x[i];
    }
    //MPI_SEND(x3(i1:i1+Nx-1),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
    charge2(Ny, Np, me-1, &j1, &jN, &q);
    //MPI_RECV(x3(jN-Nx+1:jN),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,status,statinfo)
    for(i = 0; i < iN - Nx - i1; i++){
      y[i] = D3_m[i]*x3[i-Nx]+D2_m[i]*x3[i-1]+D1[i]*x3[i]+D2_p[i]*x3[i+1]+D3_p[i]*x3[i+Nx]
    }
    for(i = iN - Nx + 1 - i1; i < iN - i1 - 1; i++){
      y[i] = D3_m[i]*x3[i-Nx]+D2_m[i]*x3[i-1]+D1[i]*x3[i]+D2_p[i]*x3[i+1]
    }
    y[Nx*Ny] = D3_m[Nx*Ny]*x3[Nx*Ny-Nx]+D2_m[Nx*Ny]*x3[Nx*Ny-1]+D1[Nx*Ny]*x3[Nx*Ny]
  } else {
    for(i = 0; i < iN - i1; i++){
      x2[i] = x[i];
    }

    Call MPI_SEND(x2(i1:i1+Nx-1),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
    Call charge2(Ny,Np,me-1,j1,jN,q)
    Call MPI_RECV(x2(jN-Nx+1:jN),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,status,statinfo)
    Call charge2(Ny,Np,me+1,j1,jN,q)
    Call MPI_RECV(x2(j1:j1+Nx-1),Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,status,statinfo)
    Call MPI_SEND(x2(iN-Nx+1:iN),Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)

    Do i=i1,iN
       y(i) = D3_m(i)*x2(i-Nx)+D2_m(i)*x2(i-1)+D1(i)*x2(i)+D2_p(i)*x2(i+1)+D3_p(i)*x2(i+Nx)
    End Do
  }
}

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
