Module matrice

  Use Function
  Implicit None


Contains

  Subroutine MatriceDF(D1,D2_m,D2_p,D3_m,D3_p,D,sx,sy,i1,iN) !Construction de la matrice DF (cas linéaire)
    Real(PR),Dimension(i1:iN),Intent(Inout) :: D1
    Real(PR),Dimension(i1:iN),Intent(Inout) :: D2_m,D2_p
    Real(PR),Dimension(i1:iN),Intent(Inout) :: D3_m,D3_p
    Real(PR),Intent(In) :: D,sx,sy
    Integer, Intent(In) :: i1,iN
    Integer :: k

    D1 = 1._PR + 2*D*(sx+sy)
    D2_m = -D*sx
    D2_p = -D*sx
    D3_m = -D*sy
    D3_p = -D*sy

    If (iN == Nx*Ny) Then
       Do k=Nx*Ny-Nx+1,Nx*Ny
          D3_p(k) = 0._PR
       End Do
    End If

    If (i1 == 1) Then
       Do k=1,Nx
          D3_m(k) = 0._PR
       End Do
    End If

    Do k=i1,iN,Nx
       D2_p(k+Nx-1) = 0._PR
       D2_m(k) = 0._PR
    End Do

  End Subroutine MatriceDF


  Subroutine sec_membre(Nx,Ny,dx,dy,dt,Lx,Ly,D,Fx,t,i1,iN)
    Integer,Intent(In) :: Nx,Ny,i1,iN
    Real(PR),Intent(In) :: dx,dy,dt,Lx,Ly,D,t
    Real(PR),Dimension(i1:iN),Intent(Out) :: Fx
    Integer :: i,j,k
    Real(PR) :: sx,sy,A,B,C

    sx = dt/(dx**2)
    sy = dt/(dy**2)
    A = 1+2*D*(sx+sy)
    B = -D*sx
    C = -D*sy

    

    If (i1 == 1) Then              ! Premier proc
       Fx(i1) = dt*f(dx,dy,t) - C*g(dx,0._PR,t)-B*h(0._PR,dy,t)
       Do k=2,Nx-1
          Fx(k) = dt*f(Reste(k,Nx)*dx,dy,t) - C*g(Reste(k,Nx)*dx,0._PR,t)
       End Do
       Fx(Nx) = dt*f(Reste(k,Nx)*dx,dy,t) - C*g(Lx-dx,0._PR,t)-B*h(Lx,dy,t)

       Do k = Nx+1,iN                ! k(i,j) = i + Nx*(j-1)
          i = Reste(k,Nx)            ! i(k) = reste de k/Nx (+ voir fonction Reste)
          j = (k-1)/Nx + 1           ! j(k) = (quotient de k-1 divis� par Nx) + 1
          Fx(k) = dt*f(i*dx,j*dy,t)
          If (Mod(k,Nx) == 1) Then
             Fx(k) = Fx(k)-B*h(0._PR,j*dy,t)
          Else If (Mod(k,Nx) == 0) Then
             Fx(k) = Fx(k)-B*h(Lx,j*dy,t)
          End If
       End Do

    Else If (iN == Nx*Ny) Then               ! Dernier proc
       Do k = i1,iN - Nx
          i = Reste(k,Nx)
          j = (k-1)/Nx + 1
          Fx(k) = dt*f(i*dx,j*dy,t)
          If (Mod(k,Nx) == 1) Then
             Fx(k) = Fx(k)-B*h(0._PR,j*dy,t)
          Else If (Mod(k,Nx) == 0) Then
             Fx(k) = Fx(k)-B*h(Lx,j*dy,t)
          End If
       End Do

       Fx(iN - Nx + 1) = dt*f(dx,Ly-dy,t)-C*g(dx,Ly,t)-B*h(0._PR,Ly-dy,t)
       Do k=iN - Nx + 2, iN-1
          Fx(k) = dt*f(Reste(k,Nx)*dx,Ly-dy,t) - C*g(Reste(k,Nx)*dx,Ly,t)
       End Do
       Fx(iN) = dt*f(Lx-dx,Ly-dy,t)-C*g(Lx-dx,Ly,t)-B*h(Lx,Ly-dy,t)



    Else                 ! Les autres procs
       Do k = i1,iN
          i = Reste(k,Nx)
          j = (k-1)/Nx + 1

          Fx(k) = dt*f(i*dx,j*dy,t)
          If (Mod(k,Nx) == 1) Then
             Fx(k) = Fx(k)-B*h(0._PR,j*dy,t)
          Else If (Mod(k,Nx) == 0) Then
             Fx(k) = Fx(k)-B*h(Lx,j*dy,t)
          End If
       End Do
    End If

  End Subroutine sec_membre




End Module matrice
