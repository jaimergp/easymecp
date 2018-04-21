       SUBROUTINE Initialize(N,H,Ea,Eb,Ga,Gb)
       implicit none
       
       INTEGER i, j, n
       DOUBLE PRECISION H(N,N), Ea, Eb, Ga(N), Gb(N)
       
       Ea = 0.d0
       Eb = 0.d0
       do i = 1, n
           Ga(i) = 0.d0
           Gb(i) = 0.d0
           H(i,i) = .7d0
C Inverse Hessian values of .7 correspond to Hessians of 1.4 Hartree/Angstrom**2 - about right
           do j = i + 1, n
               H(i,j) = 0.d0
               H(j,i) = 0.d0
           end do
       end do
       return

       END
