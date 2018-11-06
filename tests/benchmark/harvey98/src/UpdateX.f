       SUBROUTINE UpdateX(N,Nstep,FFile,X_1,X_2,X_3,G_1,G_2,HI_1,HI_2)
       implicit none
       
       ! Specially Adapted BFGS routine from Numerical Recipes
       
       integer i, j, k, n, Nstep, FFile
       double precision X_1(N), X_2(N), G_1(N), G_2(N), HI_1(N,N), X_3(N), HI_2(N,N)
C change fmc
       double precision stpmax, DelG(N), HDelG(N), ChgeX(N), DelX(N), w(N), 
     1 fac,fad, fae, sumdg, sumdx, stpl, lgstst, STPMX
       parameter (STPMX = 0.1d0)
       
       stpmax = STPMX * N
       IF ((Nstep .eq. 0) .and. (FFile .eq. 0)) THEN
           DO i = 1, n
               ChgeX(i) = -.7d0 * G_2(i)
               DO j = 1, n
                   HI_2(i,j) = HI_1(i,j)
               end do
           end do
       ELSE
           DO i = 1, n
               DelG(i) = G_2(i) - G_1(i)
               DelX(i) = X_2(i) - X_1(i)
           end do
           do i = 1, n
               HDelG(i) = 0.d0
               do j = 1, n
                   hdelg(i) = hdelg(i) + HI_1(i,j) * DelG(j)
               end do
           end do
           fac = 0.d0
           fae = 0.d0
           sumdg = 0.d0
           sumdx = 0.d0
           do i = 1, n
               fac = fac + delg(i) * delx(i)
               fae = fae + delg(i) * hdelg(i)
               sumdg = sumdg + delg(i)**2
               sumdx = sumdx + delx(i)**2
           end do
           fac = 1.d0 / fac
           fad = 1.d0 / fae
           do i = 1, n
               w(i) = fac * DelX(i) - fad * HDelG(i)
           end do
           DO I = 1, N
               do j = 1, n
                   HI_2(i,j) = HI_1(i,j) + fac * delx(i) * delx(j) -
C change fmc
C&                     fad * HDelG(I) * HDelG(j) + fae * w(i) * w(j)
     1                 fad * HDelG(I) * HDelG(j) + fae * w(i) * w(j)
               end do
           end do
           do i = 1, n
               ChgeX(i) = 0.d0
               do j = 1, n
                   ChgeX(i) = ChgeX(i) - HI_2(i,j) * G_2(j)
               end do
           end do
       END IF
       
       fac = 0.d0
       do i = 1, n
           fac = fac + ChgeX(i)**2
       end do
       stpl = SQRT(fac)
       IF (stpl .gt. stpmax) THEN
           do i = 1, n
               ChgeX(i) = ChgeX(i) / stpl * stpmax
           end do
       END IF
       lgstst = 0.d0
       do i = 1, n
            IF (ABS(ChgeX(i)) .gt. lgstst) lgstst = ABS(ChgeX(i))
       end do
       IF (lgstst .gt. STPMX) THEN
           do i = 1, n
               ChgeX(i) = ChgeX(i) / lgstst * STPMX
           end do
       END IF
       do i = 1, n
           X_3(i) = X_2(i) + ChgeX(i)
       end do
       
       END
