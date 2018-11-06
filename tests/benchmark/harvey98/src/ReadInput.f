       SUBROUTINE ReadInput(natom,nx,Ea,Ga,Eb,Gb)
       implicit none
       
       INTEGER i, j, k, kk, natom, nx, error
       DOUBLE PRECISION Ea, Eb, Ga(Nx), Gb(Nx)
       CHARACTER*60 dummy
       DOUBLE PRECISION tmpg, bohr
       PARAMETER (bohr=0.529177d0)

       OPEN(UNIT=8,FILE="ab_initio")
       
       read (8,*) dummy
       read (UNIT=8,FMT=*,IOSTAT=error) Ea
       IF (error .ne. 0) THEN
            error = -5
            CALL error_handling(error)
       END IF
       read (8,*) dummy
       DO I = 1, natom
           k = 3 * (i - 1) + 1
           READ (UNIT=8,FMT=*,IOSTAT=error) kk, (Ga(J),J=k,k+2)
           IF (error .ne. 0) THEN
                error = -5
                CALL error_handling(error)
           END IF
       END DO
       read (8,*) dummy
C The gradient output by most programs is in fact the FORCE - convert here to gradient.
C Also, convert from Hartree/Bohr to Hartree/Ang
       DO i = 1, nx
           Ga(i) = -Ga(i) / bohr
       end do

       read (UNIT=8,FMT=*,IOSTAT=error) Eb
       read (8,*) dummy
       IF (error .ne. 0) THEN
            error = -5
            CALL error_handling(error)
       END IF
       DO I = 1, natom
           k = 3 * (i - 1) + 1
           READ (UNIT=8,FMT=*,IOSTAT=error) kk, (Gb(J),J=k,k+2)
           IF (error .ne. 0) THEN
                error = -5
                CALL error_handling(error)
           END IF
       END DO
       DO i = 1, nx
           Gb(i) = -Gb(i) / bohr
       end do

       CLOSE(8)
       return

       END
