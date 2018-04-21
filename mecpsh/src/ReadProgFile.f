       SUBROUTINE ReadProgFile(natom,nx,AtNum,Nstep,FFile,X_1,X_2,HI,Ea,Eb,Ga,Gb,G)
       implicit none
       
       INTEGER natom, nx, AtNum(Natom), Nstep, i, j, k, error, FFile
       DOUBLE PRECISION X_1(Nx), X_2(Nx), HI(Nx,Nx), Ea, Eb, Ga(Nx), Gb(Nx), G(Nx)
       
       LOGICAL PGFok
       CHARACTER*80 dummy
       
       INQUIRE(FILE="ProgFile",EXIST=PGFok)
       IF (.not. PGFok) THEN
           error = -1
           CALL error_handling(error)
       END IF
       
       OPEN(UNIT=8,FILE="ProgFile")
       READ(8,*) dummy
       READ(8,*) dummy
       READ(8,*) i
       IF (i .ne. natom) then
           error = -4
           CALL error_handling(error)
       END IF
       READ(8,*) dummy
       READ(8,*) Nstep
       READ(8,*) dummy
       READ(8,*) FFile
       
       READ(8,*) dummy
       DO I = 1, natom
           k = 3 * (i - 1) + 1
           READ(8,*) AtNum(I), (X_2(J),J=k,k+2)
       END DO
       
       IF ((Nstep .eq. 0) .AND. (FFile .eq. 0)) THEN
C First step - the ProgFile is incomplete. To avoid problems,
C there are 'existence' tests in the next input lines so
C as to crash the program if Nstep > 0 yet the file is incomplete.
           CALL Initialize(Nx,HI,Ea,Eb,Ga,Gb)
           CLOSE(8)
           RETURN
       END IF
       
       READ(UNIT=8,FMT=*,IOSTAT=error) dummy
       IF (error .ne. 0) THEN
           error = -2
           CALL error_handling(error)
       END IF
       DO I = 1, Natom
           k = 3 *(I-1) + 1
           READ(UNIT=8,FMT=*,IOSTAT=error) (X_1(J),J=k,k+2)
C Geometries are in Angstrom - fine!
           IF (error .ne. 0) THEN
               error = -2
               CALL error_handling(error)
           END IF
       END DO
       READ(8,*) dummy
       READ(8,*) Ea
       READ(8,*) Eb
       READ(8,*) dummy
       DO I = 1, nx
           READ(8,*) Ga(I)
       END DO
       READ(8,*) dummy
       DO I = 1, nx
           READ(8,*) Gb(I)
       END DO
       READ(8,*) dummy
       DO I = 1, nx
           READ(8,*) G(I)
       END DO
       READ(8,*) dummy
       DO I = 1, nx
           DO J = 1, nx
               READ(8,*) HI(I,J)
           END DO
       END DO
       CLOSE(8)
       return

       END
