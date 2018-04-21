      SUBROUTINE WriteProgFile(natom,nx,AtNum,Nstep,X_2,X_3,HI,Ea,Eb,Ga,Gb,G)
      implicit none

      INTEGER natom, nx, AtNum(Natom), Nstep
      DOUBLE PRECISION X_2(Nx), X_3(Nx), HI(Nx,Nx), Ea, Eb, Ga(Nx), Gb(Nx), G(Nx)

      INTEGER I, J

      OPEN(UNIT=8,FILE="ProgFile")
      WRITE(8,*) "Progress File for MECP Optimization"
      WRITE(8,*) "Number of Atoms:"
      WRITE(8,*) Natom
      WRITE(8,*) "Number of Steps already Run"
      WRITE(8,*) Nstep
      WRITE(8,*) "Is this a full ProgFile ?"
      WRITE(8,*) 1

      WRITE(8,*) "Next Geometry to Compute:"
      DO I = 1, natom
          WRITE(8,'(I3,3F20.12)') AtNum(I), (X_3(J),J=3*(I-1)+1,3*(I-1)+3)
      END DO

      WRITE(8,*) "Previous Geometry:"
      DO I = 1, Natom
          write (8,'(3F20.12)') (X_2(J),J=3*(I-1)+1,3*(I-1)+3)
      END DO
      WRITE(8,*) "Energies of First, Second State at that Geometry:"
      WRITE(8,'(F20.12)') Ea
      WRITE(8,'(F20.12)') Eb
      WRITE(8,*) "Gradient of First State at that Geometry:"
      DO I = 1, nx
          WRITE(8,'(F20.12)') Ga(I)
      END DO
      WRITE(8,*) "Gradient of Second State at that Geometry:"
      DO I = 1, nx
          WRITE(8,'(F20.12)') Gb(I)
      END DO
      WRITE(8,*) "Effective Gradient at that Geometry:"
      DO I = 1, nx
          WRITE(8,'(F20.12)') G(I)
      END DO
      WRITE(8,*) "Approximate Inverse Hessian at that Geometry:"
      DO I = 1, nx
          DO J = 1, nx
              WRITE(8,'(F20.12)') HI(I,J)
          END DO
      END DO

      CLOSE(8)
      return

      END
