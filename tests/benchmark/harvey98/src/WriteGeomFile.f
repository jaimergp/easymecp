      SUBROUTINE WriteGeomFile(Natom,Nx,AtNum,X)
      implicit none

      INTEGER i, j, Natom, Nx, AtNum(Natom)
      DOUBLE PRECISION X(Nx)

      OPEN (UNIT=8,FILE="geom")
      DO I = 1, Natom
          write (8,'(I4,3F14.8)') AtNum(I), (X(J),J=3*(I-1)+1,3*(I-1)+3)
      END DO
      write (8,*)
      CLOSE(8)
      return

      END

