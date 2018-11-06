      SUBROUTINE TestConvergence(N,Natom,Nstep,AtNum,Ea,Eb,X_2,X_3,ParG,PerpG,G,Conv)
      implicit none

C     Checks convergence, and updates report file
C     There are 5 criteria for testing convergence
C     They are the same as in Gaussian (Except the last):
C     Av.DeltaX, Max.DeltaX, Av.Grad., Max.Grad., DeltaE

      INTEGER N, Natom, AtNum(Natom), Nstep, Conv, i, j, k
      DOUBLE PRECISION Ea, Eb, X_2(N), ParG(N), PerpG(N), X_3(N), G(N)

      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(N), DE, DXMax, DXRMS, GMax, GRMS, PpGRMS, PGRMS, TDE
      DOUBLE PRECISION TDXMax, TDxRMS, TGMax, TGRMS
      PARAMETER (TDE=5.d-5,TDXMax=4.d-3,TDXRMS=2.5d-3,TGMax=7.d-4,TGRMS=5.d-4)

      OPEN(UNIT=8,FILE="AddtoReportFile")
      IF (Nstep .eq. 0) THEN
          write (8,*) "      Geometry Optimization of an MECP"
          write (8,*) "      Program: J. N. Harvey, March 1999"
          write (8,*) "        version 2, November 2003"
          write (8,*)
          write (8,'(A)') "Initial Geometry:"
          DO I = 1, Natom
              k = 3 * (I-1) + 1
              write (8,'(I3,3F15.7)') AtNum(I), (X_2(j),j=k,k+2)
          END DO
          write (8,*)
      END IF
      DE = ABS(Ea - Eb)
      DXMax = 0.d0
      DXRMS = 0.d0
      GMax = 0.d0
      GRMS = 0.d0
      PpGRMS = 0.d0
      PGRMS = 0.d0
      DO I = 1, n
          DeltaX(i) = X_3(i) - X_2(i)
          IF (ABS(DeltaX(i)) .gt. DXMax) DXMax = ABS(DeltaX(i))
          DXRMS = DXRMS + DeltaX(i)**2
          IF (ABS(G(i)) .gt. Gmax) Gmax = ABS(G(i))
          GRMS = GRMS + G(i)**2
          PpGRMS = PpGRMS + PerpG(i)**2
          PGRMS = PGRMS + ParG(i)**2
      END DO
      DXRMS = SQRT(DXRMS / N)
      GRMS = SQRT(GRMS / N)
      PpGRMS= SQRT(PpGRMS / N)
      PGRMS = SQRT(PGRMS / N)

      Conv = 0
      do i = 1, 5
          flags(i) = " NO"
          PConv(i) = .false.
      end do

      IF (GMax .lt. TGMax) THEN
          PConv(1) = .true.
          flags(1) = "YES"
      END IF
      IF (GRMS .lt. TGRMS) THEN
          PConv(2) = .true.
          flags(2) = "YES"
      END IF
      IF (DXMax .lt. TDXMax) THEN
          PConv(3) = .true.
          flags(3) = "YES"
      END IF
      IF (DXRMS .lt. TDXRMS) THEN
          PConv(4) = .true.
          flags(4) = "YES"
      END IF
      IF (DE .lt. TDE) THEN
          PConv(5) = .true.
          flags(5) = "YES"
      END IF
      IF (PConv(1) .and. PConv(2) .and. PConv(3) .and. PConv(4) .and. PConv(5)) THEN
          Conv = 1
      ELSE
          Nstep = Nstep + 1
      END IF

      write (8,'(A,F18.10)') "Energy of First State:  ",Ea
      write (8,'(A,F18.10)') "Energy of Second State: ",Eb
      write (8,*)
      write (8,'(A)') "Convergence Check (Actual Value, then Threshold, then Status):"
      write (8,'(A,F11.6,A,F8.6,A,A)') "Max Gradient El.:", GMax," (",TGMax,")  ",flags(1)
      write (8,'(A,F11.6,A,F8.6,A,A)') "RMS Gradient El.:", GRMS," (",TGRMS,")  ",flags(2)
      write (8,'(A,F11.6,A,F8.6,A,A)') "Max Change of X: ",DXMax," (",TDXMax,")  ",flags(3)
      write (8,'(A,F11.6,A,F8.6,A,A)') "RMS Change of X: ",DXRMS," (",TDXRMS,")  ",flags(4)
      write (8,'(A,F11.6,A,F8.6,A,A)') "Difference in E: ",DE," (",TDE,")  ",flags(5)
      write (8,*)
      write (8,'(A)') "Overall Effective Gradient:"
      DO I = 1, Natom
          k = 3 * (I - 1) + 1
          write (8,'(I3,3F16.8)') I, (G(J),J=k,k+2)
      END DO
      write (8,*)
      write (8,'(A,F11.6,A)') "Difference Gradient: (RMS * DE:",PpGRMS,")"
      DO I = 1, Natom
          k = 3 * (I - 1) + 1
          write (8,'(I3,3F16.8)') I, (PerpG(J),J=k,k+2)
      END DO
      write (8,*)
      write (8,'(A,F11.6,A)') "Parallel Gradient: (RMS:",PGRMS,")"
      DO I = 1, Natom
          k = 3 * (I - 1) + 1
          write (8,'(I3,3F16.8)') I, (ParG(J),J=k,k+2)
      END DO
      write (8,*)

      IF (Conv .eq. 1) THEN
          write (8,'(A)') "The MECP Optimization has CONVERGED at that geometry !!!"
          write (8,'(A)') "Goodbye and fly with us again..."
      ELSE
          write (8,'(A,I3)') "Geometry at Step",Nstep
          DO I = 1, Natom
              k = 3 * (I - 1) + 1
              write (8,'(I3,3F15.7)') AtNum(I), (X_3(J),J=k,k+2)
          END DO
          write (8,*)
      END IF

      CLOSE(8)
      return

      END

