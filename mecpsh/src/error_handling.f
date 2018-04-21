      SUBROUTINE error_handling(e)
      implicit none

      integer e

      IF (e .eq. -1) THEN
          write (*,*) "The ProgFile does not exist for some reason."
          write (*,*) "Please create it then try again."
      END IF

      IF (e .eq. -2) THEN
          write (*,*) "Incomplete ProgFile - did you set Nstep right ?"
      END IF

      IF (e .eq. -3) THEN
          write (*,*) "Problem Reading Gaussian Output"
          OPEN(UNIT=302,FILE="AddtoReportFile")
          write (302,*) "ERROR IN GAUSSIAN STEP"
          CLOSE(302)
      END IF

      IF (e .eq. -4) THEN
           write (*,*) "Wrong Number of Atoms - Recompile !!"
      END IF

      IF (e .eq. -5) THEN
           write (*,*) "Problem with the ab initio Job."
      END IF

      STOP

      END

