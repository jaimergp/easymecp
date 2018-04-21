      PROGRAM Optimizer
      implicit none

C     New version, Nov. 2003. All gradients converted to Hartree/Ang.
C     facPp chosen so that inverse Hessian ~ diagonal 0.7 Ang*Ang/Hartree

      INTEGER Natom, Nx
      parameter (Natom=77)
C         Set Natom to the number of atoms.
      parameter (Nx=3*Natom)

      INTEGER Nstep, FFile, Conv, AtNum(Natom)
      DOUBLE PRECISION Ea_1, Ea_2, Eb_1, Eb_2, IH_1(Nx,Nx), IH_2(Nx,Nx), ParG(Nx), PerpG(Nx)
      DOUBLE PRECISION Ga_1(Nx), Ga_2(Nx), Gb_1(Nx), Gb_2(Nx), X_1(Nx), X_2(Nx), X_3(Nx)
      DOUBLE PRECISION G_1(Nx), G_2(Nx)

      CALL ReadProgFile(Natom,Nx,AtNum,Nstep,FFile,X_1,X_2,IH_1,Ea_1,Eb_1,Ga_1,Gb_1,G_1)
C          Recover data from the previous Steps
      CALL ReadInput(Natom,Nx,Ea_2,Ga_2,Eb_2,Gb_2)
C          Read in the new ab initio Energies and Gradients
      CALL Effective_Gradient(Nx,Ea_2,Eb_2,Ga_2,Gb_2,ParG,PerpG,G_2)
C          Compute the Effective Gradient fac * DeltaE * PerpG + ParG
      CALL UpdateX(Nx,Nstep,FFile,X_1,X_2,X_3,G_1,G_2,IH_1,IH_2)
C          the BFGS step
      CALL TestConvergence(Nx,Natom,Nstep,AtNum,Ea_2,Eb_2,X_2,X_3,ParG,PerpG,G_2,Conv)
C          Checks Delta E, X, Mag. Delta G. Writes Output File
      IF (Conv .ne. 1) THEN
           CALL WriteGeomFile(Natom,Nx,AtNum,X_3)
           CALL WriteProgFile(Natom,Nx,AtNum,Nstep,X_2,X_3,IH_2,Ea_2,Eb_2,Ga_2,Gb_2,G_2)
      END IF

      END
