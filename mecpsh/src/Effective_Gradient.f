       SUBROUTINE Effective_Gradient(N,Ea,Eb,Ga,Gb,ParG,PerpG,G)
       implicit none
       
C      Computes the parallel and perpendicular compenents of the Effective Gradient,
C      As well as the effective gradient itself.
       
       integer n, i
       double precision Ea, Eb, Ga(N), Gb(N), G(N), ParG(N), PerpG(N), npg, pp
       double precision facPP, facP
       parameter (facPP=140.d0,facP=1.d0)
C  These factors are only really important for the first step
C  The "difference gradient" is typically ca. 0.075 Hartree/Bohr.
C      i.e. 0.14 Hartree/Angstrom.
C  Assuming this is constant, this gives a Hessian term for the func (Ea-Eb)**2
C     of ca. 0.01 Hartree**2 / Angstrom**2  (0.14**2 / 2)
C  The Hessian term along "normal" coordinates is empirically about 1.4 Hartree / Angstrom**2
C  Using facPP ~ 140 /Hartree means that along this coordinate, too, the Hessian is about right.
       
       npg = 0.d0
       pp = 0.d0
       do i = 1, n
           PerpG(i) = Ga(i) - Gb(i)
           npg = npg + PerpG(i)**2
           pp = pp + Ga(i) * PerpG(i)
       end do
       npg = sqrt(npg)
       pp = pp / npg
       do i = 1, n
           ParG(i) = Ga(i) - PerpG(i) / npg * pp
           G(i) = (Ea - Eb) * facPP * PerpG(i) + facP * ParG(i)
       end do
       return
       
       END
       
