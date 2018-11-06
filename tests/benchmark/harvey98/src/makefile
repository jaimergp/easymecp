######  Fortran Complier  ##################################
#
FC = gfortran
#
######  Generic fortran compilation flag  ##################
#
FFLAGS= -O -ffixed-line-length-none
#
############################################################
#
OBJS = \
Effective_Gradient.o  MainProgram.o         ReadProgFile.o        WriteGeomFile.o \
ReadInput.o           TestConvergence.o     WriteProgFile.o \
Initialize.o          UpdateX.o             error_handling.o
#
############################################################
#
MECP.x  : $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o MECP.x
#
############################################################
clean	:
	rm -f *.o MECP.x




