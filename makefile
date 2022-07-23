FC = ifort 
FLAGS = -lfftw3

fourier_ps:
	$(FC) $(FLAGS) -c fourier_ps.f90 -o fourier_ps.o
	$(FC) $(FLAGS) fourier_ps.o -o fourier_ps

clean :
	rm -rf *.mod *~ *.o fourier_ps
full-clean : clean
	rm *.txt
