FC	= gfortran
FLIB    = -lgfortran
FFLAGS	= -O3
GCC	= g++
JINC    = -I/Library/Java/Home/include
JLIB    = -L/sw/lib -L/sw/lib/gcc4/lib




#############

OBJS   = dgexpv-beast.o blas.o lapack.o dgpadm.o SparseMatrixExponential.o

jni: $(OBJS) 
	 $(GCC) -o libExpoKit.jnilib -framework JavaVM -dynamiclib \
		$(OBJS) $(FLIB) $(JLIB)

.f.o:;  $(FC) $(FFLAGS) -c $<

.cpp.o:; $(GCC) $(FFLAGS) $(JINC) -c $<

clean:; rm *.o
	rm *.jnilib



