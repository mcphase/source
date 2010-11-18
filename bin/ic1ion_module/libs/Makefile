all: libblas.a liblapack.a libarpack.a#libcblaswr.a

libblas.a: 
	cd blas ; make && cp libblas.a ../

libcblaswr.a:
	cd blaswrap ; make && cp libcblaswr.a ../

liblapack.a:
	cd lapack ; make && cp liblapack.a ../

libarpack.a:
	cd arpack ; make && cp libarpack.a ../
	#cd arpack-0.92 ; ./configure && make && cp SRC/.libs/libarpack.a ../

clean:
	cd blas ; make clean
	cd blaswrap ; make clean
	cd lapack ; make clean
	#cd arpack-0.92 ; ./configure && make clean
	cd arpack ; make clean
	rm -vf libblas.a libcblaswr.a liblapack.a libarpack.a
