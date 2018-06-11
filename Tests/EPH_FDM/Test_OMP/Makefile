
all: test test_omp

test: test.cpp ../../../eph_fdm.cpp ../../../eph_fdm.h 
	mpic++ -O2 -g -std=c++11 -o test test.cpp ../../../eph_fdm.cpp -I ../../../

test_omp: test.cpp ../../../eph_fdm.cpp ../../../eph_fdm.h 
	mpic++ -DEPH_TESTING -O2 -g -fopenmp -std=c++11 -o test_omp test.cpp ../../../eph_fdm.cpp -I ../../../
clean:
	rm test
	rm test_omp
