
test: test.cpp ../../../eph_fdm.cpp ../../../eph_fdm.h ../../../eph_spline.cpp ../../../eph_spline.h 
	mpic++ -O2 -g -DEPH_OMP -DNDEBUG -fopenmp -std=c++11 -o test test.cpp ../../../eph_fdm.cpp ../../../eph_spline.cpp -I ../../../

clean:
	rm test
