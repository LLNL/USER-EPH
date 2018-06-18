
test: test.cpp ../../../eph_fdm.cpp ../../../eph_fdm.h 
	mpic++ -O2 -g -DEPH_TESTING -DNDEBUG -fopenmp -std=c++11 -o test test.cpp ../../../eph_fdm.cpp -I ../../../

clean:
	rm test
