

.PHONY: tests
tests: all
	./test

all: test.cpp ../../eph_spline.h ../../eph_kappa.h
	g++ -O2 -g -std=c++11 -o test test.cpp -I ../../

clean:
	rm test
