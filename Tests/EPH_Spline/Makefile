

.PHONY: tests
tests: all
	./test

all: test.cpp eph_spline.h 
	g++ -O2 -g -std=c++11 -o test test.cpp

clean:
	rm test
