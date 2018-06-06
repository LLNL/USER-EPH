
test: test.cpp ../../eph_spline.cpp ../../eph_spline.h
	g++ -O2 -g -std=c++11 -o test test.cpp ../../eph_spline.cpp -I ../../
	g++ -O2 -g -DNDEBUG -DEPH_UNSAFE -std=c++11 -o test_unsafe test.cpp ../../eph_spline.cpp -I ../../

clean:
	rm test
	rm test_unsafe
