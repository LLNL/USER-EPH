
EXTRA = -g -std=c++11 -march=native -mtune=native

test: test.cpp ../../eph_spline.cpp ../../eph_spline.h
	g++ -O2 -pg $(EXTRA) -o test_profile test.cpp ../../eph_spline.cpp -I ../../
	g++ -O2 $(EXTRA) -o test_O2 test.cpp ../../eph_spline.cpp -I ../../
	g++ -O3 $(EXTRA) -o test_O3 test.cpp ../../eph_spline.cpp -I ../../
	g++ -Ofast $(EXTRA) -o test_Ofast test.cpp ../../eph_spline.cpp -I ../../
	g++ -Ofast -DNDEBUG $(EXTRA) -o test_ndebug test.cpp ../../eph_spline.cpp -I ../../
	g++ -Ofast -DNDEBUG -DEPH_UNSAFE $(EXTRA) -o test_unsafe test.cpp ../../eph_spline.cpp -I ../../

clean:
	rm test_profile
	rm test_O2
	rm test_O3
	rm test_Ofast
	rm test_ndebug
	rm test_unsafe
