
test: test.cpp ../../eph_spline.cpp ../../eph_spline.h
	g++ -O2 -g -pg -std=c++11 -o test_profile test.cpp ../../eph_spline.cpp -I ../../
	g++ -O2 -g -std=c++11 -o test_O2 test.cpp ../../eph_spline.cpp -I ../../
	g++ -O3 -g -std=c++11 -o test_O3 test.cpp ../../eph_spline.cpp -I ../../
	g++ -Ofast -g -std=c++11 -o test_Ofast test.cpp ../../eph_spline.cpp -I ../../
	g++ -Ofast -g -DNDEBUG -std=c++11 -o test_ndebug test.cpp ../../eph_spline.cpp -I ../../
	g++ -Ofast -g -DNDEBUG -DEPH_UNSAFE -std=c++11 -o test_unsafe test.cpp ../../eph_spline.cpp -I ../../

clean:
	rm test_profile
	rm test_O2
	rm test_O3
	rm test_Ofast
	rm test_ndebug
	rm test_unsafe
