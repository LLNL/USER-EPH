
ROOT = ../../../
CFLAGS = -DHAVE_INLINE
LDFLAGS = -lgsl
INCFLAGS = -I $(ROOT)

test: test.cpp $(ROOT)/eph_spline.cpp $(ROOT)/eph_spline.h
	g++ -O2 -g -pg -std=c++11 -o test_profile test.cpp $(ROOT)/eph_spline.cpp $(INCFLAGS) $(LDFLAGS) $(CFLAGS)
	g++ -O2 -g -std=c++11 -o test_O2 test.cpp $(ROOT)/eph_spline.cpp $(INCFLAGS) $(LDFLAGS) $(CFLAGS)
	g++ -O3 -g -std=c++11 -o test_O3 test.cpp $(ROOT)/eph_spline.cpp $(INCFLAGS) $(LDFLAGS) $(CFLAGS)
	g++ -Ofast -g -std=c++11 -o test_Ofast test.cpp $(ROOT)/eph_spline.cpp $(INCFLAGS) $(LDFLAGS) $(CFLAGS)
	g++ -Ofast -g -DNDEBUG -std=c++11 -o test_ndebug test.cpp $(ROOT)/eph_spline.cpp $(INCFLAGS) $(LDFLAGS) $(CFLAGS)
	g++ -Ofast -g -DNDEBUG -DEPH_UNSAFE -std=c++11 -o test_unsafe test.cpp $(ROOT)/eph_spline.cpp $(INCFLAGS) $(LDFLAGS) $(CFLAGS)

clean:
	rm test_profile
	rm test_O2
	rm test_O3
	rm test_Ofast
	rm test_ndebug
	rm test_unsafe
