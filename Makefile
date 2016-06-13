test: ./test.cpp libmatcing.a
	g++ $^ -o $@ -g3 ${protoi} --std=c++11 -Wall


matching.cex.hpp: ./matching.proto ./pbdcex.core.hxx
	./pbdcexer -mMatching -pmatching.proto -I. ${protoi} --cpp_out=.

matching.pb.h: ./matching.proto ./extensions.proto
	${protoc} $^ --cpp_out=. ${protoi} -I.

clean:
	rm -f ./matching.cex.hpp
	rm -f ./matching.pb.*
	rm -f libmatcing.a
	rm -f test

libmatcing.a: matching.cpp matching.cex.hpp matching.pb.h
	g++ -c matching.cpp --std=c++11 -O2 ${protoi} -Wall
	ar -rcs $@ matching.o

install: test
