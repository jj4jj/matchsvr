
all: matching.pb.h matching.cex.hpp 

matching.cex.hpp: ./matching.proto ./pbdcex.core.hxx
	./pbdcexer -mMatching -pmatching.proto -I. -I/usr/local/include --cpp_out=.

matching.pb.h: ./matching.proto ./extensions.proto
	protoc $^ --cpp_out=. -I/usr/local/include -I.

clean:
	rm -f ./matching.cex.hpp
	rm -f ./matching.pb.h

test: ./matching.cpp ./test.cpp ./pbdcex.core.hpp
	g++ $^ -o $@ -g3 --std=c++11 -Wall

