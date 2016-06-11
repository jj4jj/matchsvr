
all: matching.pb.h matching.cex.hpp 

matching.cex.hpp: ./matching.proto
	../pbdcex/pbdcexer -mMatching -pmatching.proto -I. -I/usr/local/include --cpp_out=.

matching.pb.h: ./matching.proto
	protoc ./matching.proto --cpp_out=. -I/usr/local/include -I.

