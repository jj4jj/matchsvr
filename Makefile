PBDC_INCLUDE=../pbdcex/include

test/test: test/test.cpp lib/libmatcing.a
	g++ $^ -o $@ -g3 ${protoi} -I${PBDC_INCLUDE} -Iproto -Iinclude --std=c++11 -Wall
	mkdir -p bin
	mv test/test bin/

proto/matching.cex.hpp: proto/matching.proto
	./tools/pbdcexer -mMatching -pmatching.proto -Iproto -I${PBDC_INCLUDE} ${protoi} --cpp_out=proto

proto/matching.pb.h: proto/matching.proto ${PBDC_INCLUDE}/extensions.proto
	${protoc} $^ --cpp_out=proto ${protoi} -Iproto -I${PBDC_INCLUDE}

clean:
	rm -rf lib/ bin/
	rm -f proto/matching.cex.hpp
	rm -f proto/matching.pb.*

libmatcing.a: src/matching.cpp proto/matching.cex.hpp proto/matching.pb.h
	g++ -c src/matching.cpp --std=c++11 -O2 ${protoi} -Iproto -Wall -I${PBDC_INCLUDE}
	ar -rcs $@ matching.o
test: test/test
install: libmatcing.a test
	mkdir -p lib
	mv libmatcing.a lib/
	mkdir -p include
	cp -f src/matching.h include/
