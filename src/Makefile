default: test-exhaust bench

DEFINES := -DHAVE_CLZ -O2 -g0

fields/small.o: fields/small.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/small.cpp -c -o fields/small.o
fields/generic_1byte.o: fields/generic_common_impl.h fields/generic_1byte.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_1byte.cpp -c -o fields/generic_1byte.o
fields/generic_2bytes.o: fields/generic_common_impl.h fields/generic_2bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_2bytes.cpp -c -o fields/generic_2bytes.o
fields/generic_3bytes.o: fields/generic_common_impl.h fields/generic_3bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_3bytes.cpp -c -o fields/generic_3bytes.o
fields/generic_4bytes.o: fields/generic_common_impl.h fields/generic_4bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_4bytes.cpp -c -o fields/generic_4bytes.o
fields/generic_5bytes.o: fields/generic_common_impl.h fields/generic_5bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_5bytes.cpp -c -o fields/generic_5bytes.o
fields/generic_6bytes.o: fields/generic_common_impl.h fields/generic_6bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_6bytes.cpp -c -o fields/generic_6bytes.o
fields/generic_7bytes.o: fields/generic_common_impl.h fields/generic_7bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_7bytes.cpp -c -o fields/generic_7bytes.o
fields/generic_8bytes.o: fields/generic_common_impl.h fields/generic_8bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 fields/generic_8bytes.cpp -c -o fields/generic_8bytes.o
fields/clmul_1byte.o: fields/clmul_common_impl.h fields/clmul_1byte.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_1byte.cpp -c -o fields/clmul_1byte.o
fields/clmul_2bytes.o: fields/clmul_common_impl.h fields/clmul_2bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_2bytes.cpp -c -o fields/clmul_2bytes.o
fields/clmul_3bytes.o: fields/clmul_common_impl.h fields/clmul_3bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_3bytes.cpp -c -o fields/clmul_3bytes.o
fields/clmul_4bytes.o: fields/clmul_common_impl.h fields/clmul_4bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_4bytes.cpp -c -o fields/clmul_4bytes.o
fields/clmul_5bytes.o: fields/clmul_common_impl.h fields/clmul_5bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_5bytes.cpp -c -o fields/clmul_5bytes.o
fields/clmul_6bytes.o: fields/clmul_common_impl.h fields/clmul_6bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_6bytes.cpp -c -o fields/clmul_6bytes.o
fields/clmul_7bytes.o: fields/clmul_common_impl.h fields/clmul_7bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_7bytes.cpp -c -o fields/clmul_7bytes.o
fields/clmul_8bytes.o: fields/clmul_common_impl.h fields/clmul_8bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DNDEBUG -mpclmul -Wall -std=c++11 fields/clmul_8bytes.cpp -c -o fields/clmul_8bytes.o
minisketch.o: minisketch.cpp ../include/minisketch.h
	$(CXX) $(DEFINES) -DNDEBUG -Wall -std=c++11 minisketch.cpp -c -o minisketch.o

fields/small_test.o: fields/small.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/small.cpp -c -o fields/small_test.o
fields/generic_1byte_test.o: fields/generic_common_impl.h fields/generic_1byte.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_1byte.cpp -c -o fields/generic_1byte_test.o
fields/generic_2bytes_test.o: fields/generic_common_impl.h fields/generic_2bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_2bytes.cpp -c -o fields/generic_2bytes_test.o
fields/generic_3bytes_test.o: fields/generic_common_impl.h fields/generic_3bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_3bytes.cpp -c -o fields/generic_3bytes_test.o
fields/generic_4bytes_test.o: fields/generic_common_impl.h fields/generic_4bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_4bytes.cpp -c -o fields/generic_4bytes_test.o
fields/generic_5bytes_test.o: fields/generic_common_impl.h fields/generic_5bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_5bytes.cpp -c -o fields/generic_5bytes_test.o
fields/generic_6bytes_test.o: fields/generic_common_impl.h fields/generic_6bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_6bytes.cpp -c -o fields/generic_6bytes_test.o
fields/generic_7bytes_test.o: fields/generic_common_impl.h fields/generic_7bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_7bytes.cpp -c -o fields/generic_7bytes_test.o
fields/generic_8bytes_test.o: fields/generic_common_impl.h fields/generic_8bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 fields/generic_8bytes.cpp -c -o fields/generic_8bytes_test.o
fields/clmul_1byte_test.o: fields/clmul_common_impl.h fields/clmul_1byte.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_1byte.cpp -c -o fields/clmul_1byte_test.o
fields/clmul_2bytes_test.o: fields/clmul_common_impl.h fields/clmul_2bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_2bytes.cpp -c -o fields/clmul_2bytes_test.o
fields/clmul_3bytes_test.o: fields/clmul_common_impl.h fields/clmul_3bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_3bytes.cpp -c -o fields/clmul_3bytes_test.o
fields/clmul_4bytes_test.o: fields/clmul_common_impl.h fields/clmul_4bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_4bytes.cpp -c -o fields/clmul_4bytes_test.o
fields/clmul_5bytes_test.o: fields/clmul_common_impl.h fields/clmul_5bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_5bytes.cpp -c -o fields/clmul_5bytes_test.o
fields/clmul_6bytes_test.o: fields/clmul_common_impl.h fields/clmul_6bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_6bytes.cpp -c -o fields/clmul_6bytes_test.o
fields/clmul_7bytes_test.o: fields/clmul_common_impl.h fields/clmul_7bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_7bytes.cpp -c -o fields/clmul_7bytes_test.o
fields/clmul_8bytes_test.o: fields/clmul_common_impl.h fields/clmul_8bytes.cpp lintrans.h sketch.h sketch_impl.h int_utils.h
	$(CXX) $(DEFINES) -mpclmul -DMINISKETCH_VERIFY -Wall -std=c++11 fields/clmul_8bytes.cpp -c -o fields/clmul_8bytes_test.o
minisketch_test.o: minisketch.cpp ../include/minisketch.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 minisketch.cpp -c -o minisketch_test.o

libminisketch.a: fields/small.o fields/generic_1byte.o fields/generic_2bytes.o fields/generic_3bytes.o fields/generic_4bytes.o fields/generic_5bytes.o fields/generic_6bytes.o fields/generic_7bytes.o fields/generic_8bytes.o fields/clmul_1byte.o fields/clmul_2bytes.o fields/clmul_3bytes.o fields/clmul_4bytes.o fields/clmul_5bytes.o fields/clmul_6bytes.o fields/clmul_7bytes.o fields/clmul_8bytes.o minisketch.o
	$(AR) cr libminisketch.a fields/small.o fields/generic_1byte.o fields/generic_2bytes.o fields/generic_3bytes.o fields/generic_4bytes.o fields/generic_5bytes.o fields/generic_6bytes.o fields/generic_7bytes.o fields/generic_8bytes.o fields/clmul_1byte.o fields/clmul_2bytes.o fields/clmul_3bytes.o fields/clmul_4bytes.o fields/clmul_5bytes.o fields/clmul_6bytes.o fields/clmul_7bytes.o fields/clmul_8bytes.o minisketch.o

libminisketch_test.a: fields/small_test.o fields/generic_1byte_test.o fields/generic_2bytes_test.o fields/generic_3bytes_test.o fields/generic_4bytes_test.o fields/generic_5bytes_test.o fields/generic_6bytes_test.o fields/generic_7bytes_test.o fields/generic_8bytes_test.o fields/clmul_1byte_test.o fields/clmul_2bytes_test.o fields/clmul_3bytes_test.o fields/clmul_4bytes_test.o fields/clmul_5bytes_test.o fields/clmul_6bytes_test.o fields/clmul_7bytes_test.o fields/clmul_8bytes_test.o minisketch_test.o
	$(AR) cr libminisketch_test.a fields/small_test.o fields/generic_1byte_test.o fields/generic_2bytes_test.o fields/generic_3bytes_test.o fields/generic_4bytes_test.o fields/generic_5bytes_test.o fields/generic_6bytes_test.o fields/generic_7bytes_test.o fields/generic_8bytes_test.o fields/clmul_1byte_test.o fields/clmul_2bytes_test.o fields/clmul_3bytes_test.o fields/clmul_4bytes_test.o fields/clmul_5bytes_test.o fields/clmul_6bytes_test.o fields/clmul_7bytes_test.o fields/clmul_8bytes_test.o minisketch_test.o

test-exhaust.o: test-exhaust.cpp ../include/minisketch.h
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 test-exhaust.cpp -pthread -c -o test-exhaust.o

bench.o: bench.cpp ../include/minisketch.h
	$(CXX) $(DEFINES) -Wall -std=c++11 bench.cpp -c -o bench.o

test-exhaust: libminisketch_test.a test-exhaust.o
	$(CXX) $(DEFINES) -DMINISKETCH_VERIFY -Wall -std=c++11 test-exhaust.o libminisketch_test.a -pthread -o test-exhaust

bench: libminisketch.a bench.o
	$(CXX) $(DEFINES) -Wall -std=c++11 bench.o libminisketch.a -o bench

.PHONY: clean
clean:
	-rm fields/*.o *.o libminisketch.a test-exhaust bench
