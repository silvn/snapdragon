CC = g++

all: test_bvec_density

test_bvec_density: test_bvec_density.cpp bvec32.cpp
	$(CC) -O test_bvec_density.cpp bvec32.cpp -o test_bvec_density

test: test_bvec_density
	@ ./test_bvec_density 1000000
	
clean:
	rm test_bvec_density

.PHONY: clean