#include <sys/stat.h> // mkdir()
#include "kmerizer.h"

int main(int argc, char *argv[])
{
	// parse args
	if (argc != 4) {
		fprintf(stderr, "Usage: %s <k> <threads> <input dir>\n", argv[0]);
		return 1;
	}
	size_t k = atoi(argv[1]);
	size_t threads = atoi(argv[2]);
	char* inprefix = argv[3];

	kmerizer *counter = new kmerizer(k, threads, inprefix,'C');
	counter->load();
	counter->histogram();
	return 0;
}