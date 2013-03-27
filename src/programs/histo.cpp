#include <zlib.h>
#include <sys/stat.h> // mkdir()
#include "kseq.h"
#include "kmerizer.h"

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	// parse args
	if (argc != 7) {
		fprintf(stderr, "Usage: %s <input file> <k> <threads> <cap_bits> <mode ('A|B|C')> <output dir>\n", argv[0]);
		return 1;
	}
	gzFile fp;
	fp = gzopen(argv[1], "r");
	size_t k = atoi(argv[2]);
	size_t threads = atoi(argv[3]);
	size_t cap_bits = atoi(argv[4]);
	char* mode = argv[5];
	char* outprefix = argv[6];
	// create output directory if it doesn't exist
	mkdir(outprefix,0755);

	kmerizer *counter = new kmerizer(k, threads, outprefix, mode[0]);
	int rc = counter->allocate(1ULL<<cap_bits);
	if (rc != 0) {
		fprintf(stderr,"failed to allocate %llu bytes\n",1ULL<<cap_bits);
		exit(1);
	}
	// process each seq from input
	kseq_t *seq = kseq_init(fp);
	int length;
	while ((length = kseq_read(seq)) >= 0)
		if ((size_t)length >= k)
			counter->addSequence(seq->seq.s,length);
	kseq_destroy(seq);
	gzclose(fp);
	counter->histogram();
	return 0;
}