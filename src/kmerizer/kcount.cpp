#include <zlib.h>
#include <sys/stat.h> // mkdir()
#include "kseq.h"
#include "kmerizer.h"

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	// parse args
	if (argc != 7) {
		fprintf(stderr, "Usage: %s <input file> <k> <threads> <cap_bytes> <mode ('A|B|C')> <output dir>\n", argv[0]);
		return 1;
	}
	gzFile fp;
	fp = gzopen(argv[1], "r");
	size_t k = atoi(argv[2]);
	size_t threads = atoi(argv[3]);
	size_t cap_bytes = atoll(argv[4]);
	char* mode = argv[5];
	char* outprefix = argv[6];
	// create output directory if it doesn't exist
	mkdir(outprefix,0755);

	Kmerizer *counter = new Kmerizer(k, threads, outprefix, mode[0]);
	int rc = counter->allocate(cap_bytes);
	if (rc != 0) {
		fprintf(stderr,"failed to allocate %zi bytes\n",cap_bytes);
		exit(1);
	}
	// process each seq from input
	kseq_t *seq = kseq_init(fp);
	int length;
	while ((length = kseq_read(seq)) >= 0) {
        int offset=0;
        while (offset < length) {
            while (offset < length && seq->seq.s[offset] == 'N') offset++;
            // offset is next non-N
            int offset2=offset+1;
            while (offset2 < length && seq->seq.s[offset2] != 'N') offset2++;
            // offset2 is end of seq or next N
            if (offset2 - offset > k)
                counter->addSequence(seq->seq.s + offset, offset2 - offset);
            offset = offset2;
        }
    }
	kseq_destroy(seq);
	gzclose(fp);
	counter->save();
	return 0;
}