#include <zlib.h>
#include <sys/stat.h> // mkdir()
#include <unistd.h> // getopt()
#include "kseq.h"
#include "kmerizer.h"

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
    extern char *optarg;
    extern int optind;
    int ca=0,k=0,t=1,u=0,l=0;
    size_t s=0;
    char m,*program,*d,*f,*o;
    bool help=false;
    vector<char *> columns;
    m='A';
    program = NULL;
    while ((ca = getopt (argc, argv, "p:k:s:t:m:d:o:u:f:l:c:h?")) != -1)
        switch (ca) {
            case 'p': // program to run (count, stats, histo, filter, query)
                program = optarg;
                break;
            case 'k': // kmer length
                k = atoi(optarg);
                break;
            case 's': // size bytes to reserve for counting
                s = atoll(optarg);
                break;
            case 't': // number of threads
                t = atoi(optarg);
                break;
            case 'm': // counting mode (A|B|C)
                m = optarg[0];
                break;
            case 'd': // data directory (for reading)
                d = optarg;
                break;
            case 'o': // data directory (for writing)
                o = optarg;
                break;
            case 'f': // input file
                f = optarg;
                break;
            case 'l': // lower count
                l = atoi(optarg);
                break;
            case 'u': // upper count
                u = atoi(optarg);
                break;
            case 'c': // column name(s) for dump - defaults to seq and count
                columns.push_back(optarg);
                break;
            case 'h':
            case '?': // help required
                help = true;
                break;
        }
    if (! program) {
        fprintf(stderr,"Usage: %s <cmd> [options]\nWhere cmd is one of: count, stats, histo, filter, query\ntry %s <cmd> -?\n",argv[0],argv[0]);
        return 2;
    }
    if (help) {
        switch (program[0]) {
            case 'c': // count help
                fprintf(stderr,"count -f <fast[aq] input> -k <kmer length> -s <max memory> -t <number of threads> -m <mode (A|B|C)> -o <output dir>\n");
                break;
            case 's': // stats help
                fprintf(stderr,"stats -d <index dir> -k <kmer length>\n");
                break;
            case 'h': // histo help
                fprintf(stderr,"histo -d <index dir> -k <kmer length>\n");
                break;
            case 'f': // filter help
                fprintf(stderr,"filter -d <index dir> -k <kmer length> -l <lower count> -u <upper count> -c <columns> -t <number of threads>\n");
                break;
            case 'q': // query help
                fprintf(stderr,"query -d <index dir> -k <kmer length> -f <fast[aq] input> -c <columns>\n");
                break;
            default: // generic help
                fprintf(stderr,"Usage: %s <cmd> [options]\nWhere cmd is one of: count, stats, histo, filter, query\ntry %s <cmd> -?\n",argv[0],argv[0]);
                break;
        }
        return 1;
    }

    switch (program[0]) {
        case 'c': { // count
            Kmerizer *counter = new Kmerizer(k, t, o, m, 'C');
        	mkdir(o,0755);

        	int rc = counter->allocate(s);
        	if (rc != 0) {
        		fprintf(stderr,"failed to allocate %zi bytes\n",s);
        		exit(1);
        	}
        	// process each seq from input
            gzFile fp;
        	fp = gzopen(f, "r");
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
                    if (offset2 - offset >= k)
                        counter->addSequence(seq->seq.s + offset, offset2 - offset);
                    offset = offset2;
                }
            }
        	kseq_destroy(seq);
        	gzclose(fp);
        	counter->save();
            break;
        }
        case 's': { // stats
            Kmerizer *counter = new Kmerizer(k, t, d, m, 'Q');
            // for(int i=0;i<columns.size();i++)
            //     counter->stats(columns[i]);
            counter->stats("count");
            break;
        }
        case 'h': { // histo
            Kmerizer *counter = new Kmerizer(k, t, d, m, 'Q');
            // for(int i=0;i<columns.size();i++)
            //     counter->histo(columns[i]);
            counter->histo("count");
            break;
        }
        case 'f': { // filter
            Kmerizer *counter = new Kmerizer(k, t, d, m, 'Q');
            if (u>0)
                counter->filter("count",l,u); // TODO: arbitrary where clause
            counter->dump(columns);
            break;
        }
        case 'q': { // query
            Kmerizer *counter = new Kmerizer(k, t, d, m, 'Q');
        	int rc = counter->allocate(s);
        	if (rc != 0) {
        		fprintf(stderr,"failed to allocate %zi bytes\n",s);
        		exit(1);
        	}
        	// read seq from input and call query()
            // on each subsequence of [ACGT]*
            // multiple calls to query() OR the results (internal to kmerizer)
            // after the entire sequence has been processed, we can dump requested columns
            gzFile fp;
        	fp = gzopen(f, "r");
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
                    if (offset2 - offset >= k)
                        counter->addSequence(seq->seq.s + offset, offset2 - offset);
                    offset = offset2;
                }
            }
        	kseq_destroy(seq);
        	gzclose(fp);
            counter->save();
            counter->dump(columns);
            break;
        }
    }
    return 0;
}