kmerizer
========
Everyone's writing their own kmer counting software so here's one more. The approach is to do things very simply and hope that the result is efficient. (so far so good).

We use Heng Li's kseq.h to read sequences from fasta or fastq.
kmers are packed (2 bits per nucleotide) into 64 bit words. But before you do anything, allocate a huge chunk of memory to hold the packed kmers. As you read sequences, each one gets converted into a set of kmers. You may choose to canonicalize the kmer by choosing between the kmer and its reverse complement.

kmers are binned uniformly into 256 bins using a hash function with a random lookup table to distribute load evenly.

When the chunk of memory gets full, we take a break from reading input sequences, and sort and uniqify the kmers in each bin. qsort lets you do this by specifying a comparison function. We are using memcmp() because OS developers are probably better programmers. Once sorted, write the distinct kmers and their counts to output files. Then you reuse the huge chunk of memory and resume parsing kmers from the input sequences.

After all the sequences have been read in, another round of sort/uniq/serialization happens for the last batch of kmers. Then, if necessary, batches are merged.

Generating a histogram of kmer frequencies is done by reading the counts files from each bin. You generate a histogram for each bin and then merge them.


But how do you store the counts? bitmap index! A set of range encoded WAH bitvectors will let you quickly fetch kmers by frequency and generate histograms.