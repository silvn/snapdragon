Kmerizer
========
Kmerizer is a library that implements functions for counting k-mers for k<=256. It counts k-mers in a naive but fast way, and then stores the distinct kmers and their counts using compressed bitvectors. Without this compression, the largest component of the run time is the IO required to write distinct k-mers to files. Therefore, reducing the size of the output files will have a huge impact on the overall performance of the software.

Overview:
K-mers from an input sequence are packed (2 bits per nucleotide) into 64 bit words. Each k-mer is optionally canonicalized by comparing it to its reverse complement and selecting the minimum. Memory is allocated in advance for holding the bit packed k-mers. For parallelization, the k-mers are uniformly hashed into one of 256 bins to distribute load evenly. If one of the bins fills before the input sequences have been processed, the k-mers in the 256 buffers are counted and written to output files. Then, memory is reused for loading more k-mers. After all the sequences have been read in, another round of serialization happens for the last batch of k-mers. Then, if necessary, batches are merged.

Counting:
K-mers are counted by sorting them in place and then iterating over them to identify and count the distinct k-mers. When k>32 multiple words are required to hold the packed sequence, and a custom multi-word comparison function is used.

Compression:
The sorted distinct k-mers each occupy 64*ceil(k/32) bits. We reduce the overall run time significantly by spending some CPU cycles to create a bit-sliced bitmap index of compressed bitvectors (one per bit position.) The high-order bits compress extremely well, while low order bits only require a small amount of overhead. Similarly, the counts are converted to a range encoded bitmap index with one compressed bitvector for each k-mer frequency. In this index, bitvector f marks the distinct k-mers that occur <= f times. The bitvector caches the number of set bits, so it is easy to calculate a histogram of k-mer frequencies.
