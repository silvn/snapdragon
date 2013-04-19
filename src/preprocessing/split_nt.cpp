#include "bvec/kseq.h"
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <sstream>
#include <fstream>
#include <errno.h>
#include <sys/stat.h>

using namespace std;

// void write_batch(vector<size_t> *offsets, vector<kseq_t*> &batch) {
// 	for(vector<size_t>::iterator i = offsets->begin(); i != offsets->end(); ++i) {
// 		
// 	}
// }

inline uint32_t vsearch(vector<uint32_t> &a, vector<uint32_t> &b, uint32_t v) {
	vector<uint32_t>::iterator it = lower_bound(a.begin(), a.end(), v);
	if (*it == v)
		return b[uint32_t(it - a.begin())];
	return 0;
}

template<typename T>
inline string toString(T t) {
	stringstream s;
	s << t;
	return s.str();
}

void write_fasta(stringstream *ssp,uint32_t tax, vector<uint32_t> &nodes_tax, vector<uint32_t> &nodes_parent) {
	vector<uint32_t> path;
	path.push_back(tax);
	while (tax != 1) {
		tax = vsearch(nodes_tax,nodes_parent,tax);
		if (tax == 0) {
			fprintf(stderr,"taxonomy fail parent(%u) undefined\n",path.back());
			exit(1);
		}
		path.push_back(tax);
	}
	stringstream path_ss;
	path_ss << "./";
	for(int i=path.size()-1;i>=0;i--)
		path_ss << path[i] << "/";
	path_ss << "split_nt.fa";
	ofstream outfile;
	outfile.open(path_ss.str().c_str(), ios::out | ios::app);
	outfile.write(ssp->str().c_str(),ssp->tellp());
	outfile.close();
}

void df_traverse(vector<uint32_t> &nodes_tax, vector<vector<uint32_t>*> &nodes_children, uint32_t parent, string path_str) {
	// append parent to the path_string
	// and mkdir
	path_str += "/" + toString(parent);
	mkdir(path_str.c_str(),0755);
	// lookup child nodes of parent
	vector<uint32_t>::iterator it = lower_bound(nodes_tax.begin(), nodes_tax.end(), parent);
	if (*it == parent) {
		uint32_t offset = uint32_t(it - nodes_tax.begin());
		if (nodes_children[offset]->size() > 0)
			for(vector<uint32_t>::iterator it2 = nodes_children[offset]->begin(); it2 < nodes_children[offset]->end(); it2++)
				df_traverse(nodes_tax, nodes_children, *it2, path_str);
	}
	return;
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	if (argc != 4) {
		fprintf(stderr, "Usage: %s gi_taxid_nucl.dmp nodes.dmp nt.gz\n", argv[0]);
		return 1;
	}

	// read gi_taxid_nucl.dmp.gz into memory as two vectors
	clock_t start,end;
	start = clock();
	gzFile file;
	file = gzopen(argv[1], "rb");
	if (file == NULL) {
		fprintf(stderr, "gzopen of '%s' failed: %s.\n", argv[1], strerror(errno));
		exit(1);
	}
	int uncomprLen = 50; // how do you choose this?
	char uncompr[50];
	vector<uint32_t> gi_vec,tax_vec;
	while (! gzeof(file)) {
		gzgets(file, uncompr, uncomprLen);
		uint32_t gi,tax;
		sscanf(uncompr,"%u\t%u",&gi,&tax);
		gi_vec.push_back(gi);
		tax_vec.push_back(tax);
	}
	gzclose(file);

	end = clock();
	fprintf(stderr,"parsing gi_taxid_nucl.dmp into two vectors of size %zi took %f seconds\n",
		gi_vec.size(),double(end-start)/CLOCKS_PER_SEC);
	
	start = clock();
	for(int i=0;i<1000000;i++) {
		uint32_t j = rand() % gi_vec.size();
		uint32_t t = vsearch(gi_vec,tax_vec,gi_vec[j]);
		if (t!=tax_vec[j]) {
			fprintf(stderr,"unexpected value tax_id(%u) is %u, but found %u\n",gi_vec[j],tax_vec[j],t);
			return 1;
		}
	}
	end = clock();
	fprintf(stderr,"1000000 gi->tax lookups took %f seconds\n",
		double(end-start)/CLOCKS_PER_SEC);

	// read nodes.dmp into memory as two vectors
	start = clock();
	FILE *pFile;
	pFile = fopen(argv[2],"r");
	vector<uint32_t> nodes_tax;
	vector<uint32_t> nodes_parent;
	char etc[1024];
	while (!feof(pFile)) {
		uint32_t tax,parent;
		fscanf(pFile,"%u\t|\t%u%[^\n]",&tax,&parent,etc);
		nodes_tax.push_back(tax);
		nodes_parent.push_back(parent);
	}
	fclose(pFile);
	end = clock();
	fprintf(stderr,"parsing nodes.dmp into two vectors of size %zi took %f seconds\n",
		nodes_tax.size(),double(end-start)/CLOCKS_PER_SEC);

	start = clock();
	for(int i=0;i<1000000;i++) {
		uint32_t j = rand() % nodes_tax.size();
		uint32_t t = vsearch(nodes_tax,nodes_parent,nodes_tax[j]);
		if (t!=nodes_parent[j]) {
			fprintf(stderr,"unexpected value parent(%u) is %u, but found %u\n",nodes_tax[j],nodes_parent[j],t);
			return 1;
		}
	}
	end = clock();
	fprintf(stderr,"1000000 tax->parent lookups took %f seconds\n",
		double(end-start)/CLOCKS_PER_SEC);

	// generate a list of child nodes for each tax id
	start = clock();
	vector<vector<uint32_t>*> nodes_children;
	for(uint32_t i=0; i<nodes_tax.size(); ++i) {
		vector<uint32_t> *v = new vector<uint32_t>;
		nodes_children.push_back(v);
	}
	for(uint32_t i=0; i<nodes_tax.size(); ++i) {
		uint32_t child = nodes_tax[i];
		uint32_t parent = nodes_parent[i];
		if (child != parent) {
			vector<uint32_t>::iterator it = lower_bound(nodes_tax.begin(), nodes_tax.end(), parent);
			if (*it == parent)
				nodes_children[uint32_t(it - nodes_tax.begin())]->push_back(child);
		}
	}
	end = clock();
	fprintf(stderr,"generating tax->children vector took %f seconds\n",
		double(end-start)/CLOCKS_PER_SEC);

	// do a depth-first traversal of the taxonomy tree
	// and create the directories on the file system
	start = clock();
	df_traverse(nodes_tax,nodes_children,1,"./");
	end = clock();
	fprintf(stderr,"generating taxonomy tree on disk took %f seconds\n",
		double(end-start)/CLOCKS_PER_SEC);


	// split the nt.gz according to taxonomy
	start = clock();
	kseq_t *seq;
	int length;
	gzFile fp = gzopen(argv[3], "r");
	seq = kseq_init(fp);
	uint32_t n_seqs=0;
	// fill a buffer for each tax_id.  when a buffer gets too full, write it to disk
	vector<stringstream*> fasta;
	for(uint32_t i=0; i<nodes_tax.size(); ++i)
		fasta.push_back(new stringstream);

	int min_seq_length = 32;
	int max_buffer = 10000000;
	while ((length = kseq_read(seq)) >= 0) {
		if (length >= min_seq_length) {
			uint32_t gi;
			sscanf(seq->name.s,"gi|%u|",&gi);
			uint32_t tax = vsearch(gi_vec,tax_vec,gi);
			if (tax > 0) {
				vector<uint32_t>::iterator it = lower_bound(nodes_tax.begin(),nodes_tax.end(),tax);
				if (*it == tax) {
					n_seqs++;
					uint32_t offset = uint32_t(it - nodes_tax.begin());
					*fasta[offset] << ">" << seq->name.s << " " << seq->comment.s << endl << seq->seq.s << endl;
					// check if we should write that stringbuf to a file
					long pos = fasta[offset]->tellp();
					if (pos > max_buffer) {
						write_fasta(fasta[offset],tax,nodes_tax,nodes_parent);
						fasta[offset]->str(string()); // empty the stringstream
					}
				}
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);

	// write the remaining sequences out
	for(uint32_t i=0; i<nodes_tax.size(); ++i)
		if (fasta[i]->tellp() > 0)
			write_fasta(fasta[i],nodes_tax[i],nodes_tax,nodes_parent);

	end = clock();
	fprintf(stderr,"processing %u sequences took %f seconds\n", n_seqs,
		double(end-start)/CLOCKS_PER_SEC);

	return 0;
}
