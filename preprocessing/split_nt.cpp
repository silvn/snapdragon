#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include "kseq.h"
#include <ctime>
#include <iostream>
#include <string>
#include <errno.h>

using namespace std;

inline uint32_t vsearch(vector<uint32_t> &a, vector<uint32_t> &b, uint32_t v) {
	vector<uint32_t>::iterator it = lower_bound(a.begin(), a.end(), v);
	if (*it == v)
		return b[uint32_t(it - a.begin())];
	return 0;
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	if (argc != 3) {
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
	fprintf(stderr,"parsing gi_taxid_nucl.dmp into two vectors took %f seconds\n",
		double(end-start)/CLOCKS_PER_SEC);
	
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
	fprintf(stderr,"parsing nodes.dmp into two vectors took %s seconds\n",
		double(end-start)/CLOCKS_PER_SEC);

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

	return 0;
	// split the nt.gz according to taxonomy
	kseq_t *seq;
	map<uint32_t,string> tax_path;
	map<uint32_t,gzFile> active;
	int max_open_files = 128;
	
	int length;
	gzFile fp = gzopen(argv[3], "r");
	seq = kseq_init(fp);
	while ((length = kseq_read(seq)) >= 0) {
		uint32_t gi;
		sscanf(seq->name.s,">gi|%u|",&gi);
		uint32_t tax = vsearch(gi_vec,tax_vec,gi);

		// check if this tax id is an active file
		map<uint32_t,gzFile>::iterator it1 = active.find(tax);
		if (it1 == active.end()) {
			// not an active file
			// check if we have already prepared a path
			map<uint32_t,string>::iterator it2 = tax_path.find(tax);
			if (it2 == tax_path.end()) {
				// need to create the path
				vector<uint32_t> path;
				path.push_back(tax);
				while (tax != 1) {
					tax = vsearch(nodes_tax,nodes_parent,tax);
					path.push_back(tax);
				}
				string path_str;
				for(int i=path.size()-1;i>=0;i--)
					path_str += path[i] + "/";
				path_str += "nt.gz";
				tax_path[tax] = path_str;
			}
			if (active.size() == max_open_files) {
				// have to close a file
				
				// take it out of the active map
			}
			// open the file in append mode
		}
		// write this seq to active[tax]
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
