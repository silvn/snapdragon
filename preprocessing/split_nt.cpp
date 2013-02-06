#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include "kseq.h"
#include <ctime>
#include <iostream>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
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

	// do a dfs traversal of the taxonomy tree
	// and create the directories on the file system

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
	fprintf(stderr,"parsing nodes.dmp into two vectors took %f seconds\n",
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

	// split the nt.gz according to taxonomy
	start = clock();
	kseq_t *seq;
//	map<uint32_t,char*> tax_path;
	map<uint32_t,FILE*> active;
	map<uint32_t, uint32_t> active_used;
	int max_open_files = 8192;
	
	int length;
	gzFile fp = gzopen(argv[3], "r");
	seq = kseq_init(fp);
	uint32_t n_seqs=0;
	while ((length = kseq_read(seq)) >= 0) {
		n_seqs++;
		uint32_t gi;
		sscanf(seq->name.s,"gi|%u|",&gi);
		uint32_t tax = vsearch(gi_vec,tax_vec,gi);
		if (tax == 0) {
			fprintf(stderr,"failed to determine tax_id of gi %u\n",gi);
		}
		else if (vsearch(nodes_tax,nodes_parent,tax) > 0) {
//			fprintf(stderr,"gi|%u|tax|%u|len|%i\n",gi,tax,length);
			// check if this tax id is an active file
			map<uint32_t,FILE*>::iterator it1 = active.find(tax);
			if (it1 == active.end()) {
				// not an active file
				// check if we have already prepared a path
				// map<uint32_t,char*>::iterator it2 = tax_path.find(tax);
				// if (it2 == tax_path.end()) {
	//				fprintf(stderr,"creating path to tax\n");
					// is this map getting too large?
					// need to create the path
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
					tax = path.front();
					char path_str[1024] = "./";
					for(int i=path.size()-1;i>=0;i--) {
						sprintf(path_str,"%s%u/",path_str,path[i]);
	//					fprintf(stderr,"mkdir(%s,0755)\n",path_str);
						mkdir(path_str,0755);
					}
					sprintf(path_str,"%ssplit_nt.fa",path_str);
	//				fprintf(stderr,"file: %s\n",path_str);
//					tax_path[tax] = path_str;
//				}
				if (active.size() == max_open_files) {
					// have to close some files - how about the least recently used 64 files?
					vector<uint32_t> last_used;
					for(map<uint32_t,uint32_t>::iterator uit = active_used.begin(); uit != active_used.end(); ++uit) {
						last_used.push_back(uit->second);
					}
					sort(last_used.begin(),last_used.end());
					uint32_t mid_ru = last_used[max_open_files/8];
					vector<uint32_t> toclose;
					for(map<uint32_t,uint32_t>::iterator uit = active_used.begin(); uit != active_used.end(); ++uit)
						if(uit->second < mid_ru)
							toclose.push_back(uit->first);
					for(vector<uint32_t>::iterator it3 = toclose.begin(); it3 != toclose.end(); ++it3) {
						fclose(active[*it3]);
						active.erase(*it3);
						active_used.erase(*it3);
					}
				}
				// open the file in append mode
//				active[tax] = fopen(tax_path[tax],"a");
				active[tax] = fopen(path_str,"a");
				if (active[tax] == NULL) {
					fprintf(stderr,"failed to open %s in append mode\n",path_str);
				}
			}
			// write this seq to active[tax]
			fprintf(active[tax],">%s\n%s\n",seq->name.s,seq->seq.s);
			active_used[tax] = n_seqs;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	// close active filehandles
	for(map<uint32_t,FILE*>::iterator fpi = active.begin(); fpi != active.end(); ++fpi) {
		fclose(fpi->second);
	}
	end = clock();
	fprintf(stderr,"processing %u sequences took %f seconds\n", n_seqs,
		double(end-start)/CLOCKS_PER_SEC);
	return 0;
}
