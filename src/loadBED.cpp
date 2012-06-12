#include "ibis.h"
#include "bord.h"	// ibis::bord, ibis::table::bufferList
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <boost/thread.hpp>

using namespace std;


int threads=1;
const char *inBED;
const char *outdir;

vector<const char *> BEDlines;
vector<char *> BEDchrs;
vector<int> BEDstarts;
vector<int> BEDends;
vector<char *> BEDlabels;
vector<int> BEDscores;
vector<char> BEDstrands;

static void usage(const char* name)
{
  std::cout << "usage:\n" << name << std::endl
	    << "[-i input BED file] " << std::endl
	    << "[-o output directory]" << std::endl
		<< "[-p parallelize using p threads]" << std::endl;
}

void parse_args(int argc, char** argv)
{
	for (size_t i=1; i<argc; ++i) {
		if (*argv[i] == '-') { // arguments starting with -
			switch (argv[i][1]) {
				default:
				case 'h':
					usage(*argv);
					exit(0);
				case 'p':
					if (i+1 < argc)
						threads= atoi(argv[++i]);
					break;
				case 'i':
					if (i+1 < argc)
						inBED = argv[++i];
					break;
				case 'o':
					if (i+1 < argc)
						outdir = argv[++i];
					break;
			}
		}
	}
	if (threads <= 1)
		threads = 1;
}

void parseChunk(int from, int to)
{
	for(int i=from;i<to;i++)
	{
//		cout << "parsing line " << i << " " << BEDlines[i] << endl;
		BEDchrs[i] = new char[50];
		BEDlabels[i] = new char[50];
		sscanf(BEDlines[i],"%s\t%d\t%d\t%s\t%d\t%s",BEDchrs[i],
			&BEDstarts[i],&BEDends[i],BEDlabels[i],&BEDscores[i],&BEDstrands[i]);
	
/*
		if(CHRmap.count(chr)==0) {
//			cout << "starting new chrvec for chr " << chr << endl;
			vector<int> *chrvec = new vector<int>;
			CHRmap.insert(pair<string,vector<int>*>(chr,chrvec));
		}
		CHRmap.find(chr)->second->push_back(i);
*/
	}
}



void createPartition(string chr, vector<int>* ind) {
	int nrows = ind->size();
	size_t ncols = 5;
	ibis::table::bufferList tbuff(ncols,0);
	ibis::table::typeList ttypes(ncols);
	ibis::table::stringList colnames(ncols);
	IBIS_BLOCK_GUARD(ibis::table::freeBuffers, ibis::util::ref(tbuff), ibis::util::ref(ttypes));

	colnames[0]="start";
	colnames[1]="end";
	colnames[2]="label";
	colnames[3]="score";
	colnames[4]="strand";
	ttypes[0]=ibis::UINT;
	ttypes[1]=ibis::UINT;
	ttypes[2]=ibis::TEXT;
	ttypes[3]=ibis::UINT;
	ttypes[4]=ibis::BYTE;
	for(size_t i=0; i<ncols; i++)
		tbuff[i] = ibis::table::allocateBuffer(ttypes[i],nrows);
	ibis::array_t<uint32_t>* starts = static_cast<ibis::array_t<uint32_t>*>(tbuff[0]);
	ibis::array_t<uint32_t>* ends = static_cast<ibis::array_t<uint32_t>*>(tbuff[1]);
	vector<string>* labels = static_cast<vector<string>*>(tbuff[2]);
	ibis::array_t<uint32_t>* scores = static_cast<ibis::array_t<uint32_t>*>(tbuff[3]);
	ibis::array_t<signed char>* strands = static_cast<ibis::array_t<signed char>*>(tbuff[4]);
	for(size_t i=0;i<nrows;i++) {
		(*starts)[i]  = BEDstarts[(*ind)[i]];
		(*ends)[i]    = BEDends[(*ind)[i]];
		(*labels)[i]  = BEDlabels[(*ind)[i]];
		(*scores)[i]  = BEDscores[(*ind)[i]];
		(*strands)[i] = BEDstrands[(*ind)[i]];
	}
	
   	ibis::table *tbl = new ibis::bord("fromBED", "fromBED", nrows, tbuff, ttypes, colnames, &colnames, 0);
	// need to sort by start,end
	tbl->orderby("start,end");
	string chrdir(outdir);
	chrdir += FASTBIT_DIRSEP;
	chrdir += chr;
	ibis::part tmp(chrdir.c_str(), static_cast<const char*>(0));

	std::string mdfile = chrdir;
	mdfile += FASTBIT_DIRSEP;
	mdfile += "-part.txt";
	std::ofstream md(mdfile.c_str());
	if(! md) {
		std::cerr << "failed to open metadata file " << mdfile << std::endl;
		exit(-3);
	}

	md << "BEGIN HEADER\nName = parsed" << "\nDescription = parsed"
	   << "\nNumber_of_rows = " << nrows
	   << "\nNumber_of_columns = 5"
	   << "\nmetaTags = FBchr = " << chr;
	md << "\nEND HEADER\n";
	for(size_t i=0; i<ncols; i++) {
		md << "\nBegin Column\nname = " << colnames[i] << "\ndata_type = "
		   << ibis::TYPESTRING[(int) ttypes[i]];
	    if(i==0)
			md << "\nsorted = true";
		md << "\nEnd Column\n";

		string cnm = chrdir;
		cnm += FASTBIT_DIRSEP;
		cnm += colnames[i];
		int fdes = UnixOpen(cnm.c_str(), OPEN_WRITEADD, OPEN_FILEMODE);
		if (fdes < 0) {
			cerr << "failed to open file " << cnm << " for writing" << endl;
			exit(-4);
		}
		IBIS_BLOCK_GUARD(UnixClose, fdes);
#if defined(_WIN32) && defined(_MSC_VER)
		(void)_setmode(fdes, _O_BINARY);
#endif
		off_t pos = UnixSeek(fdes,0,SEEK_END);
		switch(ttypes[i]) {
			case ibis::BYTE:
			{
				char* vals = new char[nrows];
				uint64_t ierr = tbl->getColumnAsBytes(colnames[i],vals);
				pos = UnixWrite(fdes,vals,sizeof(char)*nrows);
				delete vals;
				break;
			}
			case ibis::UINT:
			{
				uint32_t* vals = new uint32_t[nrows];
				uint64_t ierr = tbl->getColumnAsUInts(colnames[i],vals);
				pos = UnixWrite(fdes,vals,sizeof(uint32_t)*nrows);
				delete vals;
				break;
			}
			case ibis::TEXT:
			{
				vector<string>* vals = new vector<string>();
				uint64_t ierr = tbl->getColumnAsStrings(colnames[i],*vals);
				for (uint32_t j=0; j< nrows; ++j)
					pos = UnixWrite(fdes,(*vals)[j].c_str(),sizeof(char)*((*vals)[j].length()+1));
				delete vals;
				break;
			}
		}
#if defined(FASTBIT_SYNC_WRITE)
#if _POSIX_FSYNC+0 > 0
		(void) UnixFlush(fdes); // write to disk
#elif defined(_WIN32) && defined(_MSC_VER)
		(void) _commit(fdes);
#endif
#endif
		
	}
	md.close();

}

int main(int argc, char** argv)
{
    parse_args(argc, argv);
	ifstream BEDfile (inBED);
	if (BEDfile.is_open())
	{
		while (BEDfile.good())
		{
			string *line = new string;
			getline(BEDfile,*line);
			BEDlines.push_back(line->c_str());
		}
		BEDfile.close();
	}
	else cout << "unable to open file '" << inBED << endl;

	int nlines = BEDlines.size();
	cout << "finished reading " << nlines << " lines in to memory" << endl;

	BEDchrs.resize(nlines);
	BEDstarts.resize(nlines);
	BEDends.resize(nlines);
	BEDlabels.resize(nlines);
	BEDscores.resize(nlines);
	BEDstrands.resize(nlines);
	
	int chunksize = BEDlines.size()/threads;
	
	boost::thread_group tg;
	for(int i=0;i<nlines;i+=chunksize) {
		int j = i+chunksize;
		if (j>nlines)
			j=nlines;
		cout << "starting chunk [" << i << "," << j << "]" << endl;
		tg.create_thread(boost::bind(parseChunk, i, j));
	}
	tg.join_all();

	cout << "finished parsing BEDlines" << endl;

	// index the rows by chromosome
	map<string,vector<int>*> CHRmap;
	vector<int>* chrvec;
	string prev="not likely";
	for(int i=0;i<nlines;i++) {
		string chr(BEDchrs[i]);
		if (chr.length()>0) {
			if (chr != prev) {
				if (CHRmap.count(chr) == 0) {
					chrvec = new vector<int>;
					CHRmap.insert(pair<string,vector<int>*>(chr,chrvec));
				} else {
					chrvec = CHRmap.find(chr)->second;
				}
				prev = chr;
			}
			chrvec->push_back(i);
		}
	}
	cout << "finished with CHRmap" << endl;

	// create a partition for each chromosome
	map<string,vector<int>*>::iterator it;
	int tg_size=0;
	for ( it=CHRmap.begin(); it != CHRmap.end(); it++ ) {
		cout << (*it).first << " => " << (*it).second->size() << endl;
		if (tg_size == threads) {
			tg.join_all();
			tg_size=0;
		}
		tg.create_thread(boost::bind(createPartition, (*it).first, (*it).second));
		tg_size++;
	}
	tg.join_all();
	
	return 0;
}