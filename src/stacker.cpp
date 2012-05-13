#include "ibis.h"
#include <string>
#include <boost/thread.hpp>

// had to make these global because boost::bind is limited to 9 args
const char* Afrom = 0;
const char* Aqcnd="1=1";
const char* Asel="start,end,strand";
const char* Bfrom = 0;
const char* Bqcnd="1=1";
const char* Bsel="start,end,strand";
int left=1000;
int right=1000;
int same_strand=0;
int stranded_windows=0;

// printout the usage string
static void usage(const char* name) {
  std::cout << "usage:\n" << name << std::endl
	    << "[-d1 directory containing A dataset] " << std::endl
	    << "[-s1 select string for A]" << std::endl
	    << "[-w1 where-clause for A]" << std::endl
	    << "[-d2 directory containing B dataset] " << std::endl
	    << "[-s2 select string for B]" << std::endl
	    << "[-w2 where-clause for B]" << std::endl
		<< "[-l bases to the left of A features]" << std::endl
		<< "[-r bases to the right of A features]" << std::endl
		<< "[-sm only report hits in B that overlap A on the same strand]" << std::endl
		<< "[-sw define -l and -r based on strand]" << std::endl;
} // usage

// function to parse the command line arguments
void parse_args(int argc, char** argv) {
  
  for (int i=1; i<argc; ++i) {
    if (*argv[i] == '-') { // normal arguments starting with -
      switch (argv[i][1]) {
	      default:
	      case 'h':
	      case 'H':
			usage(*argv);
			exit(0);
		  case 'l':
			  if(i+1 < argc) {
				  ++i;
				  left = atoi(argv[i]);
			  }
			  break;
		  case 'r':
			  if(i+1 < argc) {
				  ++i;
				  right = atoi(argv[i]);
			  }
			  break;
	      case 'd':
	      case 'D':
			if (i+1 < argc) {
			  if (argv[i++][2] == '1')
				  Afrom = argv[i];
			  else
				  Bfrom = argv[i];
			}
			break;
	      case 's':
			if (i+1 < argc) {
				i++;
				switch (argv[i][2]) {
					case '1':
						Asel = argv[i];
					case '2':
						Bsel = argv[i];
					case 'm':
						same_strand = 1;
					case 'w':
						stranded_windows = 1;
				}
			}
			break;
	      case 'w':
	      case 'W':
			if (i+1 < argc) {
			  if (argv[i++][2] == '1')
				  Aqcnd = argv[i];
			  else
				  Bqcnd = argv[i];
			}
			break;
      } // switch (argv[i][1])
    } // normal arguments
  } // for (inti=1; ...)
} // parse_args

/*

const char* listToStr(ibis::array_t<const char*> list) {
	std::string str = list[0];
	for(size_t i = 1; i < list.size(); ++i) {
		str.append(",");
		str.append(list[i]);
	}
	return str.c_str();
}

void pSelect(const char* sel, const char* from, const char* qcnd) {
	ibis::table* tbl = ibis::table::create(from);
	if ((qcnd == 0 || *qcnd == 0))
		qcnd = "1=1";
	if ((sel == 0 || *sel == 0))
		sel = listToStr(tbl->columnNames());

	ibis::table* res = tbl->select(sel,qcnd);
	delete tbl;
	delete res;
}

void pSelect(const char* sel, const char* from, const char* qcnd, int threads) {
	ibis::table* ftbl = ibis::table::create(from);
	if ((qcnd == 0 || *qcnd == 0))
		qcnd = "1=1";
	if ((sel == 0 || *sel == 0))
		sel = listToStr(ftbl->columnNames());
	boost::thread_group g;
	std::vector<const ibis::part*> parts;
	ftbl->getPartitions(parts);
	for(size_t i= 0; i< parts.size(); ++i) {
		std::string part_dir = from;
		part_dir.append("/");
		part_dir.append(parts[i]->getMetaTag("FBchr"));
		std::cout << "calling pSelect sel='" << sel <<"', from='" << part_dir.c_str() << "', qcnd='" << qcnd << "'" << std::endl;
//		pSelect(sel, part_dir.c_str(), qcnd);
		g.create_thread(boost::bind(pSelect, sel, part_dir.c_str(), qcnd)); 
	}
	g.join_all();
	
	delete ftbl;
}

*/

void Selector(ibis::table*& res, const char* path, const char* dir, const char* qcnd, const char* sel) {
	std::string part_dir = path;
	part_dir.append("/");
	part_dir.append(dir);
	ibis::table* tbl = ibis::table::create(part_dir.c_str());
	if (tbl == 0) {
		std::cerr << "missing partition: " << part_dir.c_str() << std::endl;
		return;
	}
	ibis::table* r = tbl->select(sel,qcnd);
	res = r;
	delete tbl;
}


void Stacker(const char* chr) {
	// fetch results from the FBchr matched pair of partitions from A and B
	ibis::table* A=0;
	boost::thread Athread(Selector, boost::ref(A), Afrom, chr, Aqcnd, Asel);
	ibis::table* B=0;
	Selector(B, Bfrom, chr, Bqcnd, Bsel);
	Athread.join();
	// iterate over A and B and do the join
	// generate a new part to hold the results 
	// sort the partition by start,end (here if doing multithreaded merge?)
}

int main(int argc, char** argv) {
    parse_args(argc, argv);

	ibis::table* A = ibis::table::create(Afrom);
	
	// for each FBchr, build a table/partition with the stacked features
	std::vector<const ibis::part*> parts;
	A->getPartitions(parts);
	boost::thread_group g;
	for(int i=0; i< parts.size(); ++i) {
		const char* chr = parts[i]->getMetaTag("FBchr");
		g.create_thread(boost::bind(Stacker, chr));
	}
	g.join_all();
	// merge the sorted result partitions using multiple threads?
	// write the table out - need to generate a FBchr name based on FBset l, r, sw, sm ?
	// if you didn't do the fancy sort just call orderby() before saving the table
	
    return 0;
} // main
