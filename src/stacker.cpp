#include "ibis.h"
#include <string>
#include <boost/thread.hpp>
using namespace std;

// printout the usage string
static void usage(const char* name) {
  cout << "usage:\n" << name << endl
	    << "[-d1 directory_containing_a_dataset] " << endl
	    << "[-s1 select string]" << endl
	    << "[-w1 where-clause]" << endl
	    << "[-d2 directory_containing_a_dataset] " << endl
	    << "[-s2 select string]" << endl
	    << "[-w2 where-clause]" << endl
		<< "[-p threads]" << endl;
} // usage

// function to parse the command line arguments
static void parse_args(int argc, char** argv,
	const char*& tbl1, const char*& qcnd1, const char*& sel1,
	const char*& tbl2, const char*& qcnd2, const char*& sel2, int* threads) {
  
  for (int i=1; i<argc; ++i) {
    if (*argv[i] == '-') { // normal arguments starting with -
      switch (argv[i][1]) {
	      default:
	      case 'h':
	      case 'H':
			usage(*argv);
			exit(0);
		  case 'p':
			  if(i+1 < argc) {
				  ++i;
				  *threads = atoi(argv[i]);
			  }
			  break;
	      case 'd':
	      case 'D':
			if (i+1 < argc) {
			  if (argv[i++][2] == '1')
				  tbl1 = argv[i];
			  else
				  tbl2 = argv[i];
			}
			break;
	      case 's':
			if (i+1 < argc) {
			  if (argv[i++][2] == '1')
				  sel1 = argv[i];
			  else
				  sel2 = argv[i];
			}
			break;
	      case 'w':
	      case 'W':
			if (i+1 < argc) {
			  if (argv[i++][2] == '1')
				  qcnd1 = argv[i];
			  else
				  qcnd2 = argv[i];
			}
			break;
      } // switch (argv[i][1])
    } // normal arguments
  } // for (inti=1; ...)

  if (tbl1 == 0 || tbl2 == 0) {
    usage(argv[0]);
    exit(-2);
  }
} // parse_args

const char* listToStr(ibis::array_t<const char*> list) {
	string str = list[0];
	for(size_t i = 1; i < list.size(); ++i) {
		str.append(",");
		str.append(list[i]);
	}
	return str.c_str();
}

ibis::table* pSelect(const char* sel, const char* from, const char* qcnd) {
	ibis::table* tbl = ibis::table::create(from);
	if ((qcnd == 0 || *qcnd == 0))
		qcnd = "1=1";
	if ((sel == 0 || *sel == 0))
		sel = listToStr(tbl->columnNames());
	cout << "select " << sel << " from " << tbl->name() << " where " << qcnd << endl;

	ibis::table* res = tbl->select(sel,qcnd);
	delete tbl;
	return res;
}

void pSelect(const char* sel, const char* from, const char* qcnd, int threads) {
	ibis::table* tbl = ibis::table::create(from);
	vector<const ibis::part*> parts;
	tbl->getPartitions(parts);
	for(size_t i= 0; i< parts.size(); ++i) {
		string part_dir = from;
		part_dir.append("/");
		part_dir.append(parts[i]->getMetaTag("FBchr"));
		ibis::table* part_res = pSelect(sel, part_dir.c_str(), qcnd);
		delete part_res;
	}
	
	delete tbl;
}

int main(int argc, char** argv) {
    const char* tbl1 = 0;
    const char* qcnd1=0;
    const char* sel1=0;
    const char* tbl2 = 0;
    const char* qcnd2=0;
    const char* sel2=0;
	int threads=1;
    parse_args(argc, argv, tbl1, qcnd1, sel1, tbl2, qcnd2, sel2, &threads);

	if (threads > 0) {
		pSelect(sel1,tbl1,qcnd1,threads);
		pSelect(sel2,tbl2,qcnd2,threads);
	} else {
		ibis::table* res1 = pSelect(sel1,tbl1,qcnd1);
		ibis::table* res2 = pSelect(sel2,tbl2,qcnd2);
		delete res1;
		delete res2;
	}
    return 0;
} // main
