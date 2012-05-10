#include "ibis.h"

// local data types
typedef std::set< const char*, ibis::lessi > qList;

// printout the usage string
static void usage(const char* name) {
  std::cout << "usage:\n" << name << std::endl
	    << "[-d1 directory_containing_a_dataset] " << std::endl
	    << "[-s1 select string]" << std::endl
	    << "[-w1 where-clause]" << std::endl
	    << "[-d2 directory_containing_a_dataset] " << std::endl
	    << "[-s2 select string]" << std::endl
	    << "[-w2 where-clause]" << std::endl
		<< "[-p threads]" << std::endl
} // usage

// function to parse the command line arguments
static void parse_args(int argc, char** argv,
	ibis::table*& tbl1, const char*& qcnd1, const char*& sel1,
	ibis::table*& tbl2, const char*& qcnd2, const char*& sel2, int* threads) {
  
  std::vector<const char*> dirs1;
  std::vector<const char*> dirs2;

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
			  ++ i;
			  if (argv[i][2] == '1')
				  dirs1.push_back(argv[i]);
			  else
				  dirs1.push_back(argv[i]);
			}
			break;
	      case 's':
			if (i+1 < argc) {
			  ++i;
			  if (argv[i][2] == '1')
				  sel1 = argv[i];
			  else
				  sel2 = argv[i];
			}
			break;
	      case 'w':
	      case 'W':
			if (i+1 < argc) {
			  ++ i;
			  if (argv[i][2] == '1')
				  qcnd1 = argv[i];
			  else
				  qcnd2 = argv[i];
			}
			break;
      } // switch (argv[i][1])
    } // normal arguments
  } // for (inti=1; ...)

  tbl1 = ibis::table::create(0);
  for (std::vector<const char*>::const_iterator it = dirs1.begin();
       it != dirs1.end(); ++ it) {
    if (tbl1 != 0)
      tbl1->addPartition(*it);
    else
      tbl1 = ibis::table::create(*it);
  }
  tbl2 = ibis::table::create(0);
  for (std::vector<const char*>::const_iterator it = dirs2.begin();
       it != dirs2.end(); ++ it) {
    if (tbl2 != 0)
      tbl2->addPartition(*it);
    else
      tbl2 = ibis::table::create(*it);
  }
  if (tbl1 == 0 || tbl2 == 0) {
    usage(argv[0]);
    exit(-2);
  }
} // parse_args

char* listToStr(ibis::array_t<const char*> list) {
	char* str;
	for(std::vector<const char*>::const_iterator it = list.begin();
	it != list.end(); ++ it) {
		if (str != 0) {
			sprintf(str,"%s,%s",*str,*it);
		} else {
			str = *it;
		}
	}
	return str;
}

int main(int argc, char** argv) {
    ibis::table* tbl1 = 0;
    const char* qcnd1=0;
    const char* sel1=0;
    ibis::table* tbl2 = 0;
    const char* qcnd2=0;
    const char* sel2=0;
	int threads=1;
    parse_args(argc, argv, tbl1, qcnd1, sel1, tbl2, qcnd2, sel2, &threads);

	if ((qcnd1 == 0 || *qcnd1 == 0))
		qcnd1 = "1=1";
	if ((qcnd2 == 0 || *qcnd2 == 0))
		qcnd2 = "1=1";

	if ((sel1 == 0 || *sel1 == 0))
		sel1 = listToStr(tbl1->columnNames());
	if ((sel2 == 0 || *sel2 == 0))
		sel2 = listToStr(tbl2->columnNames());

	std::cout << sel1 << std::endl;
	std::cout << sel2 << std::endl;

	std::vector<const ibis::part*> parts1;
	tbl1->getPartitions(parts1);

	


	std::vector<const ibis::part*> parts2;
	tbl2->getPartitions(parts2);

    delete tbl1;
	delete tbl2;
    return 0;
} // main
