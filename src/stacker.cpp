#include "ibis.h"
#include <map>
#include <string>
#include <boost/thread.hpp>

// had to make these global
const char* Afrom = 0;
const char* Aqcnd="1=1";
const char* Asel="";
const char* Astart="start";
const char* Aend = "end";
ibis::table::stringList Anames;
ibis::table::typeList Atypes;
std::map<const char*,ibis::TYPE_T> Anaty;
ibis::qExpr* Acond=0;

const char* Bfrom = 0;
const char* Bqcnd="1=1";
const char* Bsel="";
const char* Bstart="start";
const char* Bend = "end";
ibis::table::stringList Bnames;
ibis::table::typeList Btypes;
std::map<const char*,ibis::TYPE_T> Bnaty;
ibis::qExpr* Bcond=0;

int left=1000;
int right=1000;
int same_strand=0;
int stranded_windows=0;
int bins=0;

int parallelize=0;

ibis::whereClause senseWhere = ibis::whereClause("strand == 1");
ibis::qExpr* senseExpr = senseWhere.getExpr();


// printout the usage string
static void usage(const char* name) {
  std::cout << "usage:\n" << name << std::endl
	    << "[-d1 directory containing A dataset] " << std::endl
	    << "[-c1 columns from A]" << std::endl
	    << "[-w1 where-clause for A]" << std::endl
		<< "[-s1 start column from A]" << std::endl
		<< "[-e1 end column from A]" << std::endl
	    << "[-d2 directory containing B dataset] " << std::endl
	    << "[-c2 columns from B]" << std::endl
	    << "[-w2 where-clause for B]" << std::endl
		<< "[-s2 start column from B]" << std::endl
		<< "[-e2 end column from B]" << std::endl
		<< "[-l bases to the left of A features]" << std::endl
		<< "[-r bases to the right of A features]" << std::endl
		<< "[-sm only report hits in B that overlap A on the same strand]" << std::endl
		<< "[-sw define -l and -r based on strand]" << std::endl
		<< "[-b number of bins to use for normalization after stacking. default: don't normalize]" << std::endl
		<< "[-p parallelize using threads]" << std::endl;
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
		  case 'p':
			  parallelize=1;
			  break;
  		  case 'b':
  			  if(i+1 < argc) {
  				  ++i;
  				  bins = atoi(argv[i]);
  			  }
  			  break;
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
			  if (argv[i][2] == '1')
				  Afrom = argv[++i];
			  else
				  Bfrom = argv[++i];
			}
			break;
  	      case 'c':
  			if (i+1 < argc) {
  				if (argv[i][2] == '1')
					Asel = argv[++i];
				else
					Bsel = argv[++i];
				break;
			}
			break;
	      case 'e':
			if (i+1 < argc) {
				if (argv[i][2] == '1')
					Aend = argv[++i];
				else
					Bend = argv[++i];
  				break;
  			}
  			break;
	      case 's':
			switch (argv[i][2]) {
				case '1':
					Astart = argv[++i];
					break;
				case '2':
					Bstart = argv[++i];
					break;
				case 'm':
					same_strand = 1;
					break;
				case 'w':
					stranded_windows = 1;
					break;
			}
			break;
	      case 'w':
	      case 'W':
			if (i+1 < argc) {
			  if (argv[i][2] == '1')
				  Aqcnd = argv[++i];
			  else
				  Bqcnd = argv[++i];
			}
			break;
      } // switch (argv[i][1])
    } // normal arguments
  } // for (inti=1; ...)

  std::cerr << *argv << " -d1 " << Afrom << " -d2 " << Bfrom
	  << " -s1 " << Astart << " -s2 " << Bstart
	  << " -e1 " << Aend << " -e2 " << Bend
	  << " -w1 " << Aqcnd << " -w2 " << Bqcnd
	  << " -c1 " << Asel << " -c2 " << Bsel
	  << " -r " << right << " -l " << left
	  << " -b " << bins << " -p " << parallelize
	  << " -sm " << same_strand << " -sw " << stranded_windows
	  << std::endl;
} // parse_args

void Stacker(ibis::part* Apart, ibis::part* Bpart, ibis::bitvector Amask, ibis::bitvector Bmask, int before, int after) {
	// by the time this function is called it doesn't have to worry about
	// strand, just intervals
	
	// fetch [Astart-before,Aend+after]
	// fetch [Bstart,Bend]
	// iterate and report overlaps
	// consult jRange.cpp
	ibis::column *Astart_col = Apart->getColumn(Astart);
	ibis::column *Aend_col = Apart->getColumn(Aend);
	ibis::column *Bstart_col = Bpart->getColumn(Bstart);
	ibis::column *Bend_col = Bpart->getColumn(Bend);

	ibis::array_t<unsigned int> *Astart_val = Astart_col->selectUInts(Amask);
	ibis::array_t<unsigned int> *Aend_val = Aend_col->selectUInts(Amask);
	ibis::array_t<unsigned int> *Bstart_val = Bstart_col->selectUInts(Bmask);
	ibis::array_t<unsigned int> *Bend_val = Bend_col->selectUInts(Bmask);

	std::cerr << "stacker() Astart_val: " << Astart_val->size() << std::endl;
	std::cerr << "stacker() Aend_val: " << Aend_val->size() << std::endl;
	std::cerr << "stacker() Bstart_val: " << Bstart_val->size() << std::endl;
	std::cerr << "stacker() Bend_val: " << Bend_val->size() << std::endl;

}

void getHits(ibis::part* part, ibis::bitvector* mask,
	 const ibis::qExpr* cond, const char* colName) {
	if (cond != 0) {
		ibis::countQuery que(part);
		int ierr = que.setWhereClause(cond);
		ierr = que.evaluate();
		mask->copy(*que.getHitVector());
	}
	else {
		ibis::column *c = part->getColumn(colName);
		c->getNullMask(*mask);
	}
}

void splitByStrand(ibis::part* part, ibis::bitvector* mask, ibis::bitvector* plus, ibis::bitvector* minus) {
	ibis::countQuery que(part);
	int ierr = que.setWhereClause(senseExpr);
	ierr = que.evaluate();
	plus->copy(*que.getHitVector());
	minus->copy(*que.getHitVector());
	minus->flip();
	*plus &= *mask;
	*minus &= *mask;
}

void setupStacker(const char* chr) {
	// do some checks
	// count hits to A and B on this chr
	// make masks based on Aqcnd and Bqcnd
 	std::string Apart_dir = Afrom;
 	Apart_dir.append("/");
 	Apart_dir.append(chr);
 	ibis::part Apart(Apart_dir.c_str(),false);
    ibis::bitvector Amask;
	getHits(&Apart, &Amask, Acond, Astart);
	
	std::cerr << "Amask.cnt() = " << Amask.cnt() << std::endl;
	if (Amask.cnt() == 0) return;
	
 	std::string Bpart_dir = Bfrom;
 	Bpart_dir.append("/");
 	Bpart_dir.append(chr);
 	ibis::part Bpart(Bpart_dir.c_str(),false);
    ibis::bitvector Bmask;
	getHits(&Bpart, &Bmask, Bcond, Bstart);

	std::cerr << "Bmask.cnt() = " << Bmask.cnt() << std::endl;
	if (Bmask.cnt() == 0) return;
	
	if(same_strand == 1 && (Anaty.count("strand") > 0) && (Bnaty.count("strand") > 0)) {
		std::cerr << "1" << std::endl;
		// split into subproblems for each strand
		// update the masks (Aplus, Aminus, Bplus, Bminus)
		ibis::bitvector Aplus;
		ibis::bitvector Aminus;
		splitByStrand(&Apart, &Amask, &Aplus, &Aminus);
		
		ibis::bitvector Bplus;
		ibis::bitvector Bminus;
		splitByStrand(&Bpart, &Bmask, &Bplus, &Bminus);

		Stacker(&Apart, &Bpart, Aplus, Bplus, left, right);
		if (stranded_windows>0) {
			Stacker(&Apart, &Bpart, Aminus, Bminus, right, left);
		}
		else {
			Stacker(&Apart, &Bpart, Aminus, Bminus, left, right);
		}
	}
	else {
		std::cerr << "2" << std::endl;
		if (stranded_windows>0 && left != right and Anaty.count("strand")>0) {
			// split into subproblems
			// update the masks (Aplus, Aminus)
			ibis::bitvector Aplus;
			ibis::bitvector Aminus;
			splitByStrand(&Apart, &Amask, &Aplus, &Aminus);
			Stacker(&Apart, &Bpart, Aplus, Bmask, left, right);
			Stacker(&Apart, &Bpart, Aminus, Bmask, right, left);
		}
		else {
			Stacker(&Apart, &Bpart, Amask, Bmask, left, right);
		}
	}
}

int main(int argc, char** argv) {
    parse_args(argc, argv);

	// parse qExpr once so it can be reused in each partition
	ibis::whereClause Awhere = ibis::whereClause(Aqcnd);
	Acond = Awhere.getExpr();
	ibis::whereClause Bwhere = ibis::whereClause(Bqcnd);
	Bcond = Bwhere.getExpr();

	// populate column name => type map
	ibis::table* A = ibis::table::create(Afrom);
	Anames = A->columnNames();
	Atypes = A->columnTypes();
    for (size_t i = 0; i < Anames.size(); ++ i) {
		Anaty.insert(std::pair<const char*,ibis::TYPE_T>(Anames[i],Atypes[i]));
	}	
	ibis::table* B = ibis::table::create(Bfrom);
	Bnames = B->columnNames();
	Btypes = B->columnTypes();
    for (size_t i = 0; i < Bnames.size(); ++ i) {
		Bnaty.insert(std::pair<const char*,ibis::TYPE_T>(Bnames[i],Btypes[i]));
	}

	// fetch all the FBchr names from table A
	ibis::table* FBchr = A->select("FBchr,count(*)",Acond);
	const size_t nr = static_cast<size_t>(FBchr->nRows());
	std::vector<std::string>* chrList = new std::vector<std::string>();
	int64_t ierr = FBchr->getColumnAsStrings("FBchr", *chrList);
	if (ierr < 0 || ((size_t) ierr) < nr) {
		return -1;
	}

	// for each FBchr, build a partition with the stacked features
	boost::thread_group g;
	for(int i=0; i<nr; ++i) {
		std::cerr << "calling setupStacker(" << (*chrList)[i] << ")" << std::endl;
		if (parallelize > 0)
			g.create_thread(boost::bind(setupStacker, (*chrList)[i].c_str()));
		else
			setupStacker((*chrList)[i].c_str());
	}
	if (parallelize > 0)
		g.join_all();
	// concatenate each part into one table?
	// orderby()

	// write the table out - need to generate a FBchr metaTag based on FBset l, r, sw, sm
	
	delete FBchr;
	delete A;
	delete B;
    return 0;
} // main
