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
static void usage(const char* name)
{
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
}

// function to parse the command line arguments
void parse_args(int argc, char** argv)
{
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
}

// Build a table to represent the interval-join of two partitions
// using two bitvector masks and two arrays of offsets into the matching rows
// New columns for relative start and end are included (possibly binned)
// In the end the table needs to be sorted by these columns
// When Bstart == Bend you don't add an end column
ibis::table* fillResult(const ibis::part* Apart, const ibis::part* Bpart,
	ibis::bitvector Amatch, ibis::bitvector Bmatch,
	std::vector<uint32_t>& Aidx, std::vector<uint32_t>& Bidx)
{
	
}

// fetch intervals from the partitions to be joined
// identify features that overlap
// set bits in Amatch and Bmatch so we can fetch only elements that overlap
// fill two arrays of indexes that say what rows to take from selected columns of A and B
void Stacker(const ibis::part* Apart, const ibis::part* Bpart,
	ibis::bitvector Amask, ibis::bitvector Bmask, int before, int after)
{
	if (Amask.cnt() == 0 || Bmask.cnt() == 0)
		return;

	ibis::column *Astart_col = Apart->getColumn(Astart);
	ibis::column *Aend_col = Apart->getColumn(Aend);
	ibis::column *Bstart_col = Bpart->getColumn(Bstart);
	ibis::column *Bend_col = Bpart->getColumn(Bend);

	ibis::array_t<ibis::rid_t> *Arids = Apart->getRIDs(Amask);
	ibis::array_t<ibis::rid_t> *Brids = Bpart->getRIDs(Bmask);
	ibis::array_t<uint32_t> *Astart_val = Astart_col->selectUInts(Amask);
	ibis::array_t<uint32_t> *Aend_val = Aend_col->selectUInts(Amask);
	ibis::array_t<uint32_t> *Bstart_val = Bstart_col->selectUInts(Bmask);
	ibis::array_t<uint32_t> *Bend_val = Bend_col->selectUInts(Bmask);
	
	ibis::bitvector Amatch;
	ibis::bitvector Bmatch;
	const uint32_t nA = Amask.cnt();
	const uint32_t nB = Bmask.cnt();
	std::vector<uint32_t> Aidx;
	std::vector<uint32_t> Bidx;
	uint32_t iA=0;
	uint32_t iB=0;
	uint32_t mA=0;
	uint32_t mB=0;
	std::cerr << "before Amask.cnt() = " << Amask.cnt() << std::endl;
	std::cerr << "before Bmask.cnt() = " << Bmask.cnt() << std::endl;
	while (iA < nA && iB < nB) {
		uint32_t before_ = (*Astart_val)[iA] - before;
		uint32_t after_ = (*Aend_val)[iA] + after;
		if ((*Bend_val)[iB] < before_) // B comes first
			iB++;
		else if (after_ < (*Bstart_val)[iB]) // A comes first
			iA++;
		else {
			// found a match
			Aidx.push_back(mA);
			Bidx.push_back(mB);
			Amatch.setBit((const uint32_t) (*Arids)[iA].value, 1);
			Bmatch.setBit((const uint32_t) (*Brids)[iB].value, 1);
			mA = Amatch.cnt();
			mB = Bmatch.cnt();
			// report all the matches of B within this A interval
			// then increment iA
			uint32_t iiB=iB+1;
			while (iiB < nB && (*Bstart_val)[iiB] <= after_) {
				if ((*Bend_val)[iiB] >= before_) {
					// overlap
					Aidx.push_back(mA);
					Bidx.push_back(mB);
					Bmatch.setBit((const uint32_t) (*Brids)[iiB].value, 1);
					mB = Bmatch.cnt();
				}
				iiB++;
			}
			iA++;
		}
	}

	std::cerr << "after Amatch.cnt() = " << Amatch.cnt() << std::endl;
	std::cerr << "after Bmatch.cnt() = " << Bmatch.cnt() << std::endl;

	if (Amatch.cnt() == 0 || Bmatch.cnt() == 0)
		return;

	ibis::table *jtable = fillResult(Apart,Bpart,Amatch,Bmatch,Aidx,Bidx);

}

// apply the user supplied query expression to the given part
// populate a bitvector mask that corresponds to the matching rows
// return the number of rows that match
int countHits(const ibis::part* part, ibis::bitvector* mask,
	 const ibis::qExpr* cond, const char* colName)
{
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
	return mask->cnt();
}

// populate bitvectors for hits to the sense and antisense strands
// Behold the beauty of WAH bitwise operations!
void splitByStrand(const ibis::part* part, ibis::bitvector* mask, ibis::bitvector* plus, ibis::bitvector* minus)
{
	ibis::countQuery que(part);
	int ierr = que.setWhereClause(senseExpr);
	ierr = que.evaluate();
	plus->copy(*que.getHitVector());
	minus->copy(*que.getHitVector());
	minus->flip();
	*plus &= *mask;
	*minus &= *mask;
}

// intersecting genomic intervals are always on the same chromsome
// so we operate on partitions with the same FBchr
// other command line options control how intervals are selected
// so we have to handle special cases related to strand
void setupStacker(const ibis::part* Apart, const ibis::part* Bpart)
{
	// TODO: do some checks
	// count hits to A and B on this chr
	// make masks based on Aqcnd and Bqcnd
    ibis::bitvector Amask;
	int Arows = countHits(Apart, &Amask, Acond, Astart);	
	if (Arows == 0) return;
	
    ibis::bitvector Bmask;
	int Brows = countHits(Bpart, &Bmask, Bcond, Bstart);
	if (Brows == 0) return;
	
	if(same_strand == 1 && (Anaty.count("strand") > 0) && (Bnaty.count("strand") > 0)) {
		// split into subproblems for each strand
		// update the masks (Aplus, Aminus, Bplus, Bminus)
		ibis::bitvector Aplus;
		ibis::bitvector Aminus;
		splitByStrand(Apart, &Amask, &Aplus, &Aminus);
		
		ibis::bitvector Bplus;
		ibis::bitvector Bminus;
		splitByStrand(Bpart, &Bmask, &Bplus, &Bminus);

		Stacker(Apart, Bpart, Aplus, Bplus, left, right);
		if (stranded_windows>0) {
			Stacker(Apart, Bpart, Aminus, Bminus, right, left);
		}
		else {
			Stacker(Apart, Bpart, Aminus, Bminus, left, right);
		}
	}
	else {
		if (stranded_windows>0 && left != right and Anaty.count("strand")>0) {
			// split into subproblems
			// update the masks (Aplus, Aminus)
			ibis::bitvector Aplus;
			ibis::bitvector Aminus;
			splitByStrand(Apart, &Amask, &Aplus, &Aminus);
			Stacker(Apart, Bpart, Aplus, Bmask, left, right);
			Stacker(Apart, Bpart, Aminus, Bmask, right, left);
		}
		else {
			Stacker(Apart, Bpart, Amask, Bmask, left, right);
		}
	}
}

int main(int argc, char** argv)
{
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
	std::vector<const ibis::part*> Bparts;
	B->getPartitions(Bparts);

	std::map<const char*,const ibis::part*> Bpartmap;	
	for(size_t i=0; i<Bparts.size(); ++ i)
		Bpartmap.insert(std::pair<const char*,const ibis::part*>(Bparts[i]->getMetaTag("FBchr"),Bparts[i]));
	std::vector<const ibis::part*> Aparts;
	A->getPartitions(Aparts);
	boost::thread_group g;
	for(size_t i=0; i<Aparts.size(); ++ i) {
		const char* chr = Aparts[i]->getMetaTag("FBchr");
		if (Bpartmap.count(chr)>0) {
			if (parallelize > 0)
				g.create_thread(boost::bind(setupStacker, Aparts[i], Bpartmap.find(chr)->second));
			else
				setupStacker(Aparts[i],Bpartmap.find(chr)->second);
		}
	}
	if (parallelize > 0)
		g.join_all();

	// concatenate each part into one table?
	// orderby()?

	// write the table out - user supplied outdir
	// create one subdirectory partition with FBchr = stacked
	
	delete A;
	delete B;
    return 0;
}
