#include "ibis.h"
#include "bord.h"	// ibis::bord, ibis::table::bufferList
#include <map>
#include <string>
#include <boost/thread.hpp>

// had to make these global
const char* Afrom = 0;
const char* Aqcnd="1=1";
char* Asel;
const char* Astart="start";
const char* Aend = "end";
ibis::table::stringList Anames;
ibis::table::typeList Atypes;
ibis::table::stringList Acols;
std::map<const char*,ibis::TYPE_T> Anaty;
ibis::qExpr* Acond=0;

const char* Bfrom = 0;
const char* Bqcnd="1=1";
char* Bsel;
const char* Bstart="start";
const char* Bend = "end";
ibis::table::stringList Bcols;
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

const char* outdir;

ibis::whereClause senseWhere = ibis::whereClause("strand == 1");
ibis::qExpr* senseExpr = senseWhere.getExpr();

std::vector<ibis::table::bufferList*> part_results;

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
		<< "[-o output directory]" << std::endl
		<< "[-p parallelize using threads]" << std::endl;
}

// function to parse the command line arguments
void parse_args(int argc, char** argv)
{
	for (size_t i=1; i<argc; ++i) {
	    if (*argv[i] == '-') { // arguments starting with -
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
			  case 'o':
				  if (i+1 < argc)
					  outdir = argv[++i];
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
	    } // arguments starting with -
	} // for (i=1; ...)

/*
	std::cerr << *argv << " -d1 " << Afrom << " -d2 " << Bfrom
	  << " -s1 " << Astart << " -s2 " << Bstart
	  << " -e1 " << Aend << " -e2 " << Bend
	  << " -w1 " << Aqcnd << " -w2 " << Bqcnd
	  << " -c1 " << Asel << " -c2 " << Bsel
	  << " -r " << right << " -l " << left
	  << " -b " << bins << " -p " << parallelize
	  << " -sm " << same_strand << " -sw " << stranded_windows
	  << " -o " << outdir
	  << std::endl;
*/
}

// function to fill an array with values from an input array based on a list of indices
// Similar to ibis::util::reorder, but ind.size() can be > arr.size()
// because indices can be repeated.
template <typename T>
void meorder(ibis::array_t<T> &res, ibis::array_t<T> &arr, ibis::array_t<uint32_t>& ind)
{
	for (size_t i = 0; i < ind.size(); ++ i)
	    res[i] = arr[ind[i]];
}

// string version of meorder
void meorder(std::vector<std::string> &res, std::vector<std::string> &arr,
	ibis::array_t<uint32_t>& ind)
{
	for(size_t i = 0; i < ind.size(); ++ i)
		res[i] = arr[ind[i]];
}

// given a column, a bitvector mask, and an index
// fetch the relevant rows in the matching format
// and store in result
void fillColumn(ibis::column* col, ibis::bitvector* mask,
	ibis::array_t<uint32_t> &idx, void* result)
{
	switch(col->type()) {
		case ibis::BYTE:
		{
			ibis::array_t<signed char> *values = col->selectBytes(*mask);
			meorder(*static_cast<ibis::array_t<signed char>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::UBYTE:
		{
			ibis::array_t<unsigned char> *values = col->selectUBytes(*mask);
			meorder(*static_cast<ibis::array_t<unsigned char>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::SHORT:
		{
			ibis::array_t<int16_t> *values = col->selectShorts(*mask);
			meorder(*static_cast<ibis::array_t<int16_t>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::USHORT:
		{
			ibis::array_t<uint16_t> *values = col->selectUShorts(*mask);
			meorder(*static_cast<ibis::array_t<uint16_t>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::INT:
		{
			ibis::array_t<int32_t> *values = col->selectInts(*mask);
			meorder(*static_cast<ibis::array_t<int32_t>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::UINT:
		{
			ibis::array_t<uint32_t> *values = col->selectUInts(*mask);
			meorder(*static_cast<ibis::array_t<uint32_t>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::LONG:
		{
			ibis::array_t<int64_t> *values = col->selectLongs(*mask);
			meorder(*static_cast<ibis::array_t<int64_t>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::ULONG:
		{
			ibis::array_t<uint64_t> *values = col->selectULongs(*mask);
			meorder(*static_cast<ibis::array_t<uint64_t>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::FLOAT:
		{
			ibis::array_t<float> *values = col->selectFloats(*mask);
			meorder(*static_cast<ibis::array_t<float>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::DOUBLE:
		{
			ibis::array_t<double> *values = col->selectDoubles(*mask);
			meorder(*static_cast<ibis::array_t<double>*>(result),*values,idx);
			delete values;
			break;
		}
		case ibis::TEXT:
		case ibis::CATEGORY:
		{
			std::vector<std::string> *values = col->selectStrings(*mask);
			meorder(*static_cast<std::vector<std::string>*>(result),*values,idx);
			delete values;
			break;
		}
		default:
			break;
	}
}

// Build a table to represent the interval-join of two partitions
// using two bitvector masks and two arrays of offsets into the matching rows
// New columns for relative start and end are included (possibly binned)
// In the end the table needs to be sorted by these columns
// When Bstart == Bend you don't add an end column
void fillResult(const ibis::part* Apart, const ibis::part* Bpart,
	ibis::bitvector* Amatch, ibis::bitvector* Bmatch,
	ibis::array_t<uint32_t>& Aidx, ibis::array_t<uint32_t>& Bidx,
	ibis::array_t<int32_t> &relativeStart, ibis::array_t<int32_t> &relativeEnd, uint32_t thread_id)
{
	size_t nrows = Aidx.size();
	size_t ncols = Acols.size() + Bcols.size() + 2;
	ibis::table::bufferList * tbuff = new ibis::table::bufferList(ncols);
	part_results[thread_id]=tbuff;
	ibis::table::typeList ttypes(ncols);

	(*tbuff)[0] = &relativeStart;
	(*tbuff)[1] = &relativeEnd;

	boost::thread_group g1;
	size_t j=2;
	for(size_t i=0; i<Acols.size(); i++) {
		ibis::column *col = Apart->getColumn(Acols[i]);
		ttypes[j] = col->type();
		(*tbuff)[j] = ibis::table::allocateBuffer(col->type(),nrows);
		if(parallelize > 0)
			g1.create_thread(boost::bind(fillColumn,col,Amatch,Aidx,(*tbuff)[j]));
		else
			fillColumn(col,Amatch,Aidx,(*tbuff)[j]);
		j++;
	}
	for(size_t i=0; i<Bcols.size(); i++) {
		ibis::column *col = Bpart->getColumn(Bcols[i]);
		ttypes[j] = col->type();
		(*tbuff)[j] = ibis::table::allocateBuffer(col->type(),nrows);
		if(parallelize > 0)
			g1.create_thread(boost::bind(fillColumn,col,Bmatch,Bidx,(*tbuff)[j]));
		else
			fillColumn(col,Bmatch,Bidx,(*tbuff)[j]);
		j++;
	}
	if(parallelize > 0)
		g1.join_all();
}

// calculate the position of from and to relative to before_ and after_
// convert the relative start and end to bins if requested
// trim the ends if they don't fit in the binned window.
void push_relative_position(ibis::array_t<int32_t> &relativeStart,ibis::array_t<int32_t> &relativeEnd,
	bool flipped, uint32_t &before, uint32_t &after, uint32_t &from, uint32_t &to)
{
	int32_t startDiff, endDiff;
	if (flipped) {
		startDiff = after - to;
		endDiff = after - from;
	} else {
		startDiff = from - before;
		endDiff = to - before;
	}
	if (bins > 0) {
		int32_t winlength = after - before;
		startDiff = (startDiff <= 0) ? 0 : startDiff/winlength;
		endDiff = (endDiff >= winlength) ? winlength-1 : endDiff/winlength;
	}
	relativeStart.push_back(startDiff);
	relativeEnd.push_back(endDiff);
}

// fetch intervals from the partitions to be joined
// identify features that overlap
// set bits in Amatch and Bmatch so we can fetch only elements that overlap
// fill two arrays of indexes that say what rows to take from selected columns of A and B
void Stacker(const ibis::part* Apart, const ibis::part* Bpart,
	ibis::bitvector &Amask, ibis::bitvector &Bmask, int before, int after, bool flipped, uint32_t thread_id)
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


	// fill the relativeStart and relativeEnd arrays - normalized if requested
	// pass these to fillResult() so we can include them with the user-requested columns
	ibis::array_t<int32_t> * relativeStart = new ibis::array_t<int32_t>();
	ibis::array_t<int32_t> * relativeEnd = new ibis::array_t<int32_t>();
	
	ibis::bitvector* Amatch = new ibis::bitvector;
	ibis::bitvector* Bmatch = new ibis::bitvector;
	const uint32_t nA = Amask.cnt();
	const uint32_t nB = Bmask.cnt();
	ibis::array_t<uint32_t> Aidx;
	ibis::array_t<uint32_t> Bidx;
	uint32_t iA=0;
	uint32_t iB=0;
	uint32_t mA=0;
	uint32_t mB=0;

	std::map<uint32_t,bool> Bseen;
	while (iA < nA && iB < nB) {
		uint32_t before_ = (*Astart_val)[iA] - before;
		uint32_t after_ = (*Aend_val)[iA] + after;
		if ((*Bend_val)[iB] < before_) {// B comes first
			if(Bseen.count(iB)>0)
				mB++;
			iB++;
		}
		else if (after_ < (*Bstart_val)[iB]) // A comes first
			iA++;
		else {
			// found a match
			Aidx.push_back(mA);
			Bidx.push_back(mB);
			push_relative_position(*relativeStart,*relativeEnd,flipped,before_,after_,(*Bstart_val)[iB],(*Bend_val)[iB]);

			Amatch->setBit((*Arids)[iA].value, 1);
			if (Bseen.count(iB) == 0) {
				Bseen.insert(std::pair<uint32_t,bool>(iB,true));
				Bmatch->setBit((*Brids)[iB].value, 1);
			}
			// report all the matches of B within this A interval
			// then increment iA
			uint32_t iiB=iB+1;
			uint32_t mmB=mB+1;
			while (iiB < nB && (*Bstart_val)[iiB] <= after_) {
				if ((*Bend_val)[iiB] >= before_) {
					// overlap
					Aidx.push_back(mA);
					Bidx.push_back(mmB);
					push_relative_position(*relativeStart,*relativeEnd,flipped,before_,after_,(*Bstart_val)[iiB],(*Bend_val)[iiB]);
					if (Bseen.count(iiB) == 0) {
						Bseen.insert(std::pair<uint32_t,bool>(iiB,true));
						Bmatch->setBit((*Brids)[iiB].value, 1);
					}
					mmB++;
				}
				iiB++;
			}
			mA++;
			iA++;
		}
	}

	Amatch->adjustSize(0,Apart->nRows());
	Bmatch->adjustSize(0,Bpart->nRows());

	Amatch->compress();
	Bmatch->compress();
	if (Amatch->cnt() == 0 || Bmatch->cnt() == 0)
		return;

	fillResult(Apart,Bpart,Amatch,Bmatch,Aidx,Bidx,*relativeStart,*relativeEnd, thread_id);
	delete Amatch;
	delete Bmatch;
	delete Arids;
	delete Brids;
	delete Astart_val;
	delete Aend_val;
	delete Bstart_val;
	delete Bend_val;
}

// apply the user supplied query expression to the given part
// populate a bitvector mask that corresponds to the matching rows
// return the number of rows that match
int countHits(const ibis::part* part, ibis::bitvector &mask,
	 const ibis::qExpr* cond, const char* colName)
{
	if (cond != 0) {
		ibis::countQuery que(part);
		int ierr = que.setWhereClause(cond);
		ierr = que.evaluate();
		mask.copy(*que.getHitVector());
	}
	else {
		ibis::column *c = part->getColumn(colName);
		c->getNullMask(mask);
	}
	return mask.cnt();
}

// populate bitvectors for hits to the sense and antisense strands
// Behold the beauty of WAH bitwise operations!
void splitByStrand(const ibis::part* part, ibis::bitvector &mask,
	ibis::bitvector &plus, ibis::bitvector &minus)
{
	ibis::countQuery que(part);
	int ierr = que.setWhereClause(senseExpr);
	ierr = que.evaluate();
	plus.copy(*que.getHitVector());
	minus.copy(*que.getHitVector());
	minus.flip();
	plus &= mask;
	minus &= mask;
}

// intersecting genomic intervals are always on the same chromsome
// so we operate on partitions with the same FBchr
// other command line options control how intervals are selected
// so we have to handle special cases related to strand
void setupStacker(const ibis::part* Apart, const ibis::part* Bpart, size_t thread_id, uint32_t offset)
{
	// TODO: do some checks
	// count hits to A and B on this chr
	// make masks based on Aqcnd and Bqcnd
    ibis::bitvector Amask;
	int Arows = countHits(Apart, Amask, Acond, Astart);	
	if (Arows == 0) return;
	
    ibis::bitvector Bmask;
	int Brows = countHits(Bpart, Bmask, Bcond, Bstart);
	if (Brows == 0) return;
	
	if(same_strand == 1 && (Anaty.count("strand") > 0) && (Bnaty.count("strand") > 0)) {
		// split into subproblems for each strand
		// update the masks (Aplus, Aminus, Bplus, Bminus)
		ibis::bitvector Aplus;
		ibis::bitvector Aminus;
		splitByStrand(Apart, Amask, Aplus, Aminus);
		
		ibis::bitvector Bplus;
		ibis::bitvector Bminus;
		splitByStrand(Bpart, Bmask, Bplus, Bminus);

		Stacker(Apart, Bpart, Aplus, Bplus, left, right, false, thread_id);
		if (stranded_windows>0) {
			Stacker(Apart, Bpart, Aminus, Bminus, right, left, true, thread_id+offset);
		}
		else {
			Stacker(Apart, Bpart, Aminus, Bminus, left, right, false, thread_id+offset);
		}
	}
	else {
		if (stranded_windows>0 && left != right and Anaty.count("strand")>0) {
			// split into subproblems
			// update the masks (Aplus, Aminus)
			ibis::bitvector Aplus;
			ibis::bitvector Aminus;
			splitByStrand(Apart, Amask, Aplus, Aminus);
			Stacker(Apart, Bpart, Aplus, Bmask, left, right, false, thread_id);
			Stacker(Apart, Bpart, Aminus, Bmask, right, left, true, thread_id+offset);
		}
		else {
			Stacker(Apart, Bpart, Amask, Bmask, left, right, false, thread_id);
		}
	}
}

// make a list of valid column names from the given select string
void fillColumnLists(char * sel, std::map<const char*,ibis::TYPE_T> &naty,
	ibis::table::stringList &cols)
{
	cols.clear();
	char * pch;
	pch = strtok (sel, " ,.-");
	while (pch != NULL) {
		if (naty.count(pch) > 0)
			cols.push_back(pch);
		pch = strtok (NULL, " ,.-");
	}
}

// iterate over the partial results and stitch the arrays/vectors together for the
// given column number
void concatenate_column(ibis::table::bufferList& tbuff, ibis::table::typeList &ttypes, size_t col)
{
	tbuff[col] = ibis::table::allocateBuffer(ttypes[col],0);
	for(size_t i=0; i < part_results.size(); i++) {
		if (part_results[i] != 0) {
			switch(ttypes[col]) {
				case ibis::BYTE:
				{
					ibis::array_t<signed char>* b1 = static_cast<ibis::array_t<signed char>*>(tbuff[col]);
					ibis::array_t<signed char>* b2 = static_cast<ibis::array_t<signed char>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::UBYTE:
				{
					ibis::array_t<unsigned char>* b1 = static_cast<ibis::array_t<unsigned char>*>(tbuff[col]);
					ibis::array_t<unsigned char>* b2 = static_cast<ibis::array_t<unsigned char>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::SHORT:
				{
					ibis::array_t<int16_t>* b1 = static_cast<ibis::array_t<int16_t>*>(tbuff[col]);
					ibis::array_t<int16_t>* b2 = static_cast<ibis::array_t<int16_t>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::USHORT:
				{
					ibis::array_t<uint16_t>* b1 = static_cast<ibis::array_t<uint16_t>*>(tbuff[col]);
					ibis::array_t<uint16_t>* b2 = static_cast<ibis::array_t<uint16_t>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::INT:
				{
					ibis::array_t<int32_t>* b1 = static_cast<ibis::array_t<int32_t>*>(tbuff[col]);
					ibis::array_t<int32_t>* b2 = static_cast<ibis::array_t<int32_t>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::UINT:
				{
					ibis::array_t<uint32_t>* b1 = static_cast<ibis::array_t<uint32_t>*>(tbuff[col]);
					ibis::array_t<uint32_t>* b2 = static_cast<ibis::array_t<uint32_t>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::LONG:
				{
					ibis::array_t<int64_t>* b1 = static_cast<ibis::array_t<int64_t>*>(tbuff[col]);
					ibis::array_t<int64_t>* b2 = static_cast<ibis::array_t<int64_t>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::ULONG:
				{
					ibis::array_t<uint64_t>* b1 = static_cast<ibis::array_t<uint64_t>*>(tbuff[col]);
					ibis::array_t<uint64_t>* b2 = static_cast<ibis::array_t<uint64_t>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::FLOAT:
				{
					ibis::array_t<float>* b1 = static_cast<ibis::array_t<float>*>(tbuff[col]);
					ibis::array_t<float>* b2 = static_cast<ibis::array_t<float>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::DOUBLE:
				{
					ibis::array_t<double>* b1 = static_cast<ibis::array_t<double>*>(tbuff[col]);
					ibis::array_t<double>* b2 = static_cast<ibis::array_t<double>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				case ibis::TEXT:
				case ibis::CATEGORY:
				{
					std::vector<std::string>* b1 = static_cast<std::vector<std::string>*>(tbuff[col]);
					std::vector<std::string>* b2 = static_cast<std::vector<std::string>*>((*part_results[i])[col]);
					b1->insert(b1->end(),b2->begin(),b2->end());
					break;
				}
				default:
					break;
			}
		}
	}	
}

// iterate over the columns and concatenate the partial results for each
// a logical table is constructed in the end
ibis::table* concatenate_results()
{
	size_t ncols = Acols.size() + Bcols.size() + 2;
	ibis::table::bufferList tbuff(ncols,0);
	ibis::table::typeList ttypes(ncols);
	ibis::table::stringList colnames(ncols);
	std::vector<std::string> colnamesstr;
	colnamesstr.resize(ncols);
	IBIS_BLOCK_GUARD(ibis::table::freeBuffers, ibis::util::ref(tbuff), ibis::util::ref(ttypes));

	uint32_t nrows=0;
	for(size_t i=0; i < part_results.size(); i++) {
		if (part_results[i] != 0) {
			ibis::array_t<int32_t>* buff = static_cast<ibis::array_t<int32_t>*>((*part_results[i])[0]);
			nrows += buff->size();
		}
	}
	
	colnames[0] = "start";
	colnames[1] = "end";
	ttypes[0] = ibis::INT;
	ttypes[1] = ibis::INT;
	boost::thread_group g2;
	if (parallelize > 0) {
		g2.create_thread(boost::bind(concatenate_column, tbuff, ttypes, 0));
		g2.create_thread(boost::bind(concatenate_column, tbuff, ttypes, 1));
	} else {
		concatenate_column(tbuff,ttypes,0);
		concatenate_column(tbuff,ttypes,1);
	}
	size_t j=2;
	for(size_t i=0; i<Acols.size(); i++) {
		std::ostringstream oss;
		oss << "A__" << Acols[i];
		colnamesstr[j] = oss.str();
		colnames[j] = colnamesstr[j].c_str();
		ttypes[j] = Anaty.find(Acols[i])->second;
		if (parallelize > 0)
			g2.create_thread(boost::bind(concatenate_column, tbuff, ttypes,j));
		else
			concatenate_column(tbuff,ttypes,j);
		j++;
	}
	for(size_t i=0; i<Bcols.size(); i++) {
		std::ostringstream oss;
		oss << "B__" << Bcols[i];
		colnamesstr[j] = oss.str();
		colnames[j] = colnamesstr[j].c_str();
		ttypes[j] = Bnaty.find(Bcols[i])->second;
		if (parallelize > 0)
			g2.create_thread(boost::bind(concatenate_column, tbuff, ttypes,j));
		else
			concatenate_column(tbuff,ttypes,j);
		j++;
	}
	if (parallelize > 0)
		g2.join_all();
	
	// create the table now that it's filled
    return new ibis::bord("joined", "joined tables", nrows, tbuff, ttypes, colnames,&colnames,0);
}

// fetch an array of values from the table and write them to outdir/cname
// ibis::table::getColumnAs...() makes a copy so values are contiguous
// in memory.  That makes it possible to UnixWrite() the column in bulk.
// -- not sure if categorical columns are being done the right way.
// -- relying on subsequent fastbit activity to index/process as needed.
void write_column(ibis::table* tbl, const char* cname, ibis::TYPE_T ctype, size_t nr, const char* outdir)
{
	// open the output file
	std::string cnm = outdir;
	cnm += FASTBIT_DIRSEP;
	cnm += cname;
	int fdes = UnixOpen(cnm.c_str(), OPEN_WRITEADD, OPEN_FILEMODE);
	if (fdes < 0) {
		std::cerr << "faileed to open file " << cnm << " for writing" << std::endl;
		exit(-4);
	}
	IBIS_BLOCK_GUARD(UnixClose, fdes);
#if defined(_WIN32) && defined(_MSC_VER)
	(void)_setmode(fdes, _O_BINARY);
#endif
	uint64_t ierr;
	off_t pos = UnixSeek(fdes,0,SEEK_END);
	// fetch results and write
	switch (ctype) {
		case ibis::BYTE:
		{
			char* vals = new char[nr];			
			ierr = tbl->getColumnAsBytes(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(char)*nr);
			delete vals;
			break;
		}
		case ibis::UBYTE:
		{
			unsigned char* vals = new unsigned char[nr];
			ierr = tbl->getColumnAsUBytes(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(unsigned char)*nr);
			delete vals;
			break;
		}
		case ibis::SHORT:
		{
		    int16_t* vals = new int16_t[nr];
			ierr = tbl->getColumnAsShorts(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(int16_t)*nr);
			delete vals;
			break;
		}
		case ibis::USHORT:
		{
		    uint16_t* vals = new uint16_t[nr];
			ierr = tbl->getColumnAsUShorts(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(uint16_t)*nr);
			delete vals;
			break;
		}
		case ibis::INT:
		{
		    int32_t* vals = new int32_t[nr];
			ierr = tbl->getColumnAsInts(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(int32_t)*nr);
			delete vals;
			break;
		}
		case ibis::UINT:
		{
		    uint32_t* vals = new uint32_t[nr];
			ierr = tbl->getColumnAsUInts(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(uint32_t)*nr);
			delete vals;
			break;
		}
		case ibis::LONG:
		{
		    int64_t* vals = new int64_t[nr];
			ierr = tbl->getColumnAsLongs(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(int64_t)*nr);
			delete vals;
			break;
		}
		case ibis::ULONG:
		{
		    uint64_t* vals = new uint64_t[nr];
			ierr = tbl->getColumnAsULongs(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(uint64_t)*nr);
			delete vals;
			break;
		}
		case ibis::FLOAT:
		{
			float* vals = new float[nr];
			ierr = tbl->getColumnAsFloats(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(float)*nr);
			delete vals;
			break;
		}
		case ibis::DOUBLE:
		{
			double* vals = new double[nr];
			ierr = tbl->getColumnAsDoubles(cname,vals);
			pos = UnixWrite(fdes,vals,sizeof(double)*nr);
			delete vals;
			break;
		}
		case ibis::CATEGORY:
		case ibis::TEXT:
		{
		    std::vector<std::string>* vals = new std::vector<std::string>();
		    ierr = tbl->getColumnAsStrings(cname, *vals);
			for (uint32_t j = 0; j < nr; ++ j)
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

// dumps the contents of an in-memory table to a new partition in outdir
// creates metadata file outdir/-part.txt
// and writes each column to outdir/cname using multiple threads
// subsequent fastbit activity will build any necessary indexes.
void write_table(ibis::table* tbl, const char* outdir)
{
	const char* startcolname="start";
	ibis::table::stringList names = tbl->columnNames();
	ibis::table::typeList types = tbl->columnTypes();
	size_t nr = tbl->nRows();

	ibis::part tmp(outdir, static_cast<const char*>(0));

	std::string mdfile = outdir;
	mdfile += FASTBIT_DIRSEP;
	mdfile += "-part.txt";
	std::ofstream md(mdfile.c_str());
	if(! md) {
		std::cerr << "failed to open metadata file " << mdfile << std::endl;
		exit(-3);
	}

	md << "BEGIN HEADER\nName = joined" << "\nDescription = joined"
	   << "\nNumber_of_rows = " << nr
	   << "\nNumber_of_columns = " << names.size();
	md << "\nEND HEADER\n";

	boost::thread_group tg;
	for(size_t i=0; i<names.size();i++) {
		if(parallelize>0)
			tg.create_thread(boost::bind(write_column, tbl, names[i], types[i], nr, outdir));
		else
			write_column(tbl,names[i],types[i],nr,outdir);
		md << "\nBegin Column\nname = " << names[i] << "\ndata_type = "
		   << ibis::TYPESTRING[(int) types[i]];
	    if(stricmp(names[i], startcolname)==0)
			md << "\nsorted = true";
		md << "\nEnd Column\n";
	}
	md.close();
	if(parallelize>0)
		tg.join_all();
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
	Acols = A->columnNames();
    for (size_t i = 0; i < Anames.size(); ++ i) {
		Anaty.insert(std::pair<const char*,ibis::TYPE_T>(Anames[i],Atypes[i]));
	}
	ibis::table* B = ibis::table::create(Bfrom);
	Bnames = B->columnNames();
	Btypes = B->columnTypes();
	Bcols = B->columnNames();
    for (size_t i = 0; i < Bnames.size(); ++ i) {
		Bnaty.insert(std::pair<const char*,ibis::TYPE_T>(Bnames[i],Btypes[i]));
	}

	// parse Asel and Bsel columns
	if (Asel != 0 && *Asel != 0)
		fillColumnLists(Asel,Anaty,Acols);
	if (Bsel != 0 && *Bsel != 0)
		fillColumnLists(Bsel,Bnaty,Bcols);

	std::vector<const ibis::part*> Bparts;
	B->getPartitions(Bparts);

	// build a map of B's partitions so we can only join on matching FBchrs
	std::map<const char*,const ibis::part*> Bpartmap;	
	for(size_t i=0; i<Bparts.size(); ++ i)
		Bpartmap.insert(std::pair<const char*,const ibis::part*>(Bparts[i]->getMetaTag("FBchr"),Bparts[i]));
	
	std::vector<const ibis::part*> Aparts;
	A->getPartitions(Aparts);
	part_results.resize(Aparts.size()*2); // *2 because of strand
	boost::thread_group g0;
	for(size_t i=0; i<Aparts.size(); ++ i) {
		const char* chr = Aparts[i]->getMetaTag("FBchr");
		if (Bpartmap.count(chr)>0) {
			if (parallelize > 0)
				g0.create_thread(boost::bind(setupStacker, Aparts[i], Bpartmap.find(chr)->second,i,Aparts.size()));
			else
				setupStacker(Aparts[i],Bpartmap.find(chr)->second,i,Aparts.size());
		}
	}
	if (parallelize > 0)
		g0.join_all();

	// concatenate each part into one table, sort, and save
	ibis::table *res = concatenate_results();
	res->orderby("start, end");
	write_table(res,outdir);
	
	delete res;
	delete A;
	delete B;
    return 0;
}
