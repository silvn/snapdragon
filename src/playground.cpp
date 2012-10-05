#include <unistd.h>
#include "ibis.h"
#include "bitvector64.h"

using namespace std;

// http://en.wikipedia.org/wiki/Phi_coefficient
double bitvector_phi(ibis::bitvector64 &bv1, ibis::bitvector64 &bv2)
{
	if (bv1.size() < bv2.size()) {
		bv1.adjustSize(0,bv2.size());
	}
	else if (bv1.size() > bv2.size()) {
		bv2.adjustSize(0,bv1.size());
	}
	ibis::bitvector64 bvboth = ibis::bitvector64(bv1);
	bvboth &= bv2;
	double n11 = bvboth.cnt();
	double n1x = bv1.cnt();
	double nx1 = bv2.cnt();
	double n0x = bv1.size() - n1x;
	double nx0 = bv1.size() - nx1;
	double n10 = n1x - n11;
	double n01 = nx1 - n11;
	double n00 = n0x - n01;
	return (n11 * n00 - n10 * n01)/sqrt(n1x * n0x * nx0 * nx1);	
}

double bitvector_phi2(ibis::bitvector64 &bv1, ibis::bitvector64 &bv2)
{
	if (bv1.size() < bv2.size()) {
		bv1.adjustSize(0,bv2.size());
	}
	else if (bv1.size() > bv2.size()) {
		bv2.adjustSize(0,bv1.size());
	}
	ibis::bitvector64 bvu(bv1);
	bvu |= bv2;
	double n1x = bv1.cnt();
	double nx1 = bv2.cnt();
	double n0x = bvu.cnt() - n1x;
	double nx0 = bvu.cnt() - nx1;
	ibis::bitvector64 bvboth = ibis::bitvector64(bv1);
	bvboth &= bv2;
	double n11 = bvboth.cnt();
	double n10 = n1x - n11;
	double n01 = nx1 - n11;
	double n00 = n0x - n01; // 0 by definition, so using n11*n11 instead of n11*n00
	return (n11 * n11 - n10 * n01)/sqrt(n1x * n0x * nx0 * nx1);	
}

int dump_intervals(ibis::bitvector64 &bv)
{
	ibis::bitvector64::word_t n=0;
	ibis::bitvector64::indexSet index = bv.firstIndexSet();
	ibis::bitvector64::word_t from=0;
	ibis::bitvector64::word_t to=0;
    while (index.nIndices() > 0) {
		const ibis::bitvector64::word_t *idx0 = index.indices();
		if (index.isRange()) {
			if (*idx0 == to+1)
				to+= index.nIndices();
			else {
				if (n>0)
					cout << from << "\t" << to << endl;;
				from = *idx0;
				to = *idx0 + index.nIndices() - 1;
				n++;
			}
		}
		else {
		    for (ibis::bitvector64::word_t j = 0; j<index.nIndices(); ++j) {
				if (idx0[j] == to+1)
					to = idx0[j];
				else {
					if (n>0)
						cout << from << "\t" << to << endl;;
					from=idx0[j];
					to=idx0[j];
					n++;
				}
		    }
		}
		++ index;
    }
	if (n>0)
		cout << from << "\t" << to << endl;;
	return n;
}

int main(int argc, char *argv[]) {

	// generate some bitvectors
	ibis::bitvector64 bv1,bv2,bv3;
	ibis::bitvector64::word_t intervals = 100000;
	ibis::bitvector64::word_t max_gap = 1000000;
	ibis::bitvector64::word_t min_gap = 100;
	ibis::bitvector64::word_t max_run = 10000;
	ibis::bitvector64::word_t min_run = 10;
	ibis::bitvector64::word_t g,r;
	srand(time(NULL));
	for(ibis::bitvector64::word_t i=0; i < intervals; i++) {
		g = rand() % max_gap + min_gap;
		r = rand() % max_run + min_run;
		bv1.appendFill(0,g);
		bv1.appendFill(1,r);
		g = rand() % max_gap + min_gap;
		r = rand() % max_run + min_run;
		bv2.appendFill(0,g);
		bv2.appendFill(1,r);
	}
	bv3.appendFill(1,31);
	bv3.appendFill(0,3);
	bv3.appendFill(1,77);

	cout << "bv1 has " << bv1.size() << " bits, " << bv1.cnt() << " set bits, bytes= " << bv1.bytes() << endl;
	cout << "bv2 has " << bv2.size() << " bits, " << bv2.cnt() << " set bits, bytes= " << bv2.bytes() << endl;
	cout << "bv3 has " << bv3.size() << " bits, " << bv3.cnt() << " set bits, bytes= " << bv3.bytes() << endl;

	// calculate phi coefficient between bv1 and bv2
	double phi = bitvector_phi(bv1,bv2);
	cout << "phi(bv1,bv2)= " << phi << endl;

	// AND two bitvectors
	ibis::bitvector64 bv12(bv1);
	bv12 &= bv2;
	
	cout << "bv1 & bv2 has " << bv12.cnt() << " set bits, bytes= " << bv12.bytes() << endl;

	// output the positions of the set bits
	int count = dump_intervals(bv3);
	cout << "dumped " << count << " intervals" << endl;

	// generate some random (sparse) bitvectors
	ibis::bitvector64::word_t nvec = 1000;
	std::vector<ibis::bitvector64> bvectors;
    ibis::horometer timer;
    timer.start();
	intervals = 1000;
	max_run = 1;
	min_run = 1;
	max_gap = 10;
	min_gap = 1;
	ibis::bitvector64::word_t max=0;
	for(ibis::bitvector64::word_t i=0; i<nvec; i++) {
		ibis::bitvector64 bvi;
		for(ibis::bitvector64::word_t j=0; j < intervals; j++) {
			g = rand() % max_gap + min_gap;
			r = rand() % max_run + min_run;
			bvi.appendFill(0,g);
			bvi.appendFill(1,r);
		}
		bvectors.push_back(bvi);
		if (bvi.size() > max)
			max = bvi.size();
	}
	timer.stop();
	cout << "max= " << max << endl;
	cout << "generating " << nvec << " bitvectors took "
	     << timer.CPUTime() << " CPU seconds, "
	     << timer.realTime() << " elapsed seconds" << endl;
	
	// adjust their lengths to max
	timer.start();
	for(int i=0; i<nvec; i++) {
		bvectors[i].adjustSize(0,max);
	}
	timer.stop();
	cout << "resizing " << nvec << " bitvectors took "
	     << timer.CPUTime() << " CPU seconds, "
	     << timer.realTime() << " elapsed seconds" << endl;
	
	
	// OR them all together
	timer.start();
	ibis::bitvector64 bigOR;
	bigOR.adjustSize(0,max);
	for(int i=0; i<nvec; i++) {
		bigOR |= bvectors[i];
	}
	timer.stop();
	cout << "ORing " << nvec << " bitvectors took "
	     << timer.CPUTime() << " CPU seconds, "
	     << timer.realTime() << " elapsed seconds" << endl;

	// calculate phi between all pairs of vectors
	timer.start();
	double max_phi=0;
	ibis::bitvector64::word_t max_phii=0;
	ibis::bitvector64::word_t max_phij=0;
	for(ibis::bitvector64::word_t i=0; i<nvec-1; i++) {
		for(ibis::bitvector64::word_t j=i+1; j<nvec; j++) {
			double phi_ij = bitvector_phi(bvectors[i],bvectors[j]);
			if (phi_ij > max_phi) {
				max_phi = phi_ij;
				max_phii = i;
				max_phij = j;
			}
		}
	}
	timer.stop();
	cout << "max_phi[" << max_phii << "," << max_phij << "]=" << max_phi << endl;
	cout << "calculating max_phi took "
		<< timer.CPUTime() << " CPU seconds, "
		<< timer.realTime() << " elapsed seconds" << endl;

}
