#include <node.h>
#include <ibis.h>
#include <category.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// constances defined in bitvector
const unsigned MAXBITS = 8*sizeof(ibis::bitvector::word_t) - 1;
const unsigned SECONDBIT = MAXBITS - 1;
const ibis::bitvector::word_t ALLONES = ((1U << MAXBITS) - 1);
const ibis::bitvector::word_t MAXCNT = ((1U << SECONDBIT) - 1);
const ibis::bitvector::word_t FILLBIT = (1U << SECONDBIT);
const ibis::bitvector::word_t HEADER0 = (2U << SECONDBIT);
const ibis::bitvector::word_t HEADER1 = (3U << SECONDBIT);

std::string c_stringify(v8::Local<v8::Value> lv) {
	v8::String::Utf8Value p(lv);
	return std::string(*p,p.length());
}

v8::Handle<v8::Object> TabletoJs(const ibis::table &tbl)
{
	size_t nr = tbl.nRows();
	ibis::table::stringList nms = tbl.columnNames();
	ibis::table::typeList tps = tbl.columnTypes();

	v8::Handle<v8::Object> results = v8::Object::New();

	std::vector<const ibis::part*> parts;
	tbl.getPartitions(parts);
	if (parts.size() != 1) {
		std::cerr << "parts.size() = " << parts.size() << " expected 1" << std::endl;
		return results;
	}
	int ierr = 0;
	for (size_t j = 0; j < nms.size(); ++ j) {
		v8::Handle<v8::Array> column = v8::Array::New();
		switch (tps[j]) {
			case ibis::BYTE:
			case ibis::SHORT:
			case ibis::INT:
			{
				int32_t *buf = new int32_t[nr];
				ierr = tbl.getColumnAsInts(nms[j], buf);
				for (size_t i = 0; i < nr; ++ i)
					column->Set(i,v8::Int32::New(buf[i]));
				break;
			}
			case ibis::UBYTE:
			case ibis::USHORT:
			case ibis::UINT:
			{
				uint32_t *buf = new uint32_t[nr];
				ierr = tbl.getColumnAsUInts(nms[j], buf);
				ibis::column *col = parts[0]->getColumn(nms[j]);
				const ibis::dictionary *dic = col->getDictionary();
				if (dic == 0)
					for (size_t i = 0; i < nr; ++ i)
						column->Set(i,v8::Uint32::New(buf[i]));
				else
					for (size_t i = 0; i < nr; ++ i)
						column->Set(i,v8::String::New((*dic)[buf[i]]));
				break;
			}
			case ibis::FLOAT:
			case ibis::DOUBLE:
			{
				double *buf = new double[nr];
				ierr = tbl.getColumnAsDoubles(nms[j], buf);
				for (size_t i = 0; i < nr; ++ i)
					column->Set(i,v8::Number::New(buf[i]));
				break;
			}
			case ibis::TEXT:
			case ibis::CATEGORY:
			{
			    std::vector<std::string>* buf = new std::vector<std::string>();
				ierr = tbl.getColumnAsStrings(nms[j], *buf);
				for (size_t i = 0; i < nr; ++ i)
					column->Set(i,v8::String::New((*buf)[i].c_str()));
				break;
			}
			default:
			{
				break;
			}
		}
		results->Set(v8::String::New(nms[j]),column);
	}
	return results;
}

v8::Handle<v8::Value> run_query(v8::Handle<v8::Object> p, ibis::table*& res)
{
	if (! p->Has(v8::String::New("from"))) {
		v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'from' argument")));
		return v8::Undefined();
	}
	if (! p->Has(v8::String::New("select"))) {
		v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'select' argument")));
		return v8::Undefined();
	}

	std::string data_dir = c_stringify(p->Get(v8::String::New("from")));
	ibis::table* tbl = ibis::table::create(data_dir.c_str());
	
	const ibis::qExpr* query;
	if (p->Has(v8::String::New("where"))) {
		std::string query_cnd = c_stringify(p->Get(v8::String::New("where")));
		ibis::whereClause tmp = ibis::whereClause(query_cnd.c_str());
		query = tmp.getExpr()->dup();
	}
	else {
		ibis::whereClause tmp = ibis::whereClause("1=1");
		query = tmp.getExpr()->dup();
	}
	
	// check for qExpr
	
	std::string select_str = c_stringify(p->Get(v8::String::New("select")));
	res = tbl->select(select_str.c_str(), query);
	
	if (p->Has(v8::String::New("orderby"))) {
		std::string order_by = c_stringify(p->Get(v8::String::New("orderby")));
		res->orderby(order_by.c_str());
	}
	return v8::String::New("OK");
}

// histogram params
bool adaptive = false;
uint32_t nbins=25;
double begin,end,stride;
v8::Handle<v8::Value> parse_histogram_params(v8::Handle<v8::Object> p)
{
	if (p->Has(v8::String::New("adaptive")))
		if (p->Get(v8::String::New("adaptive"))->IsTrue()) {
			adaptive = true;
			if (p->Has(v8::String::New("nbins")))
				nbins = p->Get(v8::String::New("nbins"))->Uint32Value();
		}

	if (! adaptive) {
		if (!p->Has(v8::String::New("begin"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'begin' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("end"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'end' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("stride"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'stride' argument")));
			return v8::Undefined();
		}
		begin = p->Get(v8::String::New("begin"))->NumberValue();
		end = p->Get(v8::String::New("end"))->NumberValue();
		stride = p->Get(v8::String::New("stride"))->NumberValue();
	}
	return v8::String::New("OK");
}

// 1D histogram
// returns a Javascript object {bounds:[], counts:[]}
// This function uses ibis::part histogram functions
// so we preprocess them into a single in-memory table (with one partition).
//
// binning is either adaptive or uniform
v8::Handle<v8::Value> histogram(const v8::Arguments& args)
{

	// parse args
	v8::Handle<v8::Object> p = v8::Handle<v8::Object>::Cast(args[0]);
	v8::Handle<v8::Value> rc = parse_histogram_params(p);
	if (rc->IsUndefined())
		return rc;

	// preprocess
	ibis::table *res = 0;
	rc = run_query(p,res);
	if (rc->IsUndefined())
		return rc;

	// no data after querying
	if (res->nRows() <= 0)
		return v8::Undefined();

	ibis::table::stringList nms = res->columnNames();

	std::vector<const ibis::part*> parts;
	res->getPartitions(parts);
	if (parts.size() != 1) {
		v8::ThrowException(v8::Exception::Error(v8::String::New("WTF! expected one partition after preprocessing")));
		return v8::Undefined();
	}

	long ierr;
	v8::Handle<v8::Object> JSON = v8::Object::New();
	v8::Handle<v8::Array> v8bounds = v8::Array::New();
	v8::Handle<v8::Array> v8counts = v8::Array::New();
	if (adaptive) {
		std::vector<double> bounds;
		std::vector<uint32_t> counts;
		ierr = parts[0]->get1DDistribution("1=1",nms[0],nbins,bounds,counts);
		if (ierr < 0) {
			v8::ThrowException(v8::Exception::Error(v8::String::New("adaptive 1D Distribution error")));
			return v8::Undefined();
		}
		// format results
		for(size_t i=0; i < counts.size(); i++) {
			v8bounds->Set(i,v8::Number::New(bounds[i]));
			v8counts->Set(i,v8::Number::New(counts[i]));
		}
	}
	else {
		std::vector<uint32_t> counts;
		ierr = parts[0]->get1DDistribution("1=1",nms[0],begin,end,stride,counts);
		if (ierr < 0) {
			v8::ThrowException(v8::Exception::Error(v8::String::New("uniform 1D Distribution error")));
			return v8::Undefined();
		}
		// format results
		double pos = begin;
		for(size_t i=0; i < counts.size(); i++) {
			v8bounds->Set(i,v8::Number::New(pos));
			v8counts->Set(i,v8::Number::New(counts[i]));
			pos += stride;
		}
	}
	JSON->Set(v8::String::New("bounds"),v8bounds);
	JSON->Set(v8::String::New("counts"),v8counts);
	return JSON;
}

// 2D histogram params
uint32_t nbins1=25;
uint32_t nbins2=25;
double begin1,end1,stride1,begin2,end2,stride2;
v8::Handle<v8::Value> parse_scatter_params(v8::Handle<v8::Object> p)
{
	if (p->Has(v8::String::New("adaptive")))
		if (p->Get(v8::String::New("adaptive"))->IsTrue()) {
			adaptive = true;
			if (p->Has(v8::String::New("nbins1")))
				nbins1 = p->Get(v8::String::New("nbins1"))->Uint32Value();
			if (p->Has(v8::String::New("nbins2")))
				nbins2 = p->Get(v8::String::New("nbins2"))->Uint32Value();
		}

	if (! adaptive) {
		if (!p->Has(v8::String::New("begin1"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'begin1' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("end1"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'end1' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("stride1"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'stride1' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("begin2"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'begin2' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("end2"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'end2' argument")));
			return v8::Undefined();
		}
		if (!p->Has(v8::String::New("stride2"))) {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("Missing 'stride2' argument")));
			return v8::Undefined();
		}
		begin1 = p->Get(v8::String::New("begin1"))->NumberValue();
		end1 = p->Get(v8::String::New("end1"))->NumberValue();
		stride1 = p->Get(v8::String::New("stride1"))->NumberValue();
		begin2 = p->Get(v8::String::New("begin2"))->NumberValue();
		end2 = p->Get(v8::String::New("end2"))->NumberValue();
		stride2 = p->Get(v8::String::New("stride2"))->NumberValue();
	}
	return v8::String::New("OK");
}

// 2D histogram
// returns a Javascript object {bounds1:[], bounds2:[], counts:[]}
// This function uses ibis::part histogram functions
// so we preprocess them into a single in-memory table (with one partition).
//
// binning is either adaptive or uniform
v8::Handle<v8::Value> scatter(const v8::Arguments& args)
{
	// parse args
	v8::Handle<v8::Object> p = v8::Handle<v8::Object>::Cast(args[0]);
	v8::Handle<v8::Value> rc = parse_scatter_params(p);
	if (rc->IsUndefined())
		return rc;

	// preprocess
	ibis::table *res=0;
	rc = run_query(p,res);
	if (rc->IsUndefined())
		return rc;

	// no data after querying
	if (res->nRows() <= 0)
		return v8::Undefined();

	ibis::table::stringList nms = res->columnNames();

	std::vector<const ibis::part*> parts;
	res->getPartitions(parts);
	if (parts.size() != 1) {
		v8::ThrowException(v8::Exception::Error(v8::String::New("WTF! expected one partition after preprocessing")));
		return v8::Undefined();
	}

	long ierr;
	v8::Handle<v8::Object> JSON = v8::Object::New();
	v8::Handle<v8::Array> v8bounds1 = v8::Array::New();
	v8::Handle<v8::Array> v8bounds2 = v8::Array::New();
	v8::Handle<v8::Array> v8counts = v8::Array::New();
	std::vector<double> bounds1;
	std::vector<double> bounds2;
	std::vector<uint32_t> counts;
	if (adaptive) {
		ierr = parts[0]->get2DDistribution(nms[0],nms[1],nbins1,nbins2,bounds1,bounds2,counts);
		if (ierr < 0) {
			v8::ThrowException(v8::Exception::Error(v8::String::New("adaptive 2D Distribution error")));
			return v8::Undefined();
		}
		// format results
		for(size_t i=0; i < counts.size(); i++)
			v8counts->Set(i,v8::Number::New(counts[i]));
		for(size_t i=0; i < bounds1.size(); i++)
			v8bounds1->Set(i,v8::Number::New(bounds1[i]));
		for(size_t i=0; i < bounds2.size(); i++)
			v8bounds2->Set(i,v8::Number::New(bounds2[i]));
	}
	else {
		std::vector<uint32_t> counts;
		ierr = parts[0]->get2DDistribution("1=1",nms[0],begin1,end1,stride1,nms[1],begin2,end2,stride2,counts);
		if (ierr < 0) {
			v8::ThrowException(v8::Exception::Error(v8::String::New("uniform 2D Distribution error")));
			return v8::Undefined();
		}
		// format results
		for(size_t i=0; i < counts.size(); i++)
			v8counts->Set(i,v8::Number::New(counts[i]));
		double pos = begin1;
		for(size_t i=0; i < bounds1.size(); i++) {
			v8bounds1->Set(i,v8::Number::New(pos));
			pos += stride1;
		}
		pos = begin2;
		for(size_t i=0; i < bounds2.size(); i++) {
			v8bounds2->Set(i,v8::Number::New(pos));
			pos += stride2;
		}
	}
	JSON->Set(v8::String::New("counts"),v8counts);
	JSON->Set(v8::String::New("bounds1"),v8bounds1);
	JSON->Set(v8::String::New("bounds2"),v8bounds2);
	return JSON;
}

// generic simple SQL function
v8::Handle<v8::Value> SQL(const v8::Arguments& args)
{
	v8::HandleScope scope;
	if(args.Length() != 1)
	{
		v8::ThrowException(v8::Exception::TypeError(v8::String::New("Wrong number of arguments")));
		return scope.Close(v8::Undefined());
	}

	v8::Handle<v8::Object> p = v8::Handle<v8::Object>::Cast(args[0]);

	ibis::table *res=0;
	v8::Handle<v8::Value> rc = run_query(p,res);

	if (rc->IsUndefined())
		return scope.Close(rc);

	v8::Handle<v8::Object> jsObj = TabletoJs(*res);
	delete res;
	return jsObj;
}

v8::Handle<v8::Object> BVencode(ibis::bitvector* bv)
{
	ibis::array_t<uint32_t> arr;
	bv->write(arr);
	v8::Handle<v8::Array> v8_arr = v8::Array::New(arr.size());
	for(size_t i=0; i < arr.size(); i++)
		v8_arr->Set(i,v8::Uint32::New(arr[i]));
	return v8_arr;
}

void BVdecode(v8::Handle<v8::Value> lv, ibis::bitvector* bv)
{
	v8::Handle<v8::Array> lva = v8::Handle<v8::Array>::Cast(lv);
	ibis::array_t<uint32_t> arr;
	arr.reserve(lva->Length());
	for(size_t i=0; i < lva->Length(); i++)
		arr.push_back(lva->Get(i)->Uint32Value());
	ibis::bitvector BitVector(arr);
	BitVector.swap(*bv);
}

v8::Handle<v8::Value> cnt(const v8::Arguments& args)
{
	ibis::bitvector bv;
	BVdecode(args[0],&bv);
	return v8::Integer::New(bv.cnt());
}

v8::Handle<v8::Value> size(const v8::Arguments& args)
{
	ibis::bitvector bv;
	BVdecode(args[0],&bv);
	return v8::Integer::New(bv.size());
}

v8::Handle<v8::Value> bytes(const v8::Arguments& args)
{
	ibis::bitvector bv;
	BVdecode(args[0],&bv);
	return v8::Integer::New(bv.bytes());
}

v8::Handle<v8::Value> set2bvec(const v8::Arguments& args)
{
	if (! args[0]->IsArray())
		return v8::String::New("first and only argument needs to be an array");
	v8::Handle<v8::Array> set = v8::Handle<v8::Array>::Cast(args[0]);
	ibis::bitvector bv;
	uint32_t next=0;
	for(size_t i=0; i < set->Length(); i++) {
		uint32_t v = set->Get(i)->Uint32Value();
		uint32_t g = v - next;
		if (g>0)
			bv.appendFill(0,g);
		bv.appendFill(1,1);
		next = v+1;
	}
	return BVencode(&bv);
}

v8::Handle<v8::Value> intervals2bvec(const v8::Arguments& args)
{
	if (! args[0]->IsArray())
		return v8::String::New("first argument needs to be an array");
	if (! args[1]->IsArray())
		return v8::String::New("second argument needs to be an array");
	v8::Handle<v8::Array> v8_starts = v8::Handle<v8::Array>::Cast(args[0]);
	v8::Handle<v8::Array> v8_ends = v8::Handle<v8::Array>::Cast(args[1]);
	if (v8_starts->Length() != v8_ends->Length())
		return v8::String::New("starts and ends arrays have to be the same length");

	// need to sort by start position
	ibis::array_t<uint32_t> starts;
	ibis::array_t<uint32_t> ends;
	for (uint32_t i=0; i<v8_starts->Length(); i++) {
		starts.push_back(v8_starts->Get(i)->Uint32Value());
		ends.push_back(v8_ends->Get(i)->Uint32Value());
	}
	ibis::array_t<uint32_t> ind;
	starts.sort(ind);
	
	ibis::bitvector bv;
	uint32_t next = 0;
	for(size_t i=0; i < ind.size(); i++) {
		uint32_t g = starts[ind[i]] - next;
		if (g>0)
			bv.appendFill(0,g);
		uint32_t f = ends[ind[i]] - starts[ind[i]];
		bv.appendFill(1,f);
		next = ends[ind[i]];
	}
	return BVencode(&bv);
}

v8::Handle<v8::Value> bvec2set(const v8::Arguments& args)
{
	ibis::bitvector bv;
	BVdecode(args[0],&bv);
	size_t i=0;
	v8::Handle<v8::Array> set = v8::Array::New();
	ibis::bitvector::indexSet index = bv.firstIndexSet();
	while (index.nIndices() > 0) {
		const ibis::bitvector::word_t *idx0 = index.indices();
		if (index.isRange()) {
			for(size_t j = 0; j < index.nIndices(); j++) {
				set->Set(i,v8::Number::New(*idx0+j));
				i++;
			}
		}
		else {
		    for (ibis::bitvector::word_t j = 0; j<index.nIndices(); ++j) {
				set->Set(i,v8::Number::New(idx0[j]));
				i++;
		    }
		}
		++ index;
	}

	return set;
}
// make it return an array of arrays instead
v8::Handle<v8::Value> bvec2span(const v8::Arguments& args)
{
	ibis::bitvector bv;
	BVdecode(args[0],&bv);
	size_t i=0;
	v8::Handle<v8::Array> set = v8::Array::New();
	ibis::bitvector::indexSet index = bv.firstIndexSet();
	while (index.nIndices() > 0) {
		const ibis::bitvector::word_t *idx0 = index.indices();
		if (index.isRange()) {
			set->Set(i,v8::Number::New(*idx0));
			i++;
			set->Set(i,v8::Number::New(*idx0+index.nIndices()-1));
			i++;
		}
		else {
			for (ibis::bitvector::word_t j = 0; j<index.nIndices(); ++j) {
				set->Set(i,v8::Number::New(idx0[j]));
				i++;
				set->Set(i,v8::Number::New(idx0[j]));
				i++;
			}
		}
		++ index;
	}

	return set;
}

// logical operations
v8::Handle<v8::Value> logical(const v8::Arguments& args)
{
	std::string NOT = "!";
	std::string AND = "&";
	std::string OR = "|";
	std::string XOR = "^";
	
	std::string op = c_stringify(args[0]);
	ibis::bitvector left;
	BVdecode(args[1],&left);
	if (op == NOT)
		left.flip();
	if (args.Length() == 3)
	{
		ibis::bitvector right;
		BVdecode(args[2],&right);
		if (op == AND)
			left &= right;
		else if (op == OR)
			left |= right;
		else if (op == XOR)
			left ^= right;
		else {
			v8::ThrowException(v8::Exception::TypeError(v8::String::New("unknown logical operation")));	
			return v8::Undefined();
		}
	}
	return BVencode(&left);
}

// random bitvector params
uint32_t randbv_intervals=25;
uint32_t randbv_length = 100;
uint32_t randbv_max = 10000;
v8::Handle<v8::Value> parse_randbv_params(v8::Handle<v8::Object> p)
{
	if (p->Has(v8::String::New("n")))
		randbv_intervals = p->Get(v8::String::New("n"))->Uint32Value();
	if (p->Has(v8::String::New("length")))
		randbv_length = p->Get(v8::String::New("length"))->Uint32Value();
	if (p->Has(v8::String::New("max")))
		randbv_max = p->Get(v8::String::New("max"))->Uint32Value();
	return v8::String::New("OK");
}
// random bitvector generator - sorted intervals
ibis::bitvector random_intervals(uint32_t n, uint32_t max, uint32_t len)
{
	ibis::array_t<uint32_t> starts;
	starts.resize(n);
	for(size_t i = 0; i < n; i++)
		starts[i] = rand() % max;
	ibis::array_t<uint32_t> ind;
	starts.sort(ind);
	// iterate over the sorted start positions
	// and build up a bitvector
	ibis::bitvector intspan;
	uint32_t st = starts[ind[0]];
	if (st > 0)
		intspan.appendFill(0,st);
	uint32_t end = st + len;
	for(uint32_t i=1;i<n;i++) {
		if (starts[ind[i]] > end) {
			intspan.appendFill(1,end-st);
			st = starts[ind[i]];
			intspan.appendFill(0,st-end);
		}
		end = starts[ind[i]] + len;
	}
	intspan.appendFill(1,end-st);
	return intspan;
}

v8::Handle<v8::Value> randbv(const v8::Arguments& args)
{
	// parse args
	v8::Handle<v8::Object> p = v8::Handle<v8::Object>::Cast(args[0]);
	v8::Handle<v8::Value> rc = parse_randbv_params(p);
	if (rc->IsUndefined())
		return rc;

	ibis::bitvector bv = random_intervals(randbv_intervals, randbv_max, randbv_length);
	return BVencode(&bv);
}

double bitvector_phi(ibis::bitvector &bv1, ibis::bitvector &bv2)
{
	if (bv1.size() < bv2.size()) {
		bv1.adjustSize(0,bv2.size());
	}
	else if (bv1.size() > bv2.size()) {
		bv2.adjustSize(0,bv1.size());
	}
	ibis::bitvector bvboth = ibis::bitvector(bv1);
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

v8::Handle<v8::Value> phi(const v8::Arguments& args)
{
	ibis::bitvector left, right;
	BVdecode(args[0],&left);
	BVdecode(args[1],&right);
	return v8::Number::New(bitvector_phi(left,right));
}

v8::Handle<v8::Value> MC(const v8::Arguments& args)
{
	ibis::bitvector A,B;
	BVdecode(args[0],&A);
	BVdecode(args[1],&B);
	uint32_t n = args[2]->Uint32Value();
	
	// keep only the fill words from A
	ibis::array_t<uint32_t> A_arr;
	A.write(A_arr);
	ibis::array_t<uint32_t> A_1fills;
	ibis::array_t<uint32_t> A_0fills;
	ibis::array_t<uint32_t> A_fills;
	for (uint32_t i=0; i<A_arr.size(); i++)
		if (A_arr[i] >= HEADER0) {
			A_fills.push_back(A_arr[i]);
			if (A_arr[i] >= HEADER1)
				A_1fills.push_back(A_arr[i]);
			else
				A_0fills.push_back(A_arr[i]);
		}
	ibis::bitvector A_1fills_bv = ibis::bitvector(A_1fills);
	ibis::bitvector A_0fills_bv = ibis::bitvector(A_0fills);
	ibis::bitvector A_fills_bv = ibis::bitvector(A_fills);

	// std::cerr << "A.cnt() " << A.cnt()
	// 	<< "\nA_1fills_bv.cnt() " << A_1fills_bv.cnt()
	// 	<< "\nA_0fills_bv.cnt() " << A_0fills_bv.cnt()
	// 	<< "\nA_fills_bv.cnt() " << A_fills_bv.cnt()
	// 	<< "\nA_1fills.size() " << A_1fills.size()
	// 	<< "\nA_0fills.size() " << A_0fills.size()
	// 	<< "\nA_fills.size() " << A_fills.size()
	// 	<< std::endl;

	// keep only the fill words from B
	ibis::array_t<uint32_t> B_arr;
	B.write(B_arr);
	ibis::array_t<uint32_t> B_1fills;
	ibis::array_t<uint32_t> B_0fills;
	ibis::array_t<uint32_t> B_fills;
	for (uint32_t i=0; i<B_arr.size(); i++)
		if (B_arr[i] >= HEADER0) {
			B_fills.push_back(B_arr[i]);
			if (B_arr[i] >= HEADER1)
				B_1fills.push_back(B_arr[i]);
			else
				B_0fills.push_back(B_arr[i]);
		}
	ibis::bitvector B_1fills_bv = ibis::bitvector(B_1fills);
	ibis::bitvector B_0fills_bv = ibis::bitvector(B_0fills);
	ibis::bitvector B_fills_bv = ibis::bitvector(B_fills);

	// std::cerr << "B.cnt() " << B.cnt()
	// 	<< "\nB_1fills_bv.cnt() " << B_1fills_bv.cnt()
	// 	<< "\nB_0fills_bv.cnt() " << B_0fills_bv.cnt()
	// 	<< "\nB_fills_bv.cnt() " << B_fills_bv.cnt()
	// 	<< "\nB_1fills.size() " << B_1fills.size()
	// 	<< "\nB_0fills.size() " << B_0fills.size()
	// 	<< "\nB_fills.size() " << B_fills.size()
	// 	<< std::endl;

	uint32_t o = B_fills_bv.count(A_fills_bv);

	uint32_t ge=0;
	for( uint32_t i=0; i < n; i++) {
		// create random permutations of the 4 sets of fill words
		for( uint32_t j=A_1fills.size(); j>0; j-- ) {
			uint32_t k = rand() % j;
			uint32_t t = A_1fills[k];
			A_1fills[k] = A_1fills[j-1];
			A_1fills[j-1] = t;
		}
		for( uint32_t j=B_1fills.size(); j>0; j-- ) {
			uint32_t k = rand() % j;
			uint32_t t = B_1fills[k];
			B_1fills[k] = B_1fills[j-1];
			B_1fills[j-1] = t;
		}
		for( uint32_t j=A_0fills.size(); j>0; j-- ) {
			uint32_t k = rand() % j;
			uint32_t t = A_0fills[k];
			A_0fills[k] = A_0fills[j-1];
			A_0fills[j-1] = t;
		}
		for( uint32_t j=B_0fills.size(); j>0; j-- ) {
			uint32_t k = rand() % j;
			uint32_t t = B_0fills[k];
			B_0fills[k] = B_0fills[j-1];
			B_0fills[j-1] = t;
		}

		// zipper the 0- and 1- fill words
		uint32_t k0 = 0;
		uint32_t k1 = 1;
		if (rand() % 2) {
			k0=1;
			k1=0;
		}
		for( uint32_t j=0; j < A_0fills.size(); j++ )
			A_fills[2*j+k0] = A_0fills[j];
		for( uint32_t j=0; j < A_1fills.size(); j++ )
			A_fills[2*j+k1] = A_1fills[j];
		k0 = 0;
		k1 = 1;
		if (rand() % 2) {
			k0=1;
			k1=0;
		}
		for( uint32_t j=0; j < B_0fills.size(); j++ )
			B_fills[2*j+k0] = B_0fills[j];
		for( uint32_t j=0; j < B_1fills.size(); j++ )
			B_fills[2*j+k1] = B_1fills[j];

		ibis::bitvector Arand = ibis::bitvector(A_fills);
		ibis::bitvector Brand = ibis::bitvector(B_fills);
		
		uint32_t c = Arand.count(Brand);
		if (c >= o)
			ge++;
	}

	return v8::Number::New(ge);
}
void Init(v8::Handle<v8::Object> target)
{
	srand ( time(NULL) );
	target->Set(v8::String::NewSymbol("SQL"), v8::FunctionTemplate::New(SQL)->GetFunction());
	target->Set(v8::String::NewSymbol("histogram"), v8::FunctionTemplate::New(histogram)->GetFunction());
	target->Set(v8::String::NewSymbol("scatter"), v8::FunctionTemplate::New(scatter)->GetFunction());
	target->Set(v8::String::NewSymbol("logical"), v8::FunctionTemplate::New(logical)->GetFunction());
	target->Set(v8::String::NewSymbol("cnt"), v8::FunctionTemplate::New(cnt)->GetFunction());
	target->Set(v8::String::NewSymbol("size"), v8::FunctionTemplate::New(size)->GetFunction());
	target->Set(v8::String::NewSymbol("bytes"), v8::FunctionTemplate::New(bytes)->GetFunction());
	target->Set(v8::String::NewSymbol("set2bvec"), v8::FunctionTemplate::New(set2bvec)->GetFunction());
	target->Set(v8::String::NewSymbol("intervals2bvec"), v8::FunctionTemplate::New(intervals2bvec)->GetFunction());
	target->Set(v8::String::NewSymbol("bvec2set"), v8::FunctionTemplate::New(bvec2set)->GetFunction());
	target->Set(v8::String::NewSymbol("bvec2span"), v8::FunctionTemplate::New(bvec2span)->GetFunction());
	target->Set(v8::String::NewSymbol("randbvec"), v8::FunctionTemplate::New(randbv)->GetFunction());
	target->Set(v8::String::NewSymbol("phi"), v8::FunctionTemplate::New(phi)->GetFunction());
	target->Set(v8::String::NewSymbol("MC"), v8::FunctionTemplate::New(MC)->GetFunction());
}

NODE_MODULE(fb, Init)