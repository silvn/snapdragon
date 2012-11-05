#include <node.h>
#include <ibis.h>
#include <category.h>
#include <iostream>

using namespace v8;

std::string c_stringify(Local<Value> lv) {
	String::Utf8Value p(lv);
	return std::string(*p,p.length());
}

Local<Object> TabletoJs(const ibis::table &tbl)
{
	size_t nr = tbl.nRows();
	ibis::table::stringList nms = tbl.columnNames();
	ibis::table::typeList tps = tbl.columnTypes();

	Local<Object> JSON = Object::New();

	std::vector<const ibis::part*> parts;
	tbl.getPartitions(parts);
	if (parts.size() != 1) {
		std::cerr << "parts.size() = " << parts.size() << " expected 1" << std::endl;
		return JSON;
	}
	int ierr = 0;
	for (size_t j = 0; j < nms.size(); ++ j) {
		Local<Array> column = Array::New();
		switch (tps[j]) {
			case ibis::BYTE:
			case ibis::SHORT:
			case ibis::INT:
			{
				int32_t *buf = new int32_t[nr];
				ierr = tbl.getColumnAsInts(nms[j], buf);
				for (size_t i = 0; i < nr; ++ i)
					column->Set(i,Number::New(buf[i]));
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
						column->Set(i,Number::New(buf[i]));
				else
					for (size_t i = 0; i < nr; ++ i)
						column->Set(i,String::New((*dic)[buf[i]]));
				break;
			}
			case ibis::FLOAT:
			case ibis::DOUBLE:
			{
				double *buf = new double[nr];
				ierr = tbl.getColumnAsDoubles(nms[j], buf);
				for (size_t i = 0; i < nr; ++ i)
					column->Set(i,Number::New(buf[i]));
				break;
			}
			case ibis::TEXT:
			case ibis::CATEGORY:
			{
			    std::vector<std::string>* buf = new std::vector<std::string>();
				ierr = tbl.getColumnAsStrings(nms[j], *buf);
				for (size_t i = 0; i < nr; ++ i)
					column->Set(i,String::New((*buf)[i].c_str()));
				break;
			}
			default:
			{
				break;
			}
		}
		JSON->Set(String::New(nms[j]),column);
	}
	return JSON;
}

Handle<Value> run_query(Handle<Object> p, ibis::table*& res)
{
	if (! p->Has(String::New("from"))) {
		ThrowException(Exception::TypeError(String::New("Missing 'from' argument")));
		return Undefined();
	}
	if (! p->Has(String::New("select"))) {
		ThrowException(Exception::TypeError(String::New("Missing 'select' argument")));
		return Undefined();
	}

	std::string data_dir = c_stringify(p->Get(String::New("from")));
	ibis::table* tbl = ibis::table::create(data_dir.c_str());
	
	const ibis::qExpr* query;
	if (p->Has(String::New("where"))) {
		std::string query_cnd = c_stringify(p->Get(String::New("where")));
		ibis::whereClause tmp = ibis::whereClause(query_cnd.c_str());
		query = tmp.getExpr()->dup();
	}
	else {
		ibis::whereClause tmp = ibis::whereClause("1=1");
		query = tmp.getExpr()->dup();
	}
	
	// check for qExpr
	
	std::string select_str = c_stringify(p->Get(String::New("select")));
	res = tbl->select(select_str.c_str(), query);
	
	if (p->Has(String::New("orderby"))) {
		std::string order_by = c_stringify(p->Get(String::New("orderby")));
		res->orderby(order_by.c_str());
	}
	return String::New("OK");
}

// histogram params
bool adaptive = false;
uint32_t nbins=25;
double begin,end,stride;
Handle<Value> parse_histogram_params(Handle<Object> p)
{
	if (p->Has(String::New("adaptive")))
		if (p->Get(String::New("adaptive"))->IsTrue()) {
			adaptive = true;
			if (p->Has(String::New("nbins")))
				nbins = p->Get(String::New("nbins"))->Uint32Value();
		}

	if (! adaptive) {
		if (!p->Has(String::New("begin"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'begin' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("end"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'end' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("stride"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'stride' argument")));
			return Undefined();
		}
		begin = p->Get(String::New("begin"))->NumberValue();
		end = p->Get(String::New("end"))->NumberValue();
		stride = p->Get(String::New("stride"))->NumberValue();
	}
	return String::New("OK");
}

// 1D histogram
// returns a Javascript object {bounds:[], counts:[]}
// This function uses ibis::part histogram functions
// so we preprocess them into a single in-memory table (with one partition).
//
// binning is either adaptive or uniform
Handle<Value> histogram(const Arguments& args)
{
	HandleScope scope;

	// parse args
	Handle<Object> p = Handle<Object>::Cast(args[0]);
	Handle<Value> rc = parse_histogram_params(p);
	if (rc->IsUndefined())
		return scope.Close(rc);

	// preprocess
	ibis::table *res = 0;
	rc = run_query(p,res);
	if (rc->IsUndefined())
		return scope.Close(rc);

	// no data after querying
	if (res->nRows() <= 0)
		return scope.Close(Undefined());

	ibis::table::stringList nms = res->columnNames();

	std::vector<const ibis::part*> parts;
	res->getPartitions(parts);
	if (parts.size() != 1) {
		ThrowException(Exception::TypeError(String::New("WTF! expected one partition after preprocessing")));
		return scope.Close(Undefined());
	}

	long ierr;
	Local<Object> JSON = Object::New();
	Local<Array> v8bounds = Array::New();
	Local<Array> v8counts = Array::New();
	if (adaptive) {
		std::vector<double> bounds;
		std::vector<uint32_t> counts;
		ierr = parts[0]->get1DDistribution("1=1",nms[0],nbins,bounds,counts);
		if (ierr < 0) {
			ThrowException(Exception::TypeError(String::New("adaptive 1D Distribution error")));
			return scope.Close(Undefined());
		}
		// format results
		for(size_t i=0; i < counts.size(); i++) {
			v8bounds->Set(i,Number::New(bounds[i]));
			v8counts->Set(i,Number::New(counts[i]));
		}
	}
	else {
		std::vector<uint32_t> counts;
		ierr = parts[0]->get1DDistribution("1=1",nms[0],begin,end,stride,counts);
		if (ierr < 0) {
			ThrowException(Exception::TypeError(String::New("uniform 1D Distribution error")));
			return scope.Close(Undefined());
		}
		// format results
		double pos = begin;
		for(size_t i=0; i < counts.size(); i++) {
			v8bounds->Set(i,Number::New(pos));
			v8counts->Set(i,Number::New(counts[i]));
			pos += stride;
		}
	}
	JSON->Set(String::New("bounds"),v8bounds);
	JSON->Set(String::New("counts"),v8counts);
	return JSON;
}

// 2D histogram params
uint32_t nbins1=25;
uint32_t nbins2=25;
double begin1,end1,stride1,begin2,end2,stride2;
Handle<Value> parse_scatter_params(Handle<Object> p)
{
	if (p->Has(String::New("adaptive")))
		if (p->Get(String::New("adaptive"))->IsTrue()) {
			adaptive = true;
			if (p->Has(String::New("nbins1")))
				nbins1 = p->Get(String::New("nbins1"))->Uint32Value();
			if (p->Has(String::New("nbins2")))
				nbins2 = p->Get(String::New("nbins2"))->Uint32Value();
		}

	if (! adaptive) {
		if (!p->Has(String::New("begin1"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'begin1' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("end1"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'end1' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("stride1"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'stride1' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("begin2"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'begin2' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("end2"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'end2' argument")));
			return Undefined();
		}
		if (!p->Has(String::New("stride2"))) {
			ThrowException(Exception::TypeError(String::New("Missing 'stride2' argument")));
			return Undefined();
		}
		begin1 = p->Get(String::New("begin1"))->NumberValue();
		end1 = p->Get(String::New("end1"))->NumberValue();
		stride1 = p->Get(String::New("stride1"))->NumberValue();
		begin2 = p->Get(String::New("begin2"))->NumberValue();
		end2 = p->Get(String::New("end2"))->NumberValue();
		stride2 = p->Get(String::New("stride2"))->NumberValue();
	}
	return String::New("OK");
}

// 2D histogram
// returns a Javascript object {bounds1:[], bounds2:[], counts:[]}
// This function uses ibis::part histogram functions
// so we preprocess them into a single in-memory table (with one partition).
//
// binning is either adaptive or uniform
Handle<Value> scatter(const Arguments& args)
{
	HandleScope scope;

	// parse args
	Handle<Object> p = Handle<Object>::Cast(args[0]);
	Handle<Value> rc = parse_scatter_params(p);
	if (rc->IsUndefined())
		return scope.Close(rc);

	// preprocess
	ibis::table *res=0;
	rc = run_query(p,res);
	if (rc->IsUndefined())
		return scope.Close(rc);

	// no data after querying
	if (res->nRows() <= 0)
		return scope.Close(Undefined());

	ibis::table::stringList nms = res->columnNames();

	std::vector<const ibis::part*> parts;
	res->getPartitions(parts);
	if (parts.size() != 1) {
		ThrowException(Exception::TypeError(String::New("WTF! expected one partition after preprocessing")));
		return scope.Close(Undefined());
	}

	long ierr;
	Local<Object> JSON = Object::New();
	Local<Array> v8bounds1 = Array::New();
	Local<Array> v8bounds2 = Array::New();
	Local<Array> v8counts = Array::New();
	std::vector<double> bounds1;
	std::vector<double> bounds2;
	std::vector<uint32_t> counts;
	if (adaptive) {
		ierr = parts[0]->get2DDistribution(nms[0],nms[1],nbins1,nbins2,bounds1,bounds2,counts);
		if (ierr < 0) {
			ThrowException(Exception::TypeError(String::New("adaptive 2D Distribution error")));
			return scope.Close(Undefined());
		}
		// format results
		for(size_t i=0; i < counts.size(); i++)
			v8counts->Set(i,Number::New(counts[i]));
		for(size_t i=0; i < bounds1.size(); i++)
			v8bounds1->Set(i,Number::New(bounds1[i]));
		for(size_t i=0; i < bounds2.size(); i++)
			v8bounds2->Set(i,Number::New(bounds2[i]));
	}
	else {
		std::vector<uint32_t> counts;
		ierr = parts[0]->get2DDistribution("1=1",nms[0],begin1,end1,stride1,nms[1],begin2,end2,stride2,counts);
		if (ierr < 0) {
			ThrowException(Exception::TypeError(String::New("uniform 2D Distribution error")));
			return scope.Close(Undefined());
		}
		// format results
		for(size_t i=0; i < counts.size(); i++)
			v8counts->Set(i,Number::New(counts[i]));
		double pos = begin1;
		for(size_t i=0; i < bounds1.size(); i++) {
			v8bounds1->Set(i,Number::New(pos));
			pos += stride1;
		}
		pos = begin2;
		for(size_t i=0; i < bounds2.size(); i++) {
			v8bounds2->Set(i,Number::New(pos));
			pos += stride2;
		}
	}
	JSON->Set(String::New("counts"),v8counts);
	JSON->Set(String::New("bounds1"),v8bounds1);
	JSON->Set(String::New("bounds2"),v8bounds2);
	return JSON;
}

// generic simple SQL function
Handle<Value> SQL(const Arguments& args)
{
	HandleScope scope;
	if(args.Length() != 1)
	{
		ThrowException(Exception::TypeError(String::New("Wrong number of arguments")));
		return scope.Close(Undefined());
	}

	Handle<Object> p = Handle<Object>::Cast(args[0]);

	ibis::table *res=0;
	Handle<Value> rc = run_query(p,res);

	if (rc->IsUndefined())
		return scope.Close(rc);

	Local<Object> jsObj = TabletoJs(*res);
	delete res;
	return scope.Close(jsObj);
}

ibis::bitvector BVdecode(Local<Value> lv) {
	ibis::array_t<uint32_t> arr;
	Local<Array> lva = Local<Array>::Cast(lv);
	for(size_t i=0; i < lva->Length(); i++)
		arr.push_back(lva->Get(i)->Uint32Value());
	ibis::bitvector bv(arr);
	return bv;
}

Local<Value> BVencode(ibis::bitvector bv) {
	ibis::array_t<uint32_t> arr;
	bv.write(arr);
	Local<Array> v8_arr = Array::New();
	for(size_t i=0; i < arr.size(); i++)
		v8_arr->Set(i,Uint32::New(arr[i]));
	return v8_arr;
}

Handle<Value> cnt(const Arguments& args)
{
	HandleScope scope;
	ibis::bitvector bv = BVdecode(args[0]);
	return scope.Close(Integer::New(bv.cnt()));
}

Handle<Value> size(const Arguments& args)
{
	HandleScope scope;
	ibis::bitvector bv = BVdecode(args[0]);
	return scope.Close(Integer::New(bv.size()));
}

Handle<Value> set2bvec(const Arguments& args)
{
	HandleScope scope;
	if (! args[0]->IsArray())
		return scope.Close(String::New("first and only argument needs to be an array"));
	Handle<Array> set = Handle<Array>::Cast(args[0]);
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
	return scope.Close(BVencode(bv));
}

Handle<Value> bvec2set(const Arguments& args)
{
	HandleScope scope;
	ibis::bitvector bv = BVdecode(args[0]);
	size_t i=0;
	Local<Array> set = Array::New();
	ibis::bitvector::indexSet index = bv.firstIndexSet();
	while (index.nIndices() > 0) {
		const ibis::bitvector::word_t *idx0 = index.indices();
		if (index.isRange()) {
			for(size_t j = 0; j < index.nIndices(); j++) {
				set->Set(i,Number::New(*idx0+j));
				i++;
			}
		}
		else {
		    for (ibis::bitvector::word_t j = 0; j<index.nIndices(); ++j) {
				set->Set(i,Number::New(idx0[j]));
				i++;
		    }
		}
		++ index;
	}

	return scope.Close(set);
}

// logical operations
Handle<Value> logical(const Arguments& args)
{
	std::string NOT = "!";
	std::string AND = "&";
	std::string OR = "|";
	std::string XOR = "^";
	
	HandleScope scope;
	std::string op = c_stringify(args[0]);
	ibis::bitvector left = BVdecode(args[1]);
	if (op == NOT)
		left.flip();
	if (args.Length() == 3)
	{
		ibis::bitvector right = BVdecode(args[2]);
		if (op == AND)
			left &= right;
		else if (op == OR)
			left |= right;
		else if (op == XOR)
			left ^= right;
		else {
			ThrowException(Exception::TypeError(String::New("unknown logical operation")));	
			return scope.Close(Undefined());
		}
	}
	return scope.Close(BVencode(left));
}


void Init(Handle<Object> target)
{
	target->Set(String::NewSymbol("SQL"), FunctionTemplate::New(SQL)->GetFunction());
	target->Set(String::NewSymbol("histogram"), FunctionTemplate::New(histogram)->GetFunction());
	target->Set(String::NewSymbol("scatter"), FunctionTemplate::New(scatter)->GetFunction());
	target->Set(String::NewSymbol("logical"), FunctionTemplate::New(logical)->GetFunction());
	target->Set(String::NewSymbol("cnt"), FunctionTemplate::New(cnt)->GetFunction());
	target->Set(String::NewSymbol("size"), FunctionTemplate::New(size)->GetFunction());
	target->Set(String::NewSymbol("set2bvec"), FunctionTemplate::New(set2bvec)->GetFunction());
	target->Set(String::NewSymbol("bvec2set"), FunctionTemplate::New(bvec2set)->GetFunction());
}

NODE_MODULE(fb, Init)