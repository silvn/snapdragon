#include <node.h>
#include <ibis.h>
#include <category.h>
#include <iostream>

using namespace v8;

std::string c_stringify(Local<Value> lv) {
	String::Utf8Value p(lv);
	return std::string(*p);
}

// ibis::bitvector decodeBV(Local<Value> lv) {
// 	std::string s = c_stringify(lv);
// 	
// }

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

ibis::table* run_query(Handle<Object> p)
{
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
	std::string select_str = c_stringify(p->Get(String::New("select")));

	ibis::table *res = tbl->select(select_str.c_str(), query);
	
	if (p->Has(String::New("orderby"))) {
		std::string order_by = c_stringify(p->Get(String::New("orderby")));
		res->orderby(order_by.c_str());
	}
	return res;
}

Handle<Value> SQL(const Arguments& args)
{
	HandleScope scope;
	if(args.Length() != 1)
	{
		ThrowException(Exception::TypeError(String::New("Wrong number of arguments")));
		return scope.Close(Undefined());
	}
	Handle<Object> p = Handle<Object>::Cast(args[0]);
	if (! p->Has(String::New("from"))) {
		ThrowException(Exception::TypeError(String::New("Missing 'from' argument")));
		return scope.Close(Undefined());
	}
	if (! p->Has(String::New("select"))) {
		ThrowException(Exception::TypeError(String::New("Missing 'select' argument")));
		return scope.Close(Undefined());
	}


	ibis::table *res = run_query(p);

	Local<Object> jsObj = TabletoJs(*res);
	delete res;
	return scope.Close(jsObj);
}

// 1D histogram
// This function uses ibis::part histogram functions
// so if you start with multiple partitions you have to
// preprocess them into a single in-memory table (with one partition).
// Another reason to preprocess is if the select is not just a column name
//
// If the histogram is still subject to constraints, do a countQuery first to get the mask
// binning is either adaptive or uniform
// return array of counts at either bin centers or intervals
Handle<Value> histogram(const Arguments& args)
{
	HandleScope scope;
	return scope.Close(String::New("unimplemented"));
}

// 2D histogram
Handle<Value> scatter(const Arguments& args)
{
	HandleScope scope;
	return scope.Close(String::New("unimplemented"));
	
}

// logical operations
Handle<Value> logical(const Arguments& args)
{
	HandleScope scope;
	
	// ibis::bitvector bvec = decodeBV(args[1]);
	return scope.Close(String::New("unimplemented"));
	
}

Handle<Value> cnt(const Arguments& args)
{
	HandleScope scope;
	return scope.Close(String::New("unimplemented"));
	
}

Handle<Value> size(const Arguments& args)
{
	HandleScope scope;
	return scope.Close(String::New("unimplemented"));
}

void Init(Handle<Object> target)
{
	target->Set(String::NewSymbol("SQL"), FunctionTemplate::New(SQL)->GetFunction());
	target->Set(String::NewSymbol("histogram"), FunctionTemplate::New(histogram)->GetFunction());
	target->Set(String::NewSymbol("scatter"), FunctionTemplate::New(scatter)->GetFunction());
	target->Set(String::NewSymbol("logical"), FunctionTemplate::New(logical)->GetFunction());
	target->Set(String::NewSymbol("cnt"), FunctionTemplate::New(cnt)->GetFunction());
	target->Set(String::NewSymbol("size"), FunctionTemplate::New(size)->GetFunction());
}

NODE_MODULE(fb, Init)