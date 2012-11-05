var fb = require('./build/Release/fb');
var data = "/Users/olson/src/garden/iris/examples/fastbit/";
console.log(fb.SQL(
	{
		select: "GO_branch,count(*) as tally",
		from: data + "gene2GO",
		orderby: "tally desc"
	}
));
console.log(fb.SQL(
	{
		select: "GO_branch,count(*) as tally",
		from: data + "gene2GO",
		where: "49000 < gene_id < 50000",
		orderby: "tally desc"
	}
));

console.log(fb.histogram(
	{
		select: "score",
		from: data + "GWAS/3396",
		begin: 0,
		end: 10,
		stride: 1
	}
));
console.log(fb.histogram(
	{
		select: "score",
		from: data + "GWAS/3396",
		adaptive: true,
		nbins: 10
	}
));
console.log(fb.scatter(
	{
		select: "pos,score",
		from: data + "GWAS/3396/1",
		adaptive: true,
		nbins1: 10,
		nbins2: 10
	}
));
var chr_list = fb.SQL(
	{
		select: "chr,max(score) as score,max(pos) as pos",
		from: data + "GWAS/3396"
	}
);
console.log(chr_list);

var n = chr_list["chr"].length;
for( var i=0; i < n; i++) {
	var args = {
		select: "pos,score",
		from: data + "GWAS/3396/" + chr_list["chr"][i],
		begin1: 0,
		end1: chr_list["pos"][i],
		stride1: chr_list["pos"][i]/10,
		begin2: 0,
		end2: chr_list["score"][i],
		stride2: chr_list["score"][i]/10
	};
	console.log(args);
	console.log(fb.scatter(args));
}

var arr1 = [2,4,5,77,899,7654];
var bvec1 = fb.set2bvec(arr1);
var arr2 = [1,2,3,4,5,77,899,7654];
var bvec2 = fb.set2bvec(arr2);
var bvand = fb.logical('^',bvec1,bvec2);
console.log("arr1",arr1);
console.log("arr2",arr2);
console.log("xor",fb.bvec2set(bvand));
