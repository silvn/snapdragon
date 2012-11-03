var fb = require('./build/Release/fb');

console.log(fb.SQL(
	{
		select: "GO_branch,count(*) as tally",
		from: "/Users/olson/src/garden/iris/examples/fastbit/gene2GO",
		orderby: "tally desc"
	}
));
console.log(fb.SQL(
	{
		select: "GO_branch,count(*) as tally",
		from: "/Users/olson/src/garden/iris/examples/fastbit/gene2GO",
		where: "49000 < gene_id < 50000",
		orderby: "tally desc"
	}
));

console.log(fb.histogram(
	{
		select: "score",
		from: "/Users/olson/src/garden/iris/examples/fastbit/GWAS/3396",
		begin: 0,
		end: 10,
		stride: 1
	}
));
console.log(fb.histogram(
	{
		select: "score",
		from: "/Users/olson/src/garden/iris/examples/fastbit/GWAS/3396",
		adaptive: true,
		nbins: 10
	}
));
console.log(fb.scatter(
	{
		select: "pos,score",
		from: "/Users/olson/src/garden/iris/examples/fastbit/GWAS/3396/1",
		adaptive: true,
		nbins1: 10,
		nbins2: 10
	}
));
var chr_list = fb.SQL(
	{
		select: "chr,max(score) as score,max(pos) as pos",
		from: "/Users/olson/src/garden/iris/examples/fastbit/GWAS/3396"
	}
);
console.log(chr_list);

var n = chr_list["chr"].length;
for( var i=0; i < n; i++) {
	var args = {
		select: "pos,score",
		from: "/Users/olson/src/garden/iris/examples/fastbit/GWAS/3396/" + chr_list["chr"][i],
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

var bvec = fb.set2bvec([1,3,4,4,4,4,5,77,899,7654]);
console.log(bvec);
console.log(fb.cnt(bvec));
console.log(fb.size(bvec));
