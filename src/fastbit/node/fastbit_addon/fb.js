var fb = require('./build/Release/fb');
var data = "/Users/olson/src/garden/iris/examples/fastbit/";
var args = {
	select: "GO_branch,count(*) as tally",
	from: data + "gene2GO",
	orderby: "tally desc"
};
console.log("fb.SQL",args);
console.log(fb.SQL(args));

args = {
	select: "GO_branch,count(*) as tally",
	from: data + "gene2GO",
	where: "49000 < gene_id < 50000",
	orderby: "tally desc"
};
console.log("fb.SQL",args);
console.log(fb.SQL(args));

args = {
	select: "score",
	from: data + "GWAS/3396",
	begin: 0,
	end: 10,
	stride: 1
};
console.log("fb.histogram",args);
console.log(fb.histogram(args));
args = {
	select: "score",
	from: data + "GWAS/3396",
	adaptive: true,
	nbins: 10
};
console.log("fb.histogram",args);
console.log(fb.histogram(args));
args = {
	select: "pos,score",
	from: data + "GWAS/3396/1",
	adaptive: true,
	nbins1: 10,
	nbins2: 10
};
console.log("fb.scatter",args);
console.log(fb.scatter(args));

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

var arr1 = [1,2,4,5,77,78,899,7654];
var bvec1 = fb.set2bvec(arr1);
console.log("bvec1",bvec1);
var arr2 = [1,2,3,4,5,77,899,7654];
var bvec2 = fb.set2bvec(arr2);
console.log("arr1",arr1,fb.bvec2set(bvec1));
console.log("arr2",arr2);
console.log("xor",fb.bvec2set(fb.logical('^',bvec1,bvec2)));

var start = new Date;
var bvrand1 = fb.randbvec({n:1e5,max:1000000000,length:100});
var bvrand2 = fb.randbvec({n:1e5,max:1000000000,length:200});
console.log("build 2 random bitvectors",new Date - start,"ms");
console.log("bytes(bvrand1)=",fb.bytes(bvrand1));
console.log("bytes(bvrand2)=",fb.bytes(bvrand2));
console.log("cnt(bvrand1)=",fb.cnt(bvrand1));
console.log("cnt(bvrand2)=",fb.cnt(bvrand2));
start = new Date;
console.log("cnt(logical('&',bvrand1,bvrand2))=",fb.cnt(fb.logical('&',bvrand1,bvrand2)));
console.log(new Date - start,"ms");
start = new Date;
console.log("phi(bvrand1,bvrand2)=",fb.phi(bvrand1,bvrand2));
console.log(new Date - start,"ms");

bvec1 = fb.intervals2bvec([100,300,800,30],[200,400,900,50]);
bvec2 = fb.intervals2bvec([100,700,10],[400,900,50]);
var MC1 = fb.MC(bvec1,bvec2,100);
console.log(MC1);
// var iterations = [1,100];
// var sizes = [1e5,1e6,1e7];
// for (var i=0;i<sizes.length;i++) {
// 	var size = sizes[i];
// 	var bv1 = fb.randbvec({n:size,max:1000000000,length:100});
// 	var bv2 = fb.randbvec({n:size,max:1000000000,length:100});
// 	for (var j=0;j<iterations.length;j++) {
// 		var iter = iterations[j];
// 		start = new Date;
// 		var MC = fb.MC(bv1,bv2,iter);
// 		console.log(size,iter,MC,new Date - start,"ms");
// 	}
// }
