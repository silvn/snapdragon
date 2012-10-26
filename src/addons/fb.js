var fb = require('./build/Release/fb');

console.log(fb.SQL({
	select: "GO_branch,count(*) as tally",
	from: "/Users/olson/src/garden/iris/examples/fastbit/gene2GO",
	where: "gene_id < 50000",
	orderby: "tally desc"}));
