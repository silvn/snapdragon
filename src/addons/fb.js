var fb = require('./build/Release/fb');

console.log(fb.SQL({
	select: "GO_term,count(*) as tally",
	from: "/Users/olson/src/garden/iris/examples/fastbit/gene2GO",
	where: "49000 < gene_id < 50000",
	orderby: "tally desc"}));

console.log(fb.SQL({
	select: "GO_branch,count(*) as tally",
	from: "/Users/olson/src/garden/iris/examples/fastbit/gene2GO",
	where: "1=1",
	orderby: "tally desc"}));
