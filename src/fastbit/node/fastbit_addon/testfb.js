var bins = 10;
var fb = require('./build/Release/fb');
var data = "/Users/olson/Desktop/glycine_max";
var args = {
	select: "feature_type,max(start) as ms",
	from: data
};
console.log("fb.SQL",args);
var features = fb.SQL(args);
console.log("features",features);

var n = features["feature_type"].length;
for( var i=0; i < n; i++) {
    var f = features["feature_type"][i];
    var m = features["ms"][i];
    console.log(i,f,m);
    if (m>1) {
        args = {
            select: "start",
            from:   data,
            begin:  0,
            end:    m,
            stride: m/bins,
            where: 'seqid="1" and feature_type="' + f + '"'
        };
        console.log("fb.histogram",args);
        console.log(fb.histogram(args));
    }
}

// args = {
//     select: "start",
//     from: data,
//     begin: 0,
//     end: 50000000,
//     stride: 5000000,
//     where: 'seqid="1" and feature_type="gene"'
// };
// console.log("fb.histogram",args);
// console.log(fb.histogram(args));
// 
