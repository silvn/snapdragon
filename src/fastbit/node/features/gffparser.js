/* gffParser.js
	An example script that extends the node-csv-parser for converting GFF to JSON. [ or BED ]
	jermth@me.com
*/


var csv = require('csv');
var argv = require('optimist')
    .usage('Parse a GFF file.\nUsage: $0 -f input_gff_file  [ -o outputformat (\'json\' | \'bed\') ]')
    .demand('f')
	.default('o', 'json')
    .alias('f', 'file')
	.describe('o', 'json or bed. Optional defaults to json')
    .argv;
	
var file=argv.f;
var commentRegex = /^\#/;
var keyValRegex = /\S+/g;

var alterGFFColumns = function(data){
	// Convert column 9 of the GTF file into key value pairs
	var col9Array=data.col9.split("\;");
	for( var i in col9Array ){
		if (keyValRegex.test(col9Array[i])){
			var keyValPair = col9Array[i].replace(/["']/g,"").match(keyValRegex);					
			data[keyValPair[0]]=keyValPair[1];
//			console.log(keyValPair[0],": ",keyValPair[1]);						
		}
	}	
	
	// Change from 1-based coordinates to 0-based half-open
	--data['start'];
	data['end'];
  	return data;		
	
	// TODO: I wanted to separate the col9 parsing from the coordinates transformation, 
	// but for some reason I can't chain more than 1 transformation to the csv object.
	// ie: I can't do gff = csv.from.path(file).transform(col9Array).transform(coordianteChange). etc
	// Only the last transformation is performed.
	// Perhaps the data array is not passed in the chain? 
}


var gff= csv().from.path(file, {
      	columns: ['chr','source','type','start','end','score','strand','phase','col9'],
	  	delimiter: "\t"	  
    });


if ( argv.o.toLowerCase() === 'json'){
    gff.transform(alterGFFColumns).on('record', function(data,index){
          console.log(JSON.stringify(data));
    });	
}else if ( argv.o.toLowerCase() === 'bed') {
    gff.transform(alterGFFColumns).to.stream(process.stdout, {
         	columns: ['chr', 'start','end','type','score',]
    });
}else {
	console.log("Only 'bed' and 'json' accepted as output formats");
	process.exit(1);
}


/* Adding chain for converting to BED6 format

*/
