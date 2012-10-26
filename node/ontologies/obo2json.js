
var reader = require ("buffered-reader").DataReader;

var argv = require('optimist')
    .usage('Converts a OBO file into a set of json objects.\nUsage: $0')
    .demand('f')
    .alias('f', 'file')
    .argv;
	

var inStanza='';	
var inTerms=0;
var stanzaTagRegexp= /^\[\w*\]$/;

var ontologyParent={ file: argv.file };
var stanzaObject = {};
		
var fs = require('fs');
var termStream = fs.createWriteStream("terms.json");
var parentStream = fs.createWriteStream("parentOntology.json");

new reader (argv.file, { encoding: "utf8" })

		.on ("error", function (error){
			console.log (error);
		})
		.on ("line", function (line, nextByteOffset){

			if (line.length>0){
				line=line.replace(/["]/g,''); // strip out the double quotes
				
				if(stanzaTagRegexp.test(line)){ // We're hitting a stanza tag
				
					// assign the current inStanza variable to the stanza tag
					inStanza=line.replace(/[\[\]]/g,'');
					//console.log("line:", line);	 
					//reset the current stanzaObject				
					stanzaObject = {};										
					
				}else if (inStanza.length ===0){
					// we're in the header lines at the beginning of the file.
					// cpature the key values as a parent ontology object
					var keyValArray=parseKeyVal(line);			
					ontologyParent=addKeyVal(ontologyParent,keyValArray);
					

				}else{
					// After the header portion of the OBO file, everything else
					// should belong to a stanza.					

					var keyValArray=parseKeyVal(line);
					if (line.match(/^id/)){
						stanzaObject = Entry(keyValArray[1]);
					}else{
						stanzaObject = addKeyVal(stanzaObject,keyValArray);
					}
				}				
			}else{
				if (!is_empty(stanzaObject)){
					// if we already have a stanza object, print it out.
					printJSON(inStanza,stanzaObject);
				}	
				inStanza='';
			}		
		})
		.on ("end", function (){
			printJSON(inStanza,stanzaObject); // last stanzaObject
			printJSON('ontology',ontologyParent); // print the ontology parent
			termStream.end();
			parentStream.end();
		})
		.read ();
	
	
var parseKeyVal = function(line){
	var splitString=line.split(":");	
	var key=splitString.shift();
	
	var value=splitString.join(":").replace(/^\s+/,'');

	if (value.match(/ \! /)){
		value=value.split(" ! ")[0];
	}	
	return [key,value];
}

var Entry = function(id) {
	return { id: id,
			 ontology: ontologyParent.ontology};	
}

var addKeyVal = function (that,keyValPair){
	
	if (keyValPair[0] === 'xref'){
		var xrefKeyPair=parseKeyVal(keyValPair[1]);
		if (typeof(that.xref) === "undefined") {
			that.xref = {};
		}
		that.xref=addKeyVal(that.xref,xrefKeyPair);
	}else{		
		if (that.hasOwnProperty(keyValPair[0])){
			if (typeof that[keyValPair[0]]  === 'string' ){
				var tmp= that[keyValPair[0]];
				that[keyValPair[0]]=[];
				that[keyValPair[0]].push(tmp);
				that[keyValPair[0]].push(keyValPair[1]);
			}else{
				that[keyValPair[0]].push(keyValPair[1]);
			}
		}else{								
			that[keyValPair[0]] = keyValPair[1];
		}
	}	
	return that;
}

var printJSON = function (type,that){
	//TODO: print different term types to separate files
	if (type ==='Term'){
		termStream.write(JSON.stringify(that));
		termStream.write("\n");
	}else if (type ==='ontology'){
		parentStream.write(JSON.stringify(that));		
	}
}

var hasOwnProperty = Object.prototype.hasOwnProperty;

function is_empty(obj) {
    // Assume if it has a length property with a non-zero value
    // that that property is correct.
    if (obj.length && obj.length > 0)    return false;
    if (obj.length && obj.length === 0)  return true;

    for (var key in obj) {
        if (hasOwnProperty.call(obj, key))    return false;
    }

    return true;
}