## Node JS parser that converts ontology OBO flat files into JSON

- Produces 2 text files:
	1. terms.json
	- all the terms stanzas in the OBO files as JSON 
	2. parentOntology.json
	- A single JSON representing the ontology these terms belong to.

- These files can be loaded directly into a mongodb:
	- e.g: mongoimport --db DBNAME --collection COLLECTION_NAME --file terms.json

### Running with test data as an example:
$ node obo2json.js -f ../../test_data/gene_ontology_obo_100terms.obo

- this produces "parentOntology.json" and "terms.json"
