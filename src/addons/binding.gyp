{
  "targets": [
    {
      "target_name": "fb",
      "sources": [ "fb.cc" ],
	  "include_dirs": ["/Users/olson/src/garden/snapdragon/include"],
	  "libraries": ["-L/Users/olson/src/garden/snapdragon/lib",'-lfastbit'],
	  "conditions": [
		  ['OS=="mac"', {
			  'xcode_settings': {
				  'GCC_ENABLE_CPP_EXCEPTIONS': 'YES'
			  }
		  }]
	  ]
    }
  ]
}
