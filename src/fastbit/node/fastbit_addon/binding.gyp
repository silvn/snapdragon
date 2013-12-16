{
  "targets": [
    {
      "target_name": "fb",
      "sources": [ "fb.cc" ],
	  "xxxinclude_dirs": ["../../../../vendor/fastbit/src"],
	  "xxxlibraries": ["-L../../../../../vendor/fastbit/src/.libs",'-lfastbit'],
	  "include_dirs": ["/usr/local/include"],
	  "libraries": ["-L/usr/local/lib",'-lfastbit'],
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
