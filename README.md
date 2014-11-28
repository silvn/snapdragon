#Snapdragon

Version 0.0.2

[![Build Status](https://travis-ci.org/silvn/snapdragon.png)](https://travis-ci.org/silvn/snapdragon)

Biological analytics made fast and easy.

##Installation

    git submodule update --init
    ./configure
    make
    make check
    make install

##Development Notes

###Releases

We use `git flow` for managing releases. To create a new release, the following actions must be taken:

 1. `git flow release start X.Y.Z`
 2. Increment version at the top of this `README.md` file
 3. Increment version in `configure.ac` (AC_INIT directive)
 4. Run `autoreconf -vfi` to regenerated `configure`
 5. Commit changes
 6. `git flow release finish X.Y.Z`

##License

    Copyright (c) 2012-2014 Jer-Ming Chia <jermth at gmail.com>
    Copyright (c) 2012-2014 Andrew Olson <aolson at me.com>
    Copyright (c) 2012-2014 Shiran Pasternak <shiranpasternak at gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to
    deal in the Software without restriction, including without limitation the
    rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions: The above
    copyright notice and this permission notice shall be included in all copies
    or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
