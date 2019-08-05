# SLP_SQP   `slp_trust.m` and `sqp.m` MATLAB optimization functions

Copyright (c) 2015, Robert A. Canfield. All rights reserved.
* Virginia Tech and Air Force Institute of Technology
* bob.canfield@vt.edu
* <http://www.aoe.vt.edu/people/faculty/canfield.html>
* See accompanying LICENSE.txt file for conditions.

## m-Files
-  slp_trust - Sequential Linear Programming (SLP) with a Trust Region Strategy
-  sqp       - Schittkowski's Sequential Quadratic Programming method
-  foptions  - Default parameters used by the optimization routines.
-  Contents  - Text document describing files

## Documentation
-  sqp.pdf     - User's Guide for sqp and slp_trust (Adobe portable document format)
-  README.md   - ReadMe text file (markdown format)
-  LICENSE.txt - Open Source License notice

## Examples
-  run*.m    - Scripts to run example problems
-  f*.m      - Functions to evaluate objective, f, and constraints, g
-  g*.m      - Functions to evaluate gradients of f and g

## Private folder
-  Utility functions called by sqp.m and slp_trust.m

## Open Source License
Copyright (c) 2015, Robert A. Canfield. All rights reserved.

`sqp.m, slp_trust.m` MATLAB package

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal with the Software without restriction, including without 
limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so, subject to the following 
conditions:

* Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimers.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimers in the
  documentation and/or other materials provided with the distribution.

* Neither the names of Robert A. Canfield, Virginia Tech, Mark Spillman,
  Air Force Institute of Technology, nor the names of its contributors 
  may be used to endorse or promote products derived from this Software 
  without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS WITH THE SOFTWARE.