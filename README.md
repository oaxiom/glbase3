# README #

### What is glbase3? ###

glbase is a flexible and multifunctional toolkit allowing the user to perform many common analyses on ChIP-seq, microarray and RNA-seq data.

### How do I get set up? ###

I strongly recommend installing the mercurial version of glbase. glbase is updated regularly, sometimes daily to remove bugs and add features. 

You will need Python3.6, numpy, scipy,  matplotlib, networkx, graphviz, pydot and sklearn. On the Mac I use Macports or Anaconda. Linux machines should be pretty easy to do this on, more than likely you already have numpy, scipy and matplotlib (and all major linux distributions come with Python >3.x). Window is slightly more painful. I reccomend Anaconda python.

The mercurial installation is simplest. Make certain Python, Numpy, Scipy, sklearn, netwrokx, graphviz and matplotlib are all installed and correctly working. Then install Mercurial if required.

Go to the directory you wish glbase to install to (for example /Users/name/tools/).

Then at the command type:

{{{$ hg clone https://bitbucket.org/oaxiom/glbase3}}}

This will create a directory 'glbase', containing all of the code.

Next add glbase to your python path by adding this line to ~/.bash_profile (Mac) or ~/.bashrc (Unix, Linux) or equivalent:

{{{export PYTHONPATH=/Users/name/tools:$PYTHONPATH}}}

To update glbase to the most recent version in the repositories, simply:

{{{$ hg pull}}}

and if it prompts you that there has been an update:

{{{$ hg update}}}

To update to the most recent version (i.e. the version I am using).

== Running the test suite ==

These should all pass:

{{{
$ cd glbase/tests
$ python runall.py
}}}

The documentation can be found in:

{{{glbase3/docs/build/html/index.html}}}

== License ==

glbase is distributed under the MIT license:
{{{
    Copyright (C) 2009-2016 Andrew Hutchins
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
    Except as contained in this notice, the name(s) of the above copyright holders shall not be used in advertising or otherwise to promote the sale, use or other dealings in this Software without prior written authorization.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
}}}

### Author & Contributions ###

Andrew P. Hutchins. 
http://www.chrom-lab.org

You can cite glbase here:

Hutchins AP, Jauch R, Dyla M, Miranda-Saavedra D. (2014) glbase: a framework for combining, analyzing and displaying heterogeneous genomic and high-throughput sequencing data. Cell Regeneration, 3:1.

And here is the paper!

http://www.cellregenerationjournal.com/content/3/1/1/abstract

You can contact me at::
    
{{{
andrewh <at> sustc.edu.cn
}}}
    
I've worked at:
* John Innes Centre, 2000-2004.
* Genome Institute of Singapore, 2004-2010.
* Immunology Frontier Research Center, 2010-2013.
* Guangzhou Institutes of Biomedicine and Health, 2013-2015
* Southern University of Science and Technology, 2016-Present

I'm reasonably friendly and only mainly fearsome. I don't have a sense of humour anywhere near as good as the HOMER guys. After all::

    {{{Chuck Norris doesn't use HOMER for DNA motif discovery, he just scans TF binding sites with his eyes }}}
    
glbase relies upon matplotlib, numpy, scipy and sklearn. 

Development was significantly aided by the assistance and bug reports from:

* Ralf Jauch 
* Mateusz Dyla
* Diego Miranda-Saavedra
* Chu Lee Thean