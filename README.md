gene-lookup
===========

Gene Lookup for denovo mutations in a set of studies.

Install Instructions
====================

First, clone the repository:

    git clone https://github.com/michaelchess/gene-lookup [directory you want the code to go in]

And cd to the new directory:

    cd [directory you said]

Now, create a new *virtualenv*, which allows us to install the same versions of python packages that Michael used,
isolated away from the rest of your computer.

    virtualenv env

This should create a new directory "env". Now *activate* the virutal environment:

    source env/bin/activate

Your prompt should now start with (env), for example `(env)Brett-Thomas:gene-lookup bt$`

Okay, now time to actually install those packages:

    pip install -r requirements.txt

One of the packages, `rpy`, requires R to be installed.
This is easy if you have a Mac, but first you have to install [Homebrew](http://brew.sh/)
Do that, then:

    brew install R

If that doesn't work, you're probably missing a recursive dependency, the GCC compiler.
This bizarrely isn't preinstalled on mac os. Install XCode through the Mac App Store, then install Command Line Tools
[Here are some instructions](http://docwiki.embarcadero.com/RADStudio/XE4/en/Installing_the_Xcode_Command_Line_Tools_on_a_Mac)

Okay - at this point `pip install -r requirements.txt` should run without error. If it does, you can run the site:

    python GeneLookup.py

And check it out in your browser:

    http://localhost:5000/

Win!

As a side note, this process could be totally skipped if we used a Vagrant virtual machine -
but that poses its own challenges, so I think this is a good middle ground for now. 