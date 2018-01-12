====================
Simple BASH Commands
====================

BASH Basics
-----------

Write a comment, the contents of which are ignored:

.. code-block:: bash

        # nonsense goes here

Create an empty file, or update the timestamp of an existing file without changing the contents:

.. code-block:: bash

        touch [filename]

Print the date:

.. code-block:: bash

        date

Get help about a program (typically only works for core system software):

.. code-block:: bash

        man [program-name]

Display a line of text in the terminal:

.. code-block:: bash

        echo 'Something literal I want      printed'
        echo Something I would like the shell to      parse
        echo 'My $PWD is' $PWD

File Management
---------------

Copy an existing file to a new filename, preserving the original:

.. code-block:: bash

        cp [filename] [new.filename]

Make a new copy of a folder and all its contents into a target directory:

.. code-block:: bash

        cp -R [foldername] [target.directory]

Move an existing file to a new filename, removing the original (rename):

.. code-block:: bash

        mv [filename] [new.filename]

Remove a file (DANGEROUS):

.. code-block:: bash

        rm [filename]



Folder Navigation
-----------------

Print the current working directory of the shell:

.. code-block:: bash

        pwd

List the contents of the current working directory (optionally in long form):

.. code-block:: bash

        ls -l

Create a new folder in the current working directory called "dir":

.. code-block:: bash

        mkdir dir

Change directory to a subdirectory of the current working directory called "dir":

.. code-block:: bash

        cd dir

Change directory to the parent directory (choose any one):

.. code-block:: bash

        cd ../
        cd -

Remove a folder in the current working directory called "dir" (only if the directory is empty):

.. code-block:: bash

        rmdir dir

Remove a folder in the current working directory called "dir" (DANGEROUS: will destroy everything in the subdirectory):

.. code-block:: bash

        rm -R dir


Inspecting Files
----------------

List all the files in a directory:

.. code-block:: bash

        ls -l
        ls <TAB><TAB>

Inspect the top X lines of a file:

.. code-block:: bash

        head -X [filename]

Inspect the bottom X lines of a file:

.. code-block:: bash

        tail -X [filename]

Count the number of lines, words, characters in a file:

.. code-block:: bash

        wc [filename]
        wc -l [filename] #line count only
        wc -w [filename] #word count only
        wc -m [filename] #character count only


Downloading Files
-----------------
Download a URL, saving it in the working directory in file named after the last part of the URL:

.. code-block:: bash

    wget http://example.com/file.zip

Download a URL, redirecting its output to a file that you specify:

.. code-block:: bash

    wget -O - http://example.com/file.zip > my_file.zip


Download a URL with special characters in it that would confuse Bash, redirecting it to a file that you specify:

.. code-block:: bash

    wget -O - ""https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77625&format=file&file=GSE7%2Egz" > results.gz

Data Compression and Extraction
-------------------------------

Compress a file using the gzip algorithm, removing the original file and creating one with an added ".gz" suffix:

.. code-block:: bash

        gzip [filename]
        #or bzip2, ".bz2"
        #or xz, ".xz"

Extract a file using the gzip algorithm, removing the original file and creating one lacking the ".gz" suffix:

.. code-block:: bash

        gunzip [filename.gz]
        #or bunzip2, ".bz2"
        #or unxz, ".xz"

Compress a file using the gzip algorithm using maximum compression; redirect the output to a new file, maintaining the original file:

.. code-block:: bash

        gzip -9c [filename] > [out.filename.gz]
        #or bzip2, ".bz2"
        #or xz, ".bz2"

Extract a file using the gzip algorithm; redirect the output to a new file, maintaining the original file:

.. code-block:: bash

        gunzip -c [filename.gz] > [out.filename]
        #or bunzip2, ".bz2"
        #or unxz, ".xz"


Zip a file, maintaining the original file:

.. code-block:: bash

        zip -r [out.filename.zip] [filename]

Unzip a file, maintaining the original file:

.. code-block:: bash

        unzip [filename.zip]


Extract a tar archive, possibly with compression (tarball):

.. code-block:: bash

        tar xvf [filename.tar]
        tar xzvf [filename.tgz or filename.tar.gz]
        tar xjvf [filename.tbz2 or filename.tar.bz2]
        tar xJvf [filename.txz or filename.tar.xz]


Create a tar archive, possibly with compression (tarball):

.. code-block:: bash

        tar cvf [out.filename.tar] [in.targets]
        tar czvf [out.filename.tar.gz] [in.targets]
        tar cjvf [out.filename.tar.bz2] [in.targets]
        tar cJvf [out.filename.tar.xz] [in.targets]

Basic Data Management
---------------------

Sort a file:

.. code-block:: bash

        sort [filename] > [out.filename]

Sort a file, first by the first entry in each line, then by the second entry in each line which is a number:

.. code-block:: bash

        sort -k 1,1 -k 2,2g [filename] > [out.filename]

Get the unique rows of a sorted file:

.. code-block:: bash

        uniq [filename] > [out.filename]

Concatenate as many files as you want, end-to-end, in order, and write the result to a single output file:

.. code-block:: bash

        cat [first.filename] [[second.filename] ...] > [out.filename]

Combine files, line by line, in order:

.. code-block:: bash

        paste [first.filename] [[second.filename] ...] > [out.filename]

Combine two sorted files, based on the value in column X of file 1 and column Y of file 2 (equivalent to excel vlookup):

.. code-block:: bash

        join -1 X -2 Y [first.filename] [second.filename] > [out.filename]

Searching Files for Patterns
----------------------------

Get all lines containing a given character string:

.. code-block:: bash

        grep "casE SensitivE STRing" [filename]

Get all lines containing a given character string, regardless of case:

.. code-block:: bash

        grep -i "casE INSensitivE STRing" [filename]

Get all lines NOT containing a given character string:

.. code-block:: bash

        grep -v "thing to exclude" [filename]

Get all lines containing a given character string, along with the two lines preceding and three lines trailing each match:

.. code-block:: bash

        grep -B 2 -A 3 "rs9939609" [filename]




Search and Replace
------------------

Replace all instances of the word "Goofy" with the string "Very Serious":

.. code-block:: bash

        sed 's/Goofy/Very Serious/g' [filename] > [out.filename]

Replace the first instance of a tab on each line with a single space:

.. code-block:: bash

        sed 's/\t/ /' [filename] > [out.filename]


Replace all instances of the string "/home/cpalmer/Documents/stuff" with the string "/home/hacker/other/garbage.virus" while overwriting the original file:

.. code-block:: bash

        sed -i 's/\/home\/cpalmer\/Documents\/stuff/\/home\/hacker\/other\/garbage.virus/g' [filename]


Extract Individual Columns from a File
--------------------------------------

Extract the fifth row from each column in a file:

.. code-block:: bash

        awk '{print $5}' [filename] > [out.filename]

From a file, extract the third and sixth columns of each row that contains the string "PASS":

.. code-block:: bash

        awk '/PASS/ {print $3" "$6}' [filename] > [out.filename]



Pipes and Redirects
-------------------

Pipe the output of one command into the input of a second, chaining together the commands into one line:

.. code-block:: bash

        awk '/chr2L/' [filename] | sort -k 2,2g > [out.filename]

Chain together many commands into a single call, avoiding the need to write many intermediate files:

.. code-block:: bash

        grep -i "chr2L" [filename] | sort -k 2,2g | awk '! /track/' | sort | uniq | wc -l
