==============
BASH Refresher
==============

Today we are going to use another brief example of the sort of work you can do with the command line. In this case, we are going to take the most significant ChIP-seq peaks from macs2 and view them in the UCSC genome browser.

What is helix?
--------------

A remote computer system to which you can connect via "ssh," either in the terminal in Mac or in PUTTY in Windows. You log in as a user, just as you do with your personal computer, but the interface is usually text-only.

What is BASH?
-------------

BASH is the "shell" you will use. It is a piece of software, just like Word or Excel or Firefox. You write commands, in the format:

program-name [options] [files]

The BASH interpreter receives these commands, interprets them as best it can, dispatches commands to whatever programs you requested, and reports the result back to you. All command line interfaces are the simple loop: write command, press ENTER, wait as the request is processed, analyze the results.

Navigating in the Shell
-----------------------

The shell has a current working directory, where it writes files, interprets relative paths, etc. You list this directory with

.. code-block:: bash

		pwd

You list the contents of the current working directory with

.. code-block:: bash

		ls -l

You create a new folder inside the current working directory with

.. code-block:: bash

		mkdir new_folder

You navigate into the new folder with

.. code-block:: bash

		cd new_folder

Remember to use <TAB> to autocomplete program names, folder names, and filenames!


Downloading a Remote File
-------------------------

You download a file from the internet using the following:

.. code-block:: bash

		wget https://hpc.nih.gov/~palmercd/tutorial/GSE77625.tar.bz2

Confirm that the file downloaded successfully.

Decompressing Data
------------------

What format is the file you downloaded? You have to use a slightly different command to extract the data:

.. code-block:: bash

		tar xjvf GSE77625.tar.bz2

What has appeared? Navigate into the resulting directory.

Investigating a File
--------------------

What exists in this directory? Make a copy of one of these files with a nicer name:

.. code-block:: bash

		cp GSE77625_h3k27ac_chow.bed my_file.bed

Investigate the contents of this file.

.. code-block:: bash

		head my_file.bed
		tail my_file.bed
		wc my_file.bed

Finding the Number of Signals Per Chromosome
--------------------------------------------

We want to determine how many peaks exist on each chromosome. To do so, we need some new tools.

There are many ways to potentially solve this problem. One involves sorting the file by the chromosome code, and then counting the number of instances of each unique chromosome label.

.. code-block:: bash

		sort my_file.bed | head

What is the result? We want to count the number of instances of each chromosome label in the file.

.. code-block:: bash

		sort my_file.bed | cut -f 1 | head


To get the unique entries, we use the uniq utility:

.. code-block:: bash

		sort my_file.bed | cut -f 1 | uniq

That was efficient, but didn't seem to actually solve the problem. Try the man page!

.. code-block:: bash

		man uniq
		sort my_file.bed | cut -f 1 | uniq -c

Getting the most likely peaks and viewing them in UCSC
------------------------------------------------------

To extract the most interesting peaks, we can sort the file on the fifth column of the file (bigger is better):

.. code-block:: bash

		sort -k 5,5 my_file.bed | head

What happened?

.. code-block:: bash

		sort -k 5,5g my_file.bed | head

This sort appears to be in the reverse order of what we want, which has an easy fix:

.. code-block:: bash

		sort -k 5,5g my_file.bed | tail
		
Now, we select however many of these we want, and write them to file:

.. code-block:: bash

		sort -k 5,5g my_file.bed | tail -20 > my_peaks.bed

Now move the file to your local computer.
