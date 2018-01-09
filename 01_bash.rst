Introduction to Command Line Environments
=========================================

The Command Line Interface
--------------------------

A terminal window is a type of "command line interface," a traditional method of interacting with a computer. This is in contrast to modern desktop environments which use "graphical user interfaces" for interaction.

The way to interact with a command line is straightforward:

1. the user (that's you!) enters a (case-sensitive) command with the keyboard
#. the user presses "Enter" on the keyboard
#. the computer processes your request in a series of steps
#. the results of the command (if any) are reported to the screen as text


The Shell
---------

While interacting with the command line feels closer to direct interaction with the computer than a GUI, in fact the instructions you write are not direct computer instructions that are handed directly to your system. Rather, there is an intermediary: a piece of software called a "shell." This is a program that takes command line instructions and executes them for you, by interpreting your commands and dispatching the resulting tasks.

There are several types of shells in common use (much as there are different internet browsers or word processors you might choose). The most common, and the one you're using right now, is called "BASH." Though there are differences between the shells, mostly you can treat them as the same.


Filesystems and Files
=====================

Filesystems
-----------

When you look at a command line interface, you should see a program that's waiting for your input. That program (the shell) has a current "working directory" that corresponds to a location in a computer.

The computers you're working on have a "filesystem" which is typically represented as a collection of folders containing the files you create and download onto your computer, that are stored on your hard drive.

So, for example, if you open the file browser on your computer (Windows Explorer, Finder, or one of various file managers for Linux), you will probably see folders called "Documents," "Pictures," "Downloads," etc., depending on your operating system. These are all folders, which are themselves bundles of files and other folders. Note that whichever folder it's currently showing you is in fact the manager's working directory: typically "/home/cpalmer" or "C:/Users/cpalmer" or equivalent. Clicking on a folder changes the manager's working directory, and it shows you the contents of its new working directory.

The terminal you're using has a working directory: that is, a folder on a computer that is its current "location:"

- when you ask it for files, it will show you files within the working directory
- when you ask it to move to another folder, it will interpret that as relative to the working directory
- when you ask it to write a file to disk, it will do so in the working directory

Files
-----

In addition to folders, if you keep opening subfolders, you'll eventually see files: units of saved text and data that can be manipulated by other software, printed, sent by email, etc. There are two broad categories of files:

1. Plain text (or ASCII/Unicode) files, which contain human-readable characters that can be easily manipulated by the shell
2. Binary files, which store information in non-human-readable formats and require software to parse

The easiest way to determine the type of a file is to look at the "file extension," or series of characters after the final "." in the filename. If you are not on Windows, you might see something like:

- filename.pdf
- results.docx
- data.tar.gz

(In Windows, these extensions are by default hidden from you and represented differently). The extensions in this case are "pdf," "docx," and "tar.gz," referring to Adobe PDF files, Microsoft Office Word files, and a compressed data format called a gzipped archive or "tarball." These file extensions have standard formats and interpretations, and in this case none of them are human readable! If you use the shell to attempt to manually open these files, you will see illegible binary strings. Some plain text file extensions include "csv," "tsv," and "txt."

File extensions are very important to help understand files, but note that they can be changed without altering the format of the file, making them potentially misleading.


Using the shell
===============

In most cases, the shell takes commands, instructions from you, in the following standard format:
   
.. code-block:: bash

		name-of-program [options-for-program] [files-to-modify]

Here, only the name of the program is strictly necessary. Many programs can be run without additional input, but most bioinformatics tools will require (potentially very) long sets of options and filenames. So, to use a real example:

.. code-block:: bash

		pwd

If you type this into your shell and press "Enter," you'll see something get printed out, along the lines of "/home/palmercd" (with your username instead). This is the current working directory of the shell.

In this case, "pwd" is the name of a program (Print the Working Directory). The shell interpreted the first entry on the line to be a piece of software, located it on the system, and ran it with no further arguments. The program generated the string "/home/palmercd" and the shell printed that result out for you.

There are many command line programs that run in a similar fashion. For example:

.. code-block:: bash

		ls

This program, "ls," prints the contents of the current working directory. Each thing that got printed out was the name of a file or folder in the shell's current working directory.

But perhaps you want more information. You might instead enter:

.. code-block:: bash

		ls -l

Now you see very different output. In this case, you ran the same program (ls) as before, but you requested a long listing format.

Arguments to programs come in various formats. You might see things like:

- --flag-name or --flagname
- -f -l -n
- -fln
- --output my_filename.txt
- --strand forward

Different programs will expect input in different formats (and will typically generate grumpy error messages when they encounter inputs they don't understand). The stock programs you find in many terminal environments (like pwd and ls) will often accept flags of a standard format, but that's merely a convention, and bioinformatics tools will ignore these conventions the majority of the time. You will have to look up help documentation to figure out each program's usage. Remember: unless something explicitly says otherwise, assume everything is case-sensitive!

Stopping a Program
------------------

It is extremely common to execute a command on the command line, and then want to stop it before it finishes. The way to do this is CTRL-C (the control key and "c" pressed simultaneously). That sends an interrupt command to the program, and most programs will respond by terminating as soon as they can.

If this doesn't work and you're really enthusiastic about stopping the command, there is also the option to close the terminal you're using.

GETTING HELP
============

There are many ways to get help with command line interfaces.

Man(ual) Pages
--------------

The most common is built right into the shell:

.. code-block:: bash

		man ls

"Man pages" are common help guides that typically follow a standard format: they'll take some getting used to, but they're very common and worth learning how to read. In the majority of cases, this will pull up some very useful documentation for the software you're using. There are, however, some systems that lack man pages, and in that case you'll have to use a different method.

Once you're viewing a man page, you may notice that things look rather different, it's hard to navigate, and you can't escape! When you open a man page, you've secretly opened a type of text viewing software, which displays its contents in the terminal like Word would display text content in its own window. Hopefully: to navigate, use the up/down arrow keys; page up and down are 'b' and 'space' respectively; 'q' is quit.

Internal Help
-------------

Many programs will be capable of generating some help for you themselves:

.. code-block:: bash

		ls --help

Again, compliant software will respond to "--help," "-h," or may potentially emit help documentation when just run without any arguments (if arguments are required for the software to run).

Google
------

If the above methods fail or provide insufficient information, or if you do not know the exact name of the software you want information about, the internet is very likely able to help. There are examples (of varying quality) for most basic and intermediate tasks you need to perform in the command line environment.


Working with Data in the Shell
==============================

With these basics, we can now start working with data!

Setting up a workspace
----------------------

First we need a place to work. You will recall that the shell has a current working directory. That directory will often be the same exact directory each time we open the shell; however, we'd like to have a special location for the analysis we perform right now.

.. code-block:: bash

		mkdir experiment_10jan2018

This is using the program "mkdir" to make a directory for you. You are requesting that the directory be called "experiment_10jan2018" and it will create this directory in the shell's current working directory. How can you confirm that this directory now exists?

Navigating the filesystem
-------------------------

If you want to change your current directory, you can do so as follows:

.. code-block:: bash

		cd experiment_10jan2018

Now, what is your current working directory? What are the results when you list the contents of the current directory?


Absolute versus relative paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the above command, you requested that you change your working directory to a directory called "experiment_10jan2018." You specified a "relative path" : the name of a directory relative to your current working directory.

You can also specify "absolute paths" : the name of a directory relative to the very top of the filesystem. This is how the command pwd reports paths. How exactly this works varies based on the exact system you're using:

	- in Windows, an absolute path will start "C:\\Users\\cpalmer\\"...
	- in Mac OSX, an absolute path will start "/Users/cpalmer/"...
	- in Linux, an absolute path will start "/home/cpalmer/"...

Most commands in the terminal will accept either absolute or relative paths. Custom bioinformatics software has been known to be very picky about the type of path you supply, however, so note what is used in example documentation.

Directory shortcuts
~~~~~~~~~~~~~~~~~~~

Typing out paths can be very cumbersome, especially if you are manually exploring a directory tree with cd. There are some shortcuts that can be used to rapidly specify directories:

	- "./" : refers to the current working directory
	- "../" : refers to the directory one level above the current working directory
	- "~" : refers to the current user's home directory (/home/cpalmer or equivalent)
	- "cd -" : navigate to the previous working directory

What would each of the following commands (derived from http://swcarpentry.github.io/shell-novice) do?

.. code-block:: bash

	cd .
	cd /
	cd /home/cpalmer
	cd ../../
	cd ~
	cd home
	cd ~/Documents/..
	cd
	cd ..

We used ls with no arguments earlier. What happens when you run

.. code-block:: bash

	cd

How can you fix it? Find your way back to experiment_10jan2018/




But why??
+++++++++

If you're wondering why relative paths versus absolute paths matter, or why you would possibly want "./" as a way of fully specifying a path to where you currently are: the shell has a formal method of determining where it should find programs or files you're writing on the command line. BASH secretly searches a series of predetermined directories for each program you use: so for example, when you write

.. code-block:: bash

	ls -l

BASH looks in certain system directories for a program called "ls." It does not usually expect a program called "ls" to exist in your current working directory. If you had such a program, you would have trouble actually running it, but you'd very cleverly write

.. code-block:: bash

	./ls -l

and you would unambiguously tell BASH that you want the "ls" that is in your local working directory, not the one that is in the global search path.


Retrieving data
---------------

Next, we want to get data. You have likely seen in publications that people have posted their data to GEO or some other public repository, or perhaps Ryan has provided files for you to download at some point. We can get files from these public repositories using the command line:

.. code-block:: bash

		wget https://raw.githubusercontent.com/lcdb/genomics-workshop-2018/master/data/Dm.CLAMP.3prRNA.bedGraph.gz

"wget" is a program that downloads a file from the internet according to a URL provided to it. Think of it as the equivalent of "Right Click -> Save Link As" in a web browser. If the download is successful, it should save the file to the current working directory.

Inspecting the data
-------------------

File extensions
~~~~~~~~~~~~~~~

Oftentimes it's not clear what data you're actually working with. There are many ways to inspect files. The first should always be file extension, though remember: it can be changed or left off without actually changing the contents of the file. For example:

.. code-block:: bash

		cp Dm.CLAMP.3prRNA.bedGraph.gz file.pdf

The "cp" command copies the first argument (interpreted as the name of a file that exists) to the second argument (interpreted as the name of the file that will be created). Be careful! This command will destroy any file that currently exists called "file.pdf."

Now, if you inspect the contents of your current directory, you will see another file, but one that is identical to the first. But the file extensions are different!


Removing files
~~~~~~~~~~~~~~

We've just seen the "cp" command for copying one file into a different file. But now we've created a useless file with a misleading filename, so let's make it go away.

Here is the command for removing files. Be careful with this command! Type it carefully, think about what you're doing, and do not use any special characters (SHIFT-NUMBER or brackets) until you know what they do. The remove command is

.. code-block:: bash

		rm file.pdf

What is the state of the current directory?


Data compression and decompression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Real-life bioinformatics datasets tend to be huge, so most files you encounter will be under various forms of compression. Standard types of compression you might encounter are:

- zip files (.zip)
- tar archives (.tar)
- gzip files (.gz, .tar.gz, .tgz)
- bzip2 files (.bz2, .tar.bz2)
- xz files (.xz, .tar.xz)
- rar files (.rar, possibly in pieces)

Choosing which form of compression to use, and deciding how to extract these files, can be complicated. If you are using a standard desktop, I personally recommend a program called "7zip." On the command line, there are programs for extracting data:

- unzip
- untar (or tar x)
- gunzip or gzip -d; tar xzvf
- bunzip2 or bzip2 -d; tar xjvf
- unxz or xz -d; tar xJvf
- unrar (often absent in terminal)


BASH Completion
~~~~~~~~~~~~~~~

Before we move any further, you may notice that typing in the terminal is extremely slow and frustrating, and you may be wondering how anyone does this all day long.

As an example, first list all the files in your current directory:

.. code-block:: bash

		ls -l
		
This has listed everything. If on the other hand we want to just list one thing, we can request that the system provide us suggestions about what to type:

.. code-block:: bash

		ls -l <TAB><TAB>

What has been displayed? How is this different from what happens with the following:

.. code-block:: bash

		ls -l D<TAB>

This is BASH completion: BASH is trying to suggest completions for what you're typing. If you start typing and then <TAB>, it will either suggest all completions (or just fill it in for you, if there's only one candidate). If you haven't typed anything and then tab, it will suggest everything available (which can be tons of things).

It works with filenames, but it also works with software. For example,

.. code-block:: bash

		wg<TAB>

If you don't feel like googling, you can try things like:

.. code-block:: bash

		pdf<TAB>
		
Or for maximum terror:

.. code-block:: bash

		<TAB><TAB>

Note that completion of this type varies depending on your system and context: it usually works but not always. And you'll see more suggested completions in the R section of this tutorial.

So finally: how would you decompress the file you downloaded?

.. code-block:: bash

		gunzip D<TAB>


Making heads or tails of it
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can inspect the beginning or end of a file using the commands

.. code-block:: bash

		head Dm.CLAMP.3prRNA.bedGraph
		tail Dm.CLAMP.3prRNA.bedGraph

What has appeared on your screen?

Determining the size of the file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We saw earlier how to see the size of a file:

.. code-block:: bash

		ls -l

For plain text files, there is another useful command:

.. code-block:: bash

		wc Dm.CLAMP.3prRNA.bedGraph

The program "wc" prints line counts, word counts, and byte counts for a file.

Moving files around, on your computer
-------------------------------------

I'm tired of writing out this horrible filename, even with bash completion! Let's not have to do that anymore.

.. code-block:: bash

		mv Dm.CLAMP.3prRNA.bedGraph my_data.bedGraph

The program "mv" moves a file. Think of this as "cp" followed by "rm." You can use this to send a file to other folders, but in this case we've told it another filename (as with cp), so really we're just renaming the file, in this case to something shorter.

Moving files around, between computers
--------------------------------------

To transfer files between computers, we need to use a different piece of software. There are many types of software for transferring files; everyone has installed "FileZilla" for this purpose.

Lifting over chromosome annotations with UCSC utilities
-------------------------------------------------------

The file you downloaded contains chromosome locations on an old version of the Drosophila genome. Among many other problems, this will interfere with the ability to view the data in the UCSC genome browser.

We can "lift over" the chromosome and physical position data for a file using the command line utility "liftOver." It can be downloaded from the UCSC site (genome.ucsc.edu -> Downloads -> Utilities -> utilities directory), or it is available on the cluster with the following command:

.. code-block:: bash

		module load ucsc
		liftOver

If everything worked, you should see some lines of "help" content, which were emitted by liftOver when we just ran it with no arguments. It is suggesting you need something called "map.chain," or what is termed a "chain file" that links together successive builds of reference genomes.

.. code-block:: bash

		wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz

Note the syntax of the filename! It should be read as "Genome I have to Genome I want."

.. code-block:: bash

		liftOver my_data.bedGraph dm3ToDm6.over.chain.gz output.bedGraph failed.txt

What happened? 

Extracting interesting lines from files
---------------------------------------

.. code-block:: bash

		head my_data.bedGraph
		tail my_data.bedGraph

The file we inspected with head/tail has a header line, followed by an empty line, followed by reasonable looking data content. The error from liftOver suggests that it fails on "field 2" of "line 1," which in this case is "type=bedGraph."


As a first example, say we want all data from "chr2L" within this file. We can do the following:

.. code-block:: bash

		grep -w "chr2L" my_data.bedGraph

What happened? 


Piping output between commands
------------------------------

Far too much content just got emitted to the terminal. We'd like to just see the first few lines of output, to effectively preview the results of our command. To do this, we can "pipe" the output of grep into the head command:

.. code-block:: bash

		grep -w chr2L my_data.bedGraph | head
		
If instead we wanted to know how many results our grep command locates, we can instead pipe the results into wc:

.. code-block:: bash

		grep -w chr2L my_data.bedGraph | wc

What happens if you instead do the following:

.. code-block:: bash

		grep -w chr2L my_data.bedGraph | head | wc


Redirecting output to file
--------------------------

liftOver was complaining about field 2 of line 1 of our bedGraph file. Before we asked for all lines containing a search term. When needed, we can instead ask for all results NOT containing a term:

.. code-block:: bash

		grep -v type=bedGraph my_data.bedGraph | head

The program "grep" responded to our request for all lines not containing "type=bedGraph." Frequently, we want to edit a file using some piece of command line software, and then save the results of that editing to a new file that we can work with later. We can do that with a redirect:

.. code-block:: bash

		grep -v type=bedGraph my_data.bedGraph > my_data_cleaned.bedGraph


What happened? Inspect the contents of your current working directory.

With the new file that lacks the header line, we can now try to run liftOver again:

.. code-block:: bash

		liftOver my_data_cleaned.bedGraph dm3ToDm6.over.chain.gz output.bedGraph failed.txt


Viewing data in a custom track with the UCSC Genome Browser
-----------------------------------------------------------

Now that the file is in the correct genome build, we'd like to visually inspect the genomic distribution of data in the UCSC Genome Browser. But to do that, we need a local copy of the bedGraph file we just created, which we can get with FileZilla.

Finally, navigate to https://genome.ucsc.edu/ and create a custom track.
