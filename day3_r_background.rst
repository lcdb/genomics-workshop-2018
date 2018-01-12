======================
Data Subsetting with R
======================
:Authors: Meeta Mistry, Mary Piper, Radhika Khetani
:Date: Thursday, 11 Jan 2018

* This lesson has been developed by members of the teaching team at the `Harvard Chan Bioinformatics Core (HBC) <http://bioinformatics.sph.harvard.edu/>`_. These are open access materials distributed under the terms of the `Creative Commons Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

* The materials used in this lesson is adapted from work that is Copyright © Data Carpentry (http://datacarpentry.org/). All Data Carpentry instructional material is made available under the `Creative Commons Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ (CC BY 4.0).
  
See the original versions of this content and more at https://hbctraining.github.io/Intro-to-R

Reading data into R
-------------------

Regardless of the specific analysis in R we are performing, we usually need to bring data in for the analysis. The function in R we use will depend on the type of data file we are bringing in (e.g. text, Stata, SPSS, SAS, Excel, etc.) and how the data in that file are separated, or delimited.

When working with genomic data, we often have a metadata file containing information on each sample in our dataset. Let's bring in the metadata file using the `read.table` function. Check the arguments for the function to get an idea of the function options:


.. code-block:: r

		metadata <- read.table(file="data/mouse_exp_design.csv", sep=",")


Inspecting data structures
--------------------------

There are a wide selection of base functions in R that are useful for inspecting your data and summarizing it. Let's use the `metadata` file that we created to test out data inspection functions. We will only show a small number of all that are available, and they are very similar to BASH utilities we've seen.

.. code-block:: r

		head(metadata)
		tail(metadata)
		nrow(metadata)
		ncol(metadata)
		dim(metadata)
		#length(metadata)

Selecting data using indices and sequences
------------------------------------------

When analyzing data, we often want to **partition the data so that we are only working with selected columns or rows.** A data frame or data matrix is simply a collection of vectors combined together. So let's begin with vectors and how to access different elements, and then extend those concepts to dataframes.

Vectors
~~~~~~~

Selecting using indices
+++++++++++++++++++++++

If we want to extract one or several values from a vector, we must provide one or several indices using square brackets [] syntax. The **index represents the element number within a vector** (or the compartment number, if you think of the bucket analogy). R indices start at 1. Programming languages like Fortran, MATLAB, and R start counting at 1, because that's what human beings typically do. Languages in the C family (including C++, Java, Perl, and Python) count from 0 because that's simpler for computers to do.

Let's start by creating a vector called age:

.. code-block:: r

		age <- c(15, 22, 45, 52, 73, 81)


.. image:: images/vector-index.png
   :alt: vector indices


Suppose we only wanted the fifth value of this vector, we would use the following syntax:

.. code-block:: r

		age[5]


If we wanted to select more than one element we would still use the square bracket syntax, but rather than using a single value we would pass in a *vector of several index values*:

.. code-block:: r

		idx <- c(3,5,6) # create vector of the elements of interest
		age[idx]


If we wanted to be incredibly verbose (more on why we would ever want this later), we could also use a vector of logical values to specify what we want from each and every entry in a vector.
To get the third, fifth, and sixth entries of the **age** vector, as above:

.. code-block:: r

		idx <- c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE)
		age[idx]


   **Exercises** 

   1. Create a vector called alphabets with the following alphabets, C, D, X, L, F.
   2. Use the associated indices along with `[]` to do the following:
        * only display C, D and F
	* display all except X
	* display the alphabets in the opposite order (F, L, X, D, C)


Selecting using indices with logical operators
++++++++++++++++++++++++++++++++++++++++++++++

We can also use indices with logical operators. Logical operators include greater than (>), less than (<), and equal to (==). A full list of logical operators in R is displayed below:

+----------+--------------------------+
| Operator | Description              |
+==========+==========================+
| >        | greater than             |
+----------+--------------------------+
| >=       | greater than or equal to |
+----------+--------------------------+
| <        | less than                |
+----------+--------------------------+
| <=       | less than or equal to    |
+----------+--------------------------+
| ==       | equal to                 |
+----------+--------------------------+
| !=       | not equal to             |
+----------+--------------------------+
| &        | and                      |
+----------+--------------------------+
| \|       | or                       |
+----------+--------------------------+

We can use logical expressions to determine whether a particular condition is true or false. For example, let's use our age vector: 
	
.. code-block:: r

		age


If we wanted to know if each element in our age vector is greater than 50, we could write the following expression:	

.. code-block:: r

		age > 50


Returned is a vector of logical values the same length as age with TRUE and FALSE values indicating whether each element in the vector is greater than 50.

.. code-block:: r
		
		[1] FALSE FALSE FALSE  TRUE  TRUE  TRUE


We can use these logical vectors to select only the elements in a vector with TRUE values at the same position or index as in the logical vector.

Create an index with logical operators to select all values in the `age` vector over 50 **or** `age` less than 18:

.. code-block:: r

		idx <- age > 50 | age < 18	
		idx	
		age
		age[idx]




Dataframes
----------

Dataframes (and matrices) have 2 dimensions (rows and columns), so if we want to select some specific data from it we need to specify the "coordinates" we want from it. We use the same square bracket notation but rather than providing a single index, there are *two indices required*. Within the square bracket, **row numbers come first followed by column numbers (and the two are separated by a comma)**. Let's explore the `metadata` dataframe, shown below are the first six samples:

.. image:: images/metadata.png
   :alt: metadata


For example:

.. code-block:: r

		metadata[1, 1]   # element from the first row in the first column of the data frame
		metadata[1, 3]   # element from the first row in the 3rd column


Now if you only wanted to select based on rows, you would provide the index for the rows and leave the columns index blank. The key here is to include the comma, to let R know that you are accessing a 2-dimensional data structure:

.. code-block:: r

		metadata[3, ]    # vector containing all elements in the 3rd row


If you were selecting specific columns from the data frame - the rows are left blank:

.. code-block:: r
		
		metadata[ , 3]    # vector containing all elements in the 3rd column


Just like with vectors, you can select multiple rows and columns at a time. Within the square brackets, you need to provide a vector of the desired values:	

.. code-block:: r

		metadata[ , 1:2] # dataframe containing first two columns
		metadata[c(1,3,6), ] # dataframe containing first, third and sixth rows


For larger datasets, it can be tricky to remember the column number that corresponds to a particular variable. (Is celltype in column 1
or 2? oh, right... they are in column 1). In some cases, the column number for a variable can change if the script you are using adds or removes columns. It's therefore often better to use column names to refer to a particular variable, and it makes your code easier to read and your intentions clearer.

.. code-block:: r

		metadata[1:3 , "celltype"] # elements of the celltype column corresponding to the first three samples


You can do operations on a particular column, by selecting it using the `$` sign. In this case, the entire column is a vector. For instance, to extract all the genotypes from our dataset, we can use: 

.. code-block:: r

		metadata$genotype 

You can use `colnames(metadata)` or `names(metadata)` to remind yourself of the column names. We can then supply index values to select specific values from that vector. For example, if we wanted the genotype information for the first five samples in `metadata`:

.. code-block:: r

		colnames(metadata)
		metadata$genotype[1:5]


The `$` allows you to select a single column by name. To select multiple columns by name, you need to  concatenate a vector of strings that correspond to column names: 

.. code-block:: r

		metadata[, c("genotype", "celltype")]


.. code-block:: r
		
		genotype celltype
		sample1        Wt    typeA
		sample2        Wt    typeA
		sample3        Wt    typeA
		sample4        KO    typeA
		sample5        KO    typeA
		sample6        KO    typeA
		sample7        Wt    typeB
		sample8        Wt    typeB
		sample9        Wt    typeB
		sample10       KO    typeB
		sample11       KO    typeB
		sample12       KO    typeB


While there is no equivalent `$` syntax to select a row by name, you can select specific rows using the row names. To remember the names of the rows, you can use the `rownames()` function:

.. code-block:: r

		rownames(metadata)
		metadata[c("sample10", "sample12"),]


Selecting using indices with logical operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With dataframes, similar to vectors, we can use logical vectors for specific columns in the dataframe to select only the rows in a dataframe with TRUE values at the same position or index as in the logical vector. We can then use the logical vector to return all of the rows in a dataframe where those values are TRUE.

.. code-block:: r

		idx <- metadata$celltype == "typeA"
		metadata[idx, ]



Writing to file 
----------------

Everything we have done so far has only modified the data in R; the files have remained unchanged. Whenever we want to save our datasets to file, we need to use a `write` function in R. 

To write our matrix to file in comma separated format (.csv), we can use the `write.csv` function. There are two required arguments: the variable name of the data structure you are exporting, and the path and filename that you are exporting to. By default the delimiter is set, and columns will be separated by a comma:

.. code-block:: r
		
		write.csv(sub_meta, file="data/subset_meta.csv")


Similar to reading in data, there are a wide variety of functions available allowing you to export data in specific formats. Another commonly used function is `write.table`, which allows you to specify the delimiter you wish to use. This function is commonly used to create tab-delimited files.

*Note*: Sometimes when writing a dataframe with row names to file, the column names will align starting with the row names column. To avoid this, you can include the argument `col.names = NA` when writing to file to ensure all of the column names line up with the correct column values.

Writing a vector of values to file requires a different function than the functions available for writing dataframes. You can use `write()` to save a vector of values to file. For example:


.. code-block:: r

		write(glengths, file="data/genome_lengths.txt", ncolumns=1)



Plotting with ggplot
--------------------

- http://r4ds.had.co.nz/data-visualisation.html (**remainder of the time**)

---

* This lesson has been developed by members of the teaching team at the `Harvard Chan Bioinformatics Core (HBC) <http://bioinformatics.sph.harvard.edu/>`_. These are open access materials distributed under the terms of the `Creative Commons Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

* The materials used in this lesson is adapted from work that is Copyright © Data Carpentry (http://datacarpentry.org/). All Data Carpentry instructional material is made available under the `Creative Commons Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ (CC BY 4.0).
