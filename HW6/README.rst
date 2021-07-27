+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
README file for Homework 6, Math 471, Fall 2020, E.Perez
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 - The homework report is found in the /Report directory, titled HW6.pdf.
 - The code files are located in the sub-directory /Codes
 - The report files are located in the sub-directory /Report
 - To reproduce the results in the report do this:
     1. Execute the following on stampede2 supercomputer:
		hw4.job
		hw6weak.job
		hw6serweak.job
     2. For hw4.job, the results are saved in the following files:
		error20.txt
		error800.txt
		time20.txt
		time800.txt
		s_error20.txt
		s_error800.txt
		s_time20.txt
		s_time800.txt
     3. For hw4weak.exe, the results are saved in the following files:
		weakresults.txt
     4. For hw6serweak.job, the results are saved in the following files:
		serweakresults.txt
     5. Now that you have the results, reproduce the plots by:
	a) Copy and paste grid20.m, grid800.m and all result files from hw4.job into the same directory. 
	b) Run grid20.m and grid800.m in matlab. This will produce 4 figures in the current directory.
	c) Copy and paste the time and error from weakresults.txt into time and error matrices in weak.m.
	d) Copy and paste the time and error from serweakresults.txt into s_time and s_error in weak.m.
	e) Run weak.m. This will produce two plots in the current directory.
     5. Navigate to the /Report directory. Move all plots to this directory. Then run the command:
		pdflatex HW6.tex
	You will now see HW6.pdf in the /Report directory



