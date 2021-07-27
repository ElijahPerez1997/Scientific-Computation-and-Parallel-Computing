+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
README file for Homework 7, Math 471, Fall 2020, N.Chang & E.Perez
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 - The homework report is found in the /Report directory, titled HW7.pdf.
 - The code files are located in the sub-directory /Codes
 - The report files are located in the sub-directory /Report
 - To reproduce the results in the report do this:
     1. Log on to Stampede2 supercomputer.
     2. Navigate to HW7/Codes
        Run the makefile.
     3. Submit all job scripts.
     4. For each job, the results are saved in the .out file. The name corresponds to the job name.
     5. Copy and paste the error from the .out file into the y matrix in Parallel_error.m
        Run the MATLAB script. This will generate the error plots that can be found in the current directory.
     6. Repeat the process with the recorded time values, time matrix, and Strong_time.m respectively.
     7. For weak scaling, alter the grid size on line 36 for the desired scaling.
        Submit the job for the desired scaling. Here are the values we used.
            # of Processors         Grid Size
            1                       100
            2                       141
            5                       224
            8                       283
            16                      400
     8. Repeat the process in step (5) with the weak recorded time values, times matrix, and weak_runtime.m respectively.
     9. Navigate to the /Report directory. Move all plots to this directory. Then run the command:
		pdflatex HW7.tex
	You will now see HW7.pdf in the /Report directory



