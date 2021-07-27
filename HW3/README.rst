+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
README file for Homework 3, Math 471, Fall 2020, N.Chang & E.Perez
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 - The homework report is found in this subdirectory, titled HW3.pdf.
 - The codes are located in the subdirectory /Codes and are documented in the 
   appendix of the report.
 - The code files are located in the sub-directory /Codes.
 - To reproduce the results in the report do this:
     1. On your machine, navigate to the Codes subdirectory.
        Then run the command:
                make -f makefile (for MAC)
                mingw32-make (for WINDOWS)
     2. The results will be found in Codes/results.csv
     3. Import these results into MATLAB as column vectors with the following names:

	Method from Results		Name
	Trapezoid k = p	                Trap_pi1
	Trapezoid k = pi^2              Trap_pi2
	Gauss k = pi                    Gauss_pi1
	Gauss k = pi^2                  Gauss_pi2

	Run HW3.m in MATLAB.
        Copy the generated plot 'HW3_plot.png' into the Report subdirectory.
     4. Back to the command line, navigate to the Report subdirectory.
        Then run the command:
		pdflatex HW3.tex
     5. You will now see HW3.pdf in the /Report directory

