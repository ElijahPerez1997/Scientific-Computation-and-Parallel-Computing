+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
README file for Homework 4, Math 471, Fall 2020, E.Perez
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 - The homework report is found in this subdirectory, titled HW4.pdf.
 - The code files are located in the sub-directory /Codes.
 - The report files are located in the sub-directory /Reports
 - To reproduce the results in the report do this:
     1. On your machine, navigate to the HW4/Codes subdirectory.
        Then run the command:
                make -f makefile (for MAC)
                mingw32-make (for WINDOWS)
	Then run the following command to run the program:
		./hw4.exe
	This takes approximately 5 minutes to run
     4. Open and run plotgrid.m and ploterror.m in MATLAB
     5. Navigate to the /Reports directory.
	Then run the command:
		pdflatex HW4.tex
     6. You will now see HW4.pdf in the /Report directory
     7. Note: to produce results for 3 different functions, you will 
	need to manually comment/uncomment function 1, 2 or 3 
	respectively inside xycoord.f90

