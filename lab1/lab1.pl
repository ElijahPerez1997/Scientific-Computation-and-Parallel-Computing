#!/usr/apps/bin/perl

# Above is to tell the shell that this is a perl script
#
# This simple perl program opens the file lab1.tex.TEMPLATE  
# and performs some text string substitutions and saves it under a new name
# 
#
# run with : perl lab1.perl
#
#
# These are the names of the template and output files. 
$cmdFile="./lab1.tex.Template";
$outFile="./lab1_auto.tex";

$author1 = "Elijah Perez"; 
$author2 = "Natlie Chang";
$date1 = "9/3/20";
$fig1 = "lab1_.png";
$width1 = "width=0.5";


# Open the Template file and the output file. 
open(FILE,"$cmdFile") || die "cannot open file $cmdFile!" ;
open(OUTFILE,"> $outFile") || die "cannot open file $outFile!!" ;
# Setup the output file based on the template
# read one line at a time from the template file
while( $line = <FILE> )
  {
    # Replace the the stings by using substitution
    # s
    $line =~ s/AUTHOR1/$author1/;
    $line =~ s/AUTHOR2/$author2/;
    $line =~ s/today/$date1/;
    # changed the file type of the figure because the eps file would not work. Changed in Matlab using -dpng then changed it here.
    $line =~ s/lab1.eps/$fig1/;
    $line =~ s/width=0.32/$width1/;

    print OUTFILE $line;
    # You can always print to secreen to see what is happening.
    # Uncomment the next line to print the output to the screen.  
    # print $line;
  }
# Close the files
close( OUTFILE );
close( FILE );

system("diff lab1.tex lab1_auto.tex"); 

exit

