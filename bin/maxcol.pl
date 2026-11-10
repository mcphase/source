#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

unless ($#ARGV >0) 

{print STDERR << "EOF";
# maxcol - too few arguments
# Prgram maxcol is used to find local maxima in data file
#
# Usage: maxcol col *.*
#
#  col    ...   column for which the values should be maximal
#  *.*    ...  filename
# 
#  the program scans file *.* line by line and if a maximum is found 
#  in column col, the line is output to stdout
# 
#  example:
#  maxcol 3 file.dat > maximums.clc
#
EOF
exit;
}
$ARGV[0]=~s/x/*/g;$column=eval $ARGV[0];shift @ARGV;


  foreach (@ARGV)

  {
$file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "# maxima in column $column in  file <".$file.">\n";
$i=0;
while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print $line;}

       else{
       ++$i;$line=~s/D/E/g;@numbers=split(" ",$line);
      $veroldvalue=$oldvalue;$oldvalue=$value; 
      $value=$numbers[$column-1];
      if($i>2){if($oldvalue>$value&&$oldvalue>$veryoldvalue) {print $storedline;}}
      $storedline=$line;        
           }
  }
close Fin;



}