#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

unless ($#ARGV >0) 

{print " program appendfile  used to append file2 [and file 3 ..] to file1\n";
 print " usage: append file1 file2 [file3 ...] \n";
 exit 0;}

$file1=$ARGV[0];shift @ARGV;
print "appending to $file1:\n";

  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
   open (Fout, ">>$file1");
   while($line=<Fin>)
     {print Fout $line;}
      close Fin;
      close Fout;
   print ">\n";
   }

