#!/usr/bin/perl


# program to replace text in several files with another text

unless ($#ARGV >1) 

{print " program substitute used to replace text in files \n"; 

print " usage: substitute oldtext newtext *.* \n"; 

exit 0;} 



$colx=$ARGV[0];shift @ARGV; 

$coly=$ARGV[0];shift @ARGV; 

@ARGV=map{glob($_)}@ARGV;

print "substituting '".$colx."' with '".$coly."' in \n"; 



foreach (@ARGV) 

{ 

$file=$_; 

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

open (Fout, ">range.out"); 

while($line=<Fin>) 

{ 

#       if ($line=~/^\s*#/) 

$line=~s/\Q$colx\E/$coly/g; 

print Fout $line; 

} 

close Fin; 

close Fout; 

unless (rename "range.out",$file) 

{unless(open (Fout, ">$file")) 

{die "\n error:unable to write to $file\n";} 

open (Fin, "range.out"); 

while($line=<Fin>){ print Fout $line;} 

close Fin; 

close Fout; 

system "del range.out"; 

} 



print ">\n"; 



#if (!system ("sed -e 's/".$colx."/".$coly."/' $file > replace.dd"))

#   {system ("mv -f replace.dd $file");}



} 



