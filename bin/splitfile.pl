#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}


unless ($#ARGV >0) 

{print " program splitf  used to split file  *.* into files named 001*.* 002*.* ...\n";
 print " looking for a given string at the line begin\n";

 print " usage: split string *.*   \n string=string to look for\n *.* .. filenname\n";

 exit 0;}

$string=$ARGV[0];shift @ARGV;

  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "splitting $file at every line beginning with '$string'\n";   
   $i=0;$filename=sprintf("%03i%s",$i,$file);print "$filename ";
   $line=<Fin>;
   open (Fout, ">$filename");
   print Fout $line;
   while($line=<Fin>)
      {
       if ($line=~/^\s*$string/) {close Fout;++$i;
                                  $filename=sprintf("%03i%s",$i,$file);print "$filename ";
                                  open (Fout, ">$filename"); 
                                 }
       print Fout $line;
      }
      close Fin;
      close Fout;
   }
#\end{verbatim} 