#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program comment  used to comment lines from row1 to row2 with # \n";

 print " usage: comment col1 col2  *.*   \n  *.* .. filenname\n";

 exit 0;}



$command="#addc";foreach $d(@ARGV){$command.= " ".$d;}; $command.="\n";



$col1=$ARGV[0];shift @ARGV;

$col2=$ARGV[0];shift @ARGV;



  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   $i=0;

   while($line=<Fin>)

     {++$i;

       if ($i<$col1||$i>$col2||$line=~/^\s*#/) {print Fout $line;}

       else{print Fout "#".$line;}

      }

      close Fin;

 #     print Fout $command;

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

   }



#\end{verbatim} 