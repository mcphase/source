#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >0) 

{print " program newline used to create a new line  number n containing text in a data file\n";

 print " usage: newline n text -f *.*   \n n=linennumber \n text= text in new line \n *.* .. filenname(s)\n";

 exit 0;}

$ARGV[0]=~s/x/*/g;$n=eval $ARGV[0];shift @ARGV;
$commandstring=join(" ",@ARGV);
if($commandstring=~/-f/)
{@aa=split('-f',$commandstring);
$text=$aa[0];shift @aa;$files=$aa[0];@bb=split(" ",$files);
}
else
{@bb=split(' ',$commandstring);$text=$bb[0];shift @bb;
}




  foreach (@bb)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   $j=0;
   if ($n==0){print Fout $text."\n";}
   while($line=<Fin>)

     {++$j;
      if ($j==$n){print Fout $text."\n";}
      print Fout $line;

      }
      if($n>$j){print Fout $text."\n";}

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

   }



#\end{verbatim} 

