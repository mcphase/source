#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >0)

{print " program sigcol  used to calculate the signum of  a  column \n";

 print " usage: sigcol col   *.*   \n col=column \n *.* .. filenname\n";

 exit 0;}

 

$ARGV[0]=~s/x/*/g;$column=eval $ARGV[0];shift @ARGV;





  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

           	  $i=0;++$j;

		  foreach (@numbers)

		  {++$i;

		  if (($i==$column)&&($numbers[$i-1]>0)){$numbers[$i-1]=1;}
		  if (($i==$column)&&($numbers[$i-1]<0)){$numbers[$i-1]=-1;}
                  if (($i==$column)&&($numbers[$i-1]==0)){$numbers[$i-1]=0;}

		  print Fout $numbers[$i-1]." ";}     

            print Fout "\n";

           }

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

 }   



#\end{verbatim} 