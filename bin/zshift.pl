#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program zshift  used to shift coly  by  a constant such that it is zero at a specified value of colx \n";

 print " usage: zshift  constx colx coly  *.*   \n  constx=x-value, colx,coly= columns \n *.* .. filenname\n";

 exit 0;}

 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$constx=eval $ARGV[0];shift @ARGV;

$ARGV[0]=~s/x/*/g;$colx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$column=eval $ARGV[0];shift @ARGV;




  foreach (@ARGV)

  {

   $file=$_;
$delta=1e22;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   while($line=<Fin>)

   {$line=~s/D/E/g;@numbers=split(" ",$line);

        if (($d=abs($numbers[$colx-1]-$constx))<$delta)

	{$delta=$d;

	 $const=-$numbers[$column-1];

       # print $numbers[$colx-1]." $const\n";

        }    

         

   }

 

  close Fin;

   open (Fin, $file);  

   open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

           	  $i=0;++$j;

		  foreach (@numbers)

		  {++$i;

		  if ($i==$column) {$numbers[$i-1]+=$const;}

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



   print "> shifted by $const\n";

   }



#\end{verbatim} 