#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program factcolerr  used to multiply a  column and errcol with a constant\n";

 print " usage: factcolerr col errcol const  *.*   \n col=column,colerr=column with errorbars const=constant \n *.* .. filenname\n";

 exit 0;}

 

$column=$ARGV[0];shift @ARGV;
$columnerr=$ARGV[0];shift @ARGV;
if ($column==$columnerr){die "Error factcolerr: error column $columnerr = datacolumn $column\n";}
$const=$ARGV[0];shift @ARGV;



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

		  if ($i==$column) {$numbers[$i-1]*=$const;}
		  if ($i==$columnerr) {$numbers[$i-1]*=abs($const);}

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