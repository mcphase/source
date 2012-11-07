#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}

unless ($#ARGV >2) 

{print " program range used to delete (comment out) data out of specified range from a file\n";

 print " usage: range [-d] [-u] col min max *.*   \n (col=column,[min,max]..range \n *.* .. filenname\n -d ... option to delete lines from file";

 print " -u option to undo the previous range command on the same file before reapplying the new range.\n";

 exit 0;}

$column=$ARGV[0];shift @ARGV;

if ($column=~/\s*-d/) {$dd=$column;$column=$ARGV[0];shift @ARGV;}

if ($column=~/\s*-u/) {$uu=$column;$column=$ARGV[0];shift @ARGV;}
$ARGV[0]=~s/x/*/g;
$min=eval $ARGV[0];   shift @ARGV;
$ARGV[0]=~s/x/*/g;
$max=eval $ARGV[0];   shift @ARGV;

foreach (@ARGV)
{  @Lines=();
   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   while($line=<Fin>)
   {
      if($uu=~/\s*-u/) { $line =~ s/#://; }   # Removes previously applied comment

      if($line!~/^\s*#/)                      # Line is not a comment
      {
         $line =~ s/D/E/g; @numbers = split(" ",$line);                  # Splits line into numbers
         if ($numbers[$column-1]<$min || $numbers[$column-1]>$max)       # If value in columns is not in limit...
         {
            if ($dd!~/\s*-d/) { $line = "#:".$line; } else { $line=""; } # ... either comment it or delete it.
         }
      }
      push @Lines, $line;                     # Add the line to an array.
   }
   close Fin;

   unless(open (Fout, ">$file")) { die "\n error:unable to write to $file\n"; }

   foreach(@Lines) { print Fout $_; }

   close Fout;

   print ">\n";

}

#\end{verbatim} 
