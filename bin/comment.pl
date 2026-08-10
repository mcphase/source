#!/usr/bin/perl

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program comment  used to comment lines from row1 to row2 with # \n";
 print " usage: comment row1 row2  *.*   \n  *.* .. filenname\n";
 print " alternatively: comment -t string *.*\n";
 print "                      comment all lines containing the text string\n";
 print "                comment -n string *.*\n";
 print "                          comment all lines which do not contain the text string\n";
 print "                comment -cc 8  string *.*\n";
 print "                          insert line containing string before any change of data in column 8\n";
 print "                comment -b string *.*\n";
 print "                          comment all lines before finding a match to string\n";
 print "                comment -a string *.*\n";
 print "                          comment all lines after finding a match to string\n";
 exit 0;}else{print STDERR "#* $0 *";}




$row1=$ARGV[0];shift @ARGV;
$row2=$ARGV[0];shift @ARGV;
if($row1=~/-cc/){$str=$ARGV[0];shift @ARGV;}

@ARGV=map{glob($_)}@ARGV;

      unless($row1=~/-t/||$row1=~/-n/||$row1=~/-a/||$row1=~/-b/)
       { unless($row1=~/-cc/){$row1=~s/x/*/g;$row1=eval $row1;}
         $row2=~s/x/*/g;$row2=eval $row2;
           }

  foreach (@ARGV)

  {$found=0;

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   $i=0;$j=0;

   while($line=<Fin>)

     {unless($line=~/^\s*#/){++$i;}
      if($row1=~/-t/)
       {     if ($line=~/\Q$row2\E/&&!($line=~/^\s*#/)) {print Fout "#".$line;}
             else{print Fout $line;}
       }
      elsif($row1=~/-n/)
       {     if ($line=~/\Q$row2\E/||$line=~/^\s*#/) {print Fout $line;}
             else{print Fout "#".$line;}
       }
      elsif($row1=~/-a/)
       {    if($line=~/^\s*#/){print Fout $line;}
            elsif($found==1){print Fout "#".$line;} else{print Fout $line;}
            if ($line=~/\Q$row2\E/) {$found=1;}
       }
      elsif($row1=~/-b/)
       {    if($line=~/^\s*#/){print Fout $line;}
            elsif($found==0){print Fout "#".$line;} else{print Fout $line;}
            if ($line=~/\Q$row2\E/) {$found=1;}
       }
 
      elsif($row1=~/-cc/)
       {     if ($line=~/^\s*#/) {print Fout $line;}
             else{++$j;
                 @num=split(" ",$line);
                  if($j==1){$nums=$num[$row2-1];}
                  else{unless($nums==$num[$row2-1]){print Fout $str."\n";}
                       $nums=$num[$row2-1];
                      }
                   print Fout $line;
                 }
      
       }
       else
       {      if ($i<$row1||$i>$row2||$line=~/^\s*#/) {print Fout $line;}
              else{print Fout "#".$line;}
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



   print ">";

   }

print "\n";

#\end{verbatim} 