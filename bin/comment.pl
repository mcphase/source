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
 exit 0;}else{print STDERR "#* $0 *";}




$col1=$ARGV[0];shift @ARGV;
$col2=$ARGV[0];shift @ARGV;
if($col1=~/-cc/){$str=$ARGV[0];shift @ARGV;}

@ARGV=map{glob($_)}@ARGV;

      unless($col1=~/-t/||$col1=~/-n/)
       { unless($col1=~/-cc/){$col1=~s/x/*/g;$col1=eval $col1;}
         $col2=~s/x/*/g;$col2=eval $col2;
           }

  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   $i=0;$j=0;

   while($line=<Fin>)

     {++$i;
      if($col1=~/-t/)
       {     if ($line=~/\Q$col2\E/&&!($line=~/^\s*#/)) {print Fout "#".$line;}
             else{print Fout $line;}
       }
      elsif($col1=~/-n/)
       {     if ($line=~/\Q$col2\E/||$line=~/^\s*#/) {print Fout $line;}
             else{print Fout "#".$line;}
       }
      elsif($col1=~/-cc/)
       {     if ($line=~/^\s*#/) {print Fout $line;}
             else{++$j;
                 @num=split(" ",$line);
                  if($j==1){$nums=$num[$col2-1];}
                  else{unless($nums==$num[$col2-1]){print Fout $str."\n";}
                       $nums=$num[$col2-1];
                      }
                   print Fout $line;
                 }
      
       }
       else
       {      if ($i<$col1||$i>$col2||$line=~/^\s*#/) {print Fout $line;}
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