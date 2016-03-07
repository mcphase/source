#!/usr/bin/perl


#\begin{verbatim}


unless ($#ARGV >0) 

{print " program splitf  used to split file  *.* into files named 001*.* 002*.* ...\n";
 print " looking for a given string at the beginning of a line
         or 
         looking for changes of the numbers in a column
  \n";

 print " usage: splitfile -s string *.*
                splitfile -c column *.*
   
           string=string to look for
           column  column to read and look for changes (when value in this column changes
                   a new output file is started)
           *.* .. filenname\n";

 exit 0;}

$option=$ARGV[0];shift @ARGV;
$string=$ARGV[0];shift @ARGV;
@ARGV=map{glob($_)}@ARGV;
  foreach (@ARGV)
  {$file=$_;$set=0;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";} 
   if($option=~/-s/){  
   print "splitting $file at every line beginning with '$string'\n"; }
   else
   {print " splitting $file when column $string changes\n"; } 
   $i=0;$filename=sprintf("%03i%s",$i,$file);print "$filename ";
   $line=<Fin>;
   open (Fout, ">$filename");
   print Fout $line;
   if($option=~/-c/){unless ($line=~/^\s*#/) {$line=~s/D/E/g;@numbers=split(" ",$line);
    $oldvalue=$numbers[$string-1];$set=1;}}
   while($line=<Fin>)
      {
       if($option=~/-s/){if ($line=~/^\s*$string/) {close Fout;++$i;
                                  $filename=sprintf("%03i%s",$i,$file);print "$filename ";
                                  open (Fout, ">$filename"); 
                                 }
                        }
        else
            {if ($line=~/^\s*#/) {print Fout $line;}

            else{$line=~s/D/E/g;@numbers=split(" ",$line);
                 $newvalue=$numbers[$string-1];if($set==0){$oldvalue=$newvalue;$set=1;}
             if ($newvalue!=$oldvalue) {close Fout;++$i;
                                  $filename=sprintf("%03i%s",$i,$file);print "$filename ";
                                  open (Fout, ">$filename"); $oldvalue=$newvalue;
                                 } 
                }
           }
       print Fout $line;
      }
      close Fin;
      close Fout;
   }
#\end{verbatim} 