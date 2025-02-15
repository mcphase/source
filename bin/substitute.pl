#!/usr/bin/perl


# program to replace text in several files with another text

unless ($#ARGV >1) 

{print " program substitute used to replace text in files \n"; 

print " usage: substitute [option] oldtext newtext *.* 
      
       Option:  -f    only replace the first occurance of oldtext
                -n 21    only replace the 21st occurance of oldtext
\n"; 

exit 0;} else{print STDERR "#* $0 *";}
$first=0;
$_=$ARGV[0];
if(/-f/)
  {$first=1;shift @ARGV; 
  }
if(/-n/)
  {shift @ARGV; $first=$ARGV[0];shift @ARGV; 
  }
$colx=$ARGV[0];shift @ARGV; 

$coly=$ARGV[0];shift @ARGV; 

@ARGV=map{glob($_)}@ARGV;

print "substituting '".$colx."' with '".$coly."' in \n"; 



foreach (@ARGV) 

{ $found=0;

$file=$_; 

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

open (Fout, ">range.out"); 

while($line=<Fin>) 

{ 

#       if ($line=~/^\s*#/) 

if ($first==0)
{$line=~s/\Q$colx\E/$coly/g; 
}else
{if ($line=~/\Q$colx\E/&&$found<$first){++$found;if($found==$first){$line=~s/\Q$colx\E/$coly/;}}
}


print Fout $line; 

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



#if (!system ("sed -e 's/".$colx."/".$coly."/' $file > replace.dd"))

#   {system ("mv -f replace.dd $file");}



} 

print "\n";



