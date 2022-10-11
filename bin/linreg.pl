#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}



# perl packages to use
use FileHandle;
use PDL;





#get filename from command line if given

if ($#ARGV<1)

{print STDERR "program linreg: calculates linear regression of n columns in file

    use as linreg col n filename

          col     ...... column containing y_k values followed by
          n       ...... |n| columns containing x_ik (i=1 to n)
	  filename... filename

    the program calculates the linear regression, i.e. the best values
    of coefficients ai such that y_k~sum_i a_i*x_ik for every data line k
    in the file. The n linear regression equations solved to determin a_i
    are (i,j=1 ...n):

     sum_k x_jk y_k = sum_i a_i (sum_k x_ik * x_jk)

   Output: - sdtoud: best coefficients a_i  and standard deviation
             sta=sum_k (y_k-sum_i a_i*x_ik)^2
           - environment variable MCPHASE_STA: standard deviation sta
           - file: new column col+n+1 contining sum_i a_i*x_ik
\n";



exit(0);

}else{print STDERR "#* $0 *\n";}

$ARGV[0]=~s/x/*/g;$col=eval $ARGV[0]; shift @ARGV;--$col;
$ARGV[0]=~s/x/*/g;$n=eval $ARGV[0]; shift @ARGV;

 foreach (@ARGV)
  { $file=$_;print "echo \"<".$file."\"\n";
#     sum_k x_jk y_k = sum_i a_i (sum_k x_ik * x_jk)
#     rewritten as b_j= sum_i a_i c_ij
# read data file 2D reads b_j and c_ij as PDL piddles
   my ($b,$c) = read_data_file_2D($_);
# here determine the best coefficients ai in PDL $a:
     $a=($b x inv($c))->clump(-1) ;
# here do the printout to the file and screen
    open (Fout, ">range.out");
  if(open(Fin,$file))
  {      while($line=<Fin>)
     {if ($line=~/^\s*#/){print Fout $line;}
      else
      {$line=~s/D/E/g;@numbers=split(" ",$line);
	  $x=new PDL(@numbers);
        #  print $x;print $a;print $line;
          $x=$x->slice(($col+1).":".($col+$n));
          $ycalc=sum($x*$a);
          $sta+=($ycalc-$numbers[$col])*($ycalc-$numbers[$col]);
          $i=0;
		  foreach (@numbers)
		  {print Fout $numbers[$i]." ";
                   if ($i==$col+$n) {print Fout $ycalc." ";} # print out results
                   ++$i;
 		  }
            print Fout "\n";++$i;
      }
     }  
   close Fin;
  }
  print Fout "# Result of linear regression col".($col+1)." ~ sum_i=".($col+2)."...".($col+$n+1)." ai coli \n";
  for($i=1;$i<=$n;++$i)
 { print Fout "#! a".($col+$i+1)."=".$a->at($i-1)."\n"; }
  print Fout "#! sta=$sta\n";

  print  "echo \"# Result of linear regression col".($col+1)." ~ sum_i=".($col+2)."...".($col+$n+1)." ai coli \"\n";
  for($i=1;$i<=$n;++$i)
 { print "echo \"#! a".($col+$i+1)."=".$a->at($i-1)."\"\n"; }
  print "echo \"#! sta=$sta\"\n";
    
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
      print "echo \">\"\n";
  }

# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_STA=$sta\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_STA=$sta\n";
#close Fout;
 
if ($^O=~/MSWin/){
print "set MCPHASE_STA=$sta\n";
                 }
                 else
                  {
print "export MCPHASE_STA=$sta\n";
                  }


# Read 2D numeric data, skipping comment lines.
#     sum_k x_jk y_k = sum_i a_i (sum_k x_ik * x_jk)
#     rewritten as b_j= sum_i a_i c_ij
# read data file 2D reads b_j and c_ij as PDL piddles
sub read_data_file_2D {
    my ($file) = @_;
    # initialize $b
    $b=zeroes($n);
    $c=zeroes($n,$n);
    my $hh = new FileHandle;
#    my @xlist = ();
  # input data int piddle
  if(open($hh,$file))
  {      while(<$hh>)
     {if (/^\s*#/){}
      else
      {$x=new PDL(split " ");
       $y=$x->slice("$col:$col");
       $x=$x->slice(($col+1).":".($col+$n));
     # print $x; print transpose($x)*$x;
       $b+=$x *$y;
       $c+=transpose($x)*$x;#print $c;
      }
     } 
  
   close $hh;
#   return (cat @blist,cat @clist);
   return ($b,$c);
   } else {
	print STDOUT "Warning: failed to read data file \"$file\"\n";
	return undef;
   }
}



