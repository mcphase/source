#! /usr/bin/perl



# perl packages to use
use FileHandle;
use PDL;





#get filename from command line if given

if ($#ARGV<1)

{print "program sortf: sorts file according to ascending values in column

    use as sortf column filename

          column  ...... column containing values to be sorted

	  filename... filename to be read 
\n";



exit(0);

}

$ARGV[0]=~s/x/*/g;$col=eval $ARGV[0]; shift @ARGV;--$col;
 foreach (@ARGV)
  { $file=$_;print "<".$file;
    open (Fout, ">range.out");
   my ($data) = read_data_file_2D($_);
   $ix = qsorti $data;
    for($i=0;$i<$ix->getdim(0);++$i)
    {$row=$ix->at($i);
     print Fout $lines[$row];
    }
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

# Read 2D numeric data, skipping comment lines.
sub read_data_file_2D {
    my ($file) = @_;
    my $hh = new FileHandle;
     my @xlist = ();
     @lines= ();
  # input data int piddle
  if(open($hh,$file))
  {      while(<$hh>)
     {if (/^\s*#/){print Fout $_;}
      else
      {$line=$_;
       @numbers=split(" ",$line);      
       $x=new PDL($numbers[$col]);
       push(@xlist,$x);
       push(@lines,$line);
      }
     }
  
   close $hh;
    return (cat @xlist);
   } else {
	print STDOUT "Warning: failed to read data file \"$file\"\n";
	return undef;
   }
}



