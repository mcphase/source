#! /usr/bin/perl



# perl packages to use

use FileHandle;

use PDL;





#get filename from command line if given

if ($#ARGV<2)

{print "program spectrum: creates spectrum as xy table from file mcdisp.qom or mcdisp.dsigma

   use as spectrum h k l [filename(s)] 

          hkl  ...... scattering vector (millerindices)
	  filename... filename to be read (results/mcdisp.qom as default)

   ouput: to stdout - energy vs intensity

\n";
exit(0);
}else{print STDERR "#* $0 *\n";}



#if(($in)=($ARGV[0]=~m|\A(\d)\Z|)){shift @ARGV;}
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$h=eval $ARGV[0]; shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$k=eval $ARGV[0]; shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$l=eval $ARGV[0]; shift @ARGV;

$eps=0.00001;

my (@file) = @ARGV;



 unless ($file[0])

{$file[0] = "./results/mcdisp.qom";}


 

    my ($data,$Limits) = get_detector_data_2D($file[0]);

    my ($intcol,$Ecol) = @{$Limits};

    

    foreach(@file)

    { print "#using data from file ".$_."\n";
      print "#E(meV) Imag[barn/sr/f.u.]\n";
     my ($data,$Limits) = get_detector_data_2D($_);

    $int=$data->slice("$intcol,:");



    for ($kk=0;$kk<(($data->dims)[1]);++$kk)

    {

     if(($data->at(4,$kk))>$h-$eps&&($data->at(4,$kk))<$h+$eps){

     if(($data->at(5,$kk))>$k-$eps&&($data->at(5,$kk))<$k+$eps){

     if(($data->at(6,$kk))>$l-$eps&&($data->at(6,$kk))<$l+$eps){

          print $data->at($Ecol,$kk)." ".$int->at(0,$kk)."\n";

     }}}

    }

}

# Get numerical data, reading it from file 

sub get_detector_data_2D {

    my ($info) = @_;

    {

    

    ($info->{'Numeric Data'}) = read_data_file_2D($info);

     $Ecol=7;

     $intcol=8;

     

     ($info->{'Limits'})=[$intcol,$Ecol];

    

     } 

    return ($info->{'Numeric Data'},$info->{'Limits'});

}





# Read 2D numeric data, skipping comment lines.

sub read_data_file_2D {

    my ($file) = @_;

    my $hh = new FileHandle;

    my @v=();

     my @xlist = ();

  # input data int piddle

  if(open($hh,$file))

  { 

     while(<$hh>)

     {next if /^\s*#/;

      # detect > to see where int data starts [Ha Hb Hc T h k l en1 en2 en3 ... > int1 int2 int3 ...

      if (/^.*\Q>\E/)

       {($a)=(/^(.*)\Q>\E/);($b)=(/^.*\Q>\E(.*)\n/);

       $x=new PDL(split " ",$a);

       if (($x->dims)[0]>8){$xx=new PDL(split " ",$b);}

       else {$xx=new PDL($b,0);}

       }

      else

      {$x=new PDL(split " ");if (($x->dims)[0]>7) {$xx=xvals(($x->dims)[0]-7)+1;}

      } 

      # extract en and int

      for ($kk=7;$kk<(($x->dims)[0]);++$kk)

      {$y=$x->slice("0:6");

       #add xyz to piddle

       $z=$x->slice($kk.":".$kk);

       $int=$xx->slice(($kk-7).":".($kk-7));

       push(@xlist,$y->append($z)->append($int));

       $out=$y->append($z)->append($int);$out=~s/\[/ /g;$out=~s/\]/ /g;   

       }

     }

     

     close $hh; 

     return (cat @xlist);

    } else {

	print STDOUT "Warning: failed to read data file \"$file\"\n";

	return undef;

    }

}



