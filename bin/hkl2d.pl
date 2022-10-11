#! /usr/bin/perl -w



# perl packages to use

use FileHandle;

use PDL;

use PDL::Graphics::PGPLOT;

use PGPLOT;



#get filename from command line if given

my ($file) = @ARGV;

if ($file&&$file=~/-.*h/)

{print "program hkl2d [options] [file] 

 produces neutron intensity graphic and postscript-file hkl.ps

 Options: -n .... plot reflex number n

          -h .... help

 If no file is given the program uses file mcphas.hkl.\n";exit;

}else{print STDERR "#* $0 *\n";}

if ($file&&$file=~/-/)

{$n=-$file;shift @ARGV;($file)=@ARGV;

$xtitle="x";$ytitle="y";

}

else

{# input the labelling of the axes

 print "x-axis title: <T[K]>";$xtitle=($_=<STDIN>);

  $xtitle="T[K]" unless /./;$xtitle=~s/\n//;	

 print "y-axis title: <\\gm\\d0\\uH[T]>";$ytitle=($_=<STDIN>);

  $ytitle='\gm\d0\uH[T]' unless /./;$ytitle=~s/\n//;

}



$file = "./mcphas.hkl" unless $file;





     ($data1,$Limits1) = get_detector_data_2D($file);



    overview_plot("/xserv", $data1,$Limits1); # plot on screen

    overview_plot("./hkl2d.ps/cps", $data1,$Limits1); #then plot on psfile

    print "Wrote postscript file 'hkl2d.ps'\n";





sub overview_plot {

    my ($devspec,$data,$Limits) = @_;

#open page on pgplot server to be able to disply/print

    my $dev = dev("$devspec");

    die "PGOPEN failed!" unless $dev > 0;

#plot data

    plot_array_2D($data,$Limits);

#close page on pgplot server

#	pgclos;

	return ();

}



sub plot_array_2D {

    my ($data,$Limits) = @_;

    my ($x0,$x1,$y0,$y1,$dx,$dy,$title) = @{$Limits};



    my $tr = pdl [ $x0 , $dx, 0, $y0 , 0, $dy ];

    my ($min, $max) = (min($data), max($data));

    if ($min == $max) {

	if($min == 0) {

	    $max = 1;

	} else {

	    $min = 0.9*$min;

	    $max = 0.9*$max;

	}

    }

    my $numcol = 64;

    my $ramp = pdl [[ 0,  1/8,  3/8,  5/8,  7/8,  8/8],

		    [ 0,    0,    0,    1,    1,   .5],

		    [ 0,    0,    1,    1,    0,    0],

		    [.5,    1,    1,    0,    0,    0]];

#     pgpage;

     pgbbuf;  #begin a buffer

     pgsci(1); # set color index

    hold;

    pgvstd;# set standard (default) viewport

    pgswin ($x0,$x1,$y0,$y1); #@{$info->{'Limits'}}; #set window

    pgscir(16,16+$numcol-1); #set color representation

    pgbox("BCNSTI", 0.0, 0.0, "BCNSTI", 0.0, 0.0); #draw labeled frame arount viewport

    ctab $ramp; # set color ramp



# the next 3 lines have been inserted to give contours

#    my $ct = pdl [ $max,$min+($max-$min)*0.25,$min+($max-$min)*0.5,

#                   $min+($max-$min)*0.75, $min ];

#    cont ($data,{CONTOURS=>$ct,TRANSFORM=>$tr});



    imag $data, $min, $max, $tr;

    pgwedg("RI", 0.5, 3.0, $min, $max, ' '); #draw color list on left

    pglab($xtitle, $ytitle, "");

    pgmtxt("T", 2.5, 0.5, 0.5, $title);

#    pgmtxt("T", 1.0, 0.5, 0.5, "[$info->{'Filename'}]");

    pgebuf; #end the buffer and

    release; #display it

}





# Get numerical data, reading it from file 

sub get_detector_data_2D {

    my ($info) = @_;

     unless ($info->{'file read'})

    {($info->{'Numeric Data'},$ytext,$info->{'file read'}) = read_data_file_2D($info);

     # let user select axes

     #print "\n 1. x\n 2. y\n 3. T[K]\n 4. |H|[T]\n 5. Ha[T]\n 6. Hb[T]\n 7. Hc[T]\n";

     #print "select x axis (1-7) ? ";

     # $in=<STDIN>;

     $xcol=0;

     $ycol=1;

     $t=$info->{'Numeric Data'};@dim=$t->dims;

     $x0=min(($info->{'Numeric Data'})->slice("$xcol,:"));         

     $x1=max(($info->{'Numeric Data'})->slice("$xcol,:"));         

     $y0=min(($info->{'Numeric Data'})->slice("$ycol,:"));         

     $y1=max(($info->{'Numeric Data'})->slice("$ycol,:"));         

     # calculate  step intervals

     $dx=0;for ($i=0;$dx==0&&$i<$dim[1]-1;++$i)

              {$dx=abs(-$t->at($xcol,$i)+$t->at($xcol,$i+1));}

     if ($dx==0) {$dx=1;$x1=$x0+1;}

     $dy=0;for ($i=0;$dy==0&&$i<$dim[1]-1;++$i)

              {$dy=abs(-$t->at($ycol,$i)+$t->at($ycol,$i+1));}

     if ($dy==0) {$dy=1;$y1=$y0+1;}

     $nn=int(($x1-$x0)/$dx-0.000001)+2;

     $m=int(($y1-$y0)/$dy-0.000001)+2;

     $b=$t->slice("7,:");

     $t=$b->clump(2);

     reshape $t,$m,$nn;

     $b=$t->mv(0,1);

     $info->{'Numeric Data'}=$b;

     

     ($info->{'Limits'})=[$x0,$x1,$y0,$y1,$dx,$dy,$ytext];

     } 

    return ($info->{'Numeric Data'},$info->{'Limits'});

}





# Read 2D numeric data, skipping comment lines.

sub read_data_file_2D {

    my ($file) = @_;

    my $h = new FileHandle;

    my $l = new FileHandle;

    my @v=();

    my @imax=();



format STDOUT =

@<<. @<<<<<<<<<<<<<<<<<<<<<<<<<<<< imax=@#####.##### 

$nn,$_,$imax[$nn-1]       

.

    if(open($h, $file)) {# let user select reflection

     while(<$h>)

     {next if /^\s*#/;  $x=new PDL(split " ");

      for ($k=7;$k<(($x->dims)[0]-1);$k+=4)

      {$y=$x->slice($k.":".($k+2));

       $ok=0;

       $nn=0;foreach (@v){++$nn;if (sum($y==$_)==3){

         if ($imax[$nn-1]<$x->at($k+3)){$imax[$nn-1]=$x->at($k+3);}

         $ok=1;}}

       if ($ok==0)

       {$v[($#v-$[+1)]=($y);   # addy $y to list $v

        $imax[($#imax-$[+1)]=$x->at($k+3);

       }

      }

     }

     unless ($n)

     {$nn=0;foreach (@v){++$nn;write STDOUT;}

      print "reflection(1-".$nn.")? ";$n=<STDIN>;

     }

     --$n;

     close $h;

     my @xlist = ();

     # input data int piddle

     open($h,$file);

     open($l,">./hkl.asc");

     print $l ("# x y T[K] |H| Ha Hb Hc [T]  vs Intensity of  ".$v[$n]."vs sqrt(int)\n");

     while(<$h>)

     {next if /^\s*#/;  $x=new PDL(split " ");

      $ok=0;

      for ($k=7;$k<(($x->dims)[0]-1);$k+=4)

      {$y=$x->slice($k.":".($k+2));

       if (sum($y==$v[$n])==3){

       #add xyz to piddle

       $y=$x->slice("0:6");

       $z=$x->slice(($k+3).":".($k+3));

       push(@xlist,$y->append($z));

#       print $y->append($z)."\n";

       $out=(($y->append($z*4))->append(sqrt($z)*2));

       $out=~s/\[/ /g;   

       $out=~s/\]/ /g;   

       print $l ($out."\n");

       $ok=1;

                              }

      } 

      if ($ok==0)

      { $y=$x->slice("0:6");

        $z=pdl[0];

	push(@xlist,$y->append($z));

      } 

     }

     close $h; close $l;

     return ((cat @xlist),($v[$n]),1);

    } else {

	print STDOUT "Error program hkl2d: failed to read data file \"$file\"\n";exit 1;

	return undef;

    }

}



