#! /usr/bin/perl -w



# program to plot the phasediagram



# use some perl modules

use FileHandle;

use PDL;

use PDL::Graphics::PGPLOT;

use PGPLOT;



# get filename from command line

my ($file) = @ARGV;

$file = "./mcphas.xyt" unless $file;





print "reading file <".$file.">\n";

print "info: set character hight \\f set font \\g greek \\d subscript \\u normal again\n";



# input the labelling of the axes

print "plot title: <Phasediagram>";$title=($_=<STDIN>);

 $title="Phasediagram" unless /./;$title=~s/\n//;	

print "x-axis title: <T(K)>";$xtitle=($_=<STDIN>);

 $xtitle="T(K)" unless /./;$xtitle=~s/\n//;	

print "y-axis title: <\\gm\\d0\\uH(T)>";$ytitle=($_=<STDIN>);

 $ytitle='\gm\d0\uH(T)' unless /./;$ytitle=~s/\n//;





print <<END;

Click plot for displaying magnetic structure figures.

Drag and drop colors to change colors - (result stored in phased.xyt).

Type 'U' (in graphics window) to undo drag and drop

Type 'P' (in graphics window) for hardcopy, 'Q' or right mouse button to quit.

END



system("cp -f ".$file." phased.xyt");

    $file="phased.xyt";

    

for(;;) {

    my ($cc,$cx,$cy,$idx);    

    ($bx,$by,$cx,$cy,$cc) = overview_plot("/xserv", $file, 1);

    last if $cc =~ /[xq]/i;	# Quit?

    if($cc =~ /p/i) {		# Hardcopy?

	overview_plot("./phased.ps/cps", $file, 0);

	print "Wrote postscript file 'phased.ps'\n";

	next;

                     }

    if($cc =~ /u/i) {		# Hardcopy?

         system "mv -f phased.xyt.old ".$file;

	next;

                     }



    # show magnetic structure ?

    ($T,$ha,$hb,$hc,$id1,$id2)=read_data_th($file,$bx,$by);

    ($T1,$ha1,$hb1,$hc1,$dd,$dd1)=read_data_th($file,$cx,$cy);

    if ($T==$T1&&$ha==$ha1&&$hb==$hb1&&$hc==$hc1)

    {system "spins ".$T." ".$ha." ".$hb." ".$hc;

     system "gv spins3dab.eps &";

     system "gv spins3dac.eps &";

     system "gv spins3dbc.eps &";

     system "gv spinsab.eps &";

     system "gv spinsac.eps &";

     system "gv spinsbc.eps &";

    }

    else    

    # drag and drop ?

    {drag_drop($file,$id1,$id2,$dd,$dd1);

    }

}

system "rm -f phased.xyt.old";





sub overview_plot {

    my ($devspec, $datalist, $interactive) = @_;

    my $dev = dev("$devspec");

    die "PGOPEN failed!" unless $dev > 0;

    plot_array_2D($datalist);

    if($interactive) {

	# Wait for user to select two points.

	my ($ax,$ay,$bx,$by,$cx,$cy,$cc) = (0,0,0,0,0,0,"");

	pgband(0, 0, $ax, $ay, $bx, $by, $cc);

	pgband(1, 0, $bx, $by, $cx, $cy, $cc);

#	pgclos;

	return ($bx,$by,$cx,$cy,$cc);

    } else {

#	pgclos;

close;

	return ();

    }

}





sub plot_array_2D {

    my ($info) = @_;

    my ($data,$Limits) = get_detector_data_2D($info);

    my ($x0,$x1,$y0,$y1,$dx,$dy) = @{$Limits};

#    my $tr = pdl [ $x0 + $dx/2, $dx, 0, $y0 + $dy/2, 0, $dy ];

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

#     pgpage;

     pgbbuf;  #begin a buffer

     pgsci(1); # set color index

    hold;

    pgvstd;# set standard (default) viewport

    pgswin ($x0,$x1,$y0,$y1); #@{$info->{'Limits'}}; #set window

    pgscir(16,16+$numcol-1); #set color representation

    pgbox("BCNSTI", 0.0, 0.0, "BCNSTI", 0.0, 0.0); #draw labeled frame arount viewport



    my $ramp = pdl [[ 0,  1/8,  3/8,  5/8,  7/8,  8/8],

		    [ 0,    0,    0,    1,    1,   .5],

		    [ 0,    0,    1,    1,    0,    0],

		    [.5,    1,    1,    0,    0,    0]];

    ctab($ramp); # set color ramp



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



# Get numerical data for 2D detector, reading it from file if necessary

# the first time.

sub get_detector_data_2D {

    my ($info) = @_;

#     unless ($info->{'Numeric Data'})

    {($info->{'Numeric Data'}) = read_data_file_2D($info);

     # calculate limits

     

     #1. take temperature slice (column 0)

     $xcol=0;

     #2. take field slice (column 1)

     $ycol=1;

     

     $t=$info->{'Numeric Data'};@dim=$t->dims;

     $x0=min(($info->{'Numeric Data'})->slice("$xcol,:"));         

     $x1=max(($info->{'Numeric Data'})->slice("$xcol,:"));         

     $y0=min(($info->{'Numeric Data'})->slice("$ycol,:"));         

     $y1=max(($info->{'Numeric Data'})->slice("$ycol,:"));         

     # calculate  step intervals

     $dx=0;for ($i=0;

     $dx==0&&$i<$dim[1]-1;++$i)

              {$dx=abs(-$t->at($xcol,$i)+$t->at($xcol,$i+1));}

     if ($dx==0) {$dx=1;$x1=$x0+1;}

     $dy=0;for ($i=0;$dy==0&&$i<$dim[1]-1;++$i)

              {$dy=abs(-$t->at($ycol,$i)+$t->at($ycol,$i+1));}

     if ($dy==0) {$dy=1;$y1=$y0+1;}

     $n=int(($x1-$x0)/$dx-0.000001)+2;

     $m=int(($y1-$y0)/$dy-0.000001)+2;

     $b=$t->slice("7,:");

     $t=$b->clump(2);

     reshape $t,$m,$n;

     $b=$t->mv(0,1);

     $info->{'Numeric Data'}=$b;

     

     ($info->{'Limits'})=[$x0,$x1,$y0,$y1,$dx,$dy];

     } 

    return ($info->{'Numeric Data'},$info->{'Limits'});

}





# Read 2D numeric data, skipping comment lines.

sub read_data_file_2D {

    my ($file) = @_;

    my $h = new FileHandle;

    if(open($h, $file)) {

	my @list = ();

	while(<$h>) {

	    next if /^\s*#/;

            $_=~s|p||g;

	    

	    push(@list, new PDL (split " "));

	}

	close $h;

        return (cat @list);

    } else {

	print STDOUT "Warning: failed to read data file \"$file\"\n";

	return undef;

    }

}



# Read 2D numeric data, skipping comment lines.

sub read_data_th {

    my ($file,$cx,$cy) = @_;

    my $h = new FileHandle;

    

    my ($T,$ha,$hb,$hc,$id1,$id2,$dist)=(0,0,0,0,10000);

    $dist=10000;

    if (open($h, $file)) {

            while(<$h>) {

	    next if /^\s*#/;

            $_=~s|p||g;

	    @numbers=split(" ");

           	  $d=sqrt(($numbers[0]-$cx)*($numbers[0]-$cx)+($numbers[1]-$cy)* ($numbers[1]-$cy));

		  if ($d<$dist)

		  {$dist=$d;

		   $T=$numbers[2];

		   $ha=$numbers[4];

		   $hb=$numbers[5];

		   $hc=$numbers[6];		    

		   $id1=$numbers[7];

		   $id2=$numbers[8];

		  }

	}

	close $h;

        return ($T,$ha,$hb,$hc,$id1,$id2);

    } else {

	print STDOUT "Warning: failed to read data file \"$file\"\n";

	return undef;

    }

}



# Read 2D numeric data, skipping comment lines.

sub drag_drop {

    my ($file,$id1,$id2,$d1,$d2) = @_;

    my $h = new FileHandle;

    my $k = new FileHandle;

    

    if (open($h, $file)) {open($k, ">phased.xyt.new");

            while(<$h>) {

	     if (/^\s*#/){print $k ($_);next;}

            $_=~s|p||g;

            @numbers=split(" ");

	    if ($numbers[7]==$d1) {$numbers[7]=$id1;

	      if ($#numbers>10){ $numbers[8]=$id2;} # this line is forthe new file format

	     }

                  $i=0;

		  foreach (@numbers)

		  {++$i;print $k $numbers[$i-1]." ";}     



	     print $k "\n";

	}

	close $h;close $k;

        system "cp -f ".$file." phased.xyt.old";

        system "mv -f phased.xyt.new ".$file;

	return;

    } else {

	print STDOUT "Warning: failed to read data file \"$file\"\n";

	return undef;

    }

}

