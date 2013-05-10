#! /usr/bin/perl -w
 
 # perl packages to use
 use FileHandle;
 use PDL;
 #use PDL::Graphics::PGPLOT;
 #use PGPLOT;
 
 #get filename from command line if given
 my ($file) = @ARGV;
 if ($file&&$file=~/-.*h/) # look if user wants help message (option -h)
 {print "program hkl [options] [file] 
  produces neutron intensity for one reflection from results/mcphas*.hkl file
  Options: -n .... plot reflex number n
           -h .... help
  If no file is given the program uses file results/mcphas.hkl.\n";exit;
 }
 $xcol=0;$ycol=7;$n=0;  # set default axes
 if ($file&&$file=~/-/) # look if user gives reflection number (option -n)
 {$n=-$file;shift @ARGV;($file)=@ARGV;}
 
 $file = "./results/mcphas.hkl"  unless $file;
 
 
 
 
     overview_plot("/xserv", $file); # plot on screen
 #    overview_plot("./hkl.ps/cps", $file); #then plot on psfile
     print "Wrote output file 'results/hkl.asc'\n";
 
 
 sub overview_plot {
     my ($devspec, $datalist) = @_;
 
     my ($data,$Limits) = get_detector_data_2D($datalist);
 #open page on pgplot server to be able to disply/print
 #    my $dev = dev("$devspec");
 #    die "PGOPEN failed!" unless $dev > 0;
 #unless ($title)
 #{# input plot title
 #print "plot title: <Neutron Intensity>";$title=($_=<STDIN>);
 #   $title="Neutron Intensity"  unless /./;$title=~s/\n//;	
 #}
 
 #plot data
 #    plot_array_2D($data,$Limits);
 #close page on pgplot server
 #	pgclos;
 #close;
 	return ();
 }
 
 sub plot_array_2D {
     my ($data,$Limits) = @_;
     my ($x0,$x1,$y0,$y1,$xcol,$ycol,$xtext,$ytext) = @{$Limits};
 #    pgvstd; # set standard (default) viewport
 #    pgpap(6.5,1.0);
 #    hold;   # wait with display until release command is given
 #    pgswin ($x0,$x1,$y0,$y1); #@{$info->{'Limits'}}; #set window
 #    pgbox("BCNSTI", 0.0, 0.0, "BCNSTI", 0.0, 0.0); #draw labeled frame arount viewport
 
     #draw datapoints from slice of piddle "$data"
 #    points $data->slice("$xcol,:"),$data->slice("$ycol,:") ;
 #    pglab($xtext,$ytext, "");
     print "#".$xtext." ".$ytext."\n";
     print $data->slice("$xcol,:")." ".$data->slice("$ycol,:")."\n" ;
     #label plot
 
 
 #    pgmtxt("T", 2.5, 0.5, 0.5, $title);
 #    pgebuf; #end the buffer and
 #    release; #display it
 }
 
 # Get numerical data, reading it from file 
 sub get_detector_data_2D {
     my ($info) = @_;
      unless ($info->{'file read'})
     {($info->{'Numeric Data'},$ytext,$info->{'file read'}) = read_data_file_2D($info);
      $in=1;
 #     unless ($n)
      unless (1)
      {# let user select axes
      print "\n 1. x\n 2. y\n 3. T(K)\n 4. |H|(T)\n 5. Ha(T)\n 6. Hb(T)\n 7. Hc(T)\n";
      print "select x axis (1-7) ? ";
       $in=<STDIN>;
      $xcol=$in-1;}
      $ycol=7;
      if ($in==1) {$xtext="x";}
      if ($in==2) {$xtext="y";}
      if ($in==3) {$xtext="T(K)";}
      if ($in==4) {$xtext='\gm\d0\u|H|(T)';}
      if ($in==5) {$xtext='\gm\d0\uHa(T)';}
      if ($in==6) {$xtext='\gm\d0\uHb(T)';}
      if ($in==7) {$xtext='\gm\d0\uHc(T)';}
      
      # calculate limits
      $x0=min(($info->{'Numeric Data'})->slice("$xcol,:"));         
      $x1=max(($info->{'Numeric Data'})->slice("$xcol,:"));         
      $y0=min(($info->{'Numeric Data'})->slice("$ycol,:"));         
      $y1=max(($info->{'Numeric Data'})->slice("$ycol,:"));         
      ($info->{'Limits'})=[$x0,$x1,$y0,$y1,$xcol,$ycol,$xtext,$ytext];
      
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
      {next if /^\s*#/;  
      $x=new PDL(split " ");
 
       for ($k=7;$k<(($x->dims)[0]-1);$k+=4)
       {
 
        $y=$x->slice($k.":".($k+2));
        $ok=0;
        $nn=0;foreach (@v){++$nn;if (sum($y==$_)==3){
          if ($imax[$nn-1]<$x->at($k+3)){$imax[$nn-1]=$x->at($k+3);}
                                            $ok=1;}       
                          }
        if ($ok==0){
        $v[($#v-$[+1)]=($y);   # addy $y to list $v
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
      open($l,">./results/hkl.asc");
      print $l ("# x y T[K] |H| Ha Hb Hc [T]  vs  values of Intensity(or |m(Q)|) of ".$v[$n]."vs sqrt(int or |m(Q)|) as read from file $file\n");
      while(<$h>)
      {next if /^\s*#/;  
      $x=new PDL(split " ");
       for ($k=7;$k<(($x->dims)[0]-1);$k+=4)
       {$y=$x->slice($k.":".($k+2));
        if (sum(abs($y-$v[$n])<1e-5)==3){
        #add xyz to piddle
        $y=$x->slice("0:6");
        $z=$x->slice(($k+3).":".($k+3));
        push(@xlist,$y->append($z));
        $out=(($y->append($z))->append(sqrt($z)));
        $out=~s/\[/ /g;   
        $out=~s/\]/ /g;   
        print $l ($out."\n");
        }
       }
      }
      close $h; close $l;
      return ((cat @xlist),($v[$n]),1);
     } else {
 	print STDERR "Error program hkl: failed to read data file \"$file\"\n";
        print "\n usage:  hkl\n            or \n         hkl mcphas.hkl\n\n extracts data from reflection list mcphas.hkl";exit 1;
 	return undef;
     }
 }
 
 