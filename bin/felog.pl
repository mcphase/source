#! /usr/bin/perl 



# Perl / PDL program to plot the contents of mcphas.log 



# declare perl modules to be used

use FileHandle;

use PDL;

use PDL::Graphics::TriD;

use PDL::Graphics::TriD::Image;

use PDL::IO::Pic;



# turn verbose off

$PDL::Graphics::TriD::verbose=0;



# some global variables

my ($header);

my ($c0);

my($c1);

my ($h);

my ($hsave);



my ($file) = @ARGV;  # read filename from commandline

$file = "./mcphas.log" unless $file; # default filename



print <<ENDE;

 Rotate the image by pressing mouse button one and

 dragging in the graphics window.

 Zoom in/out by pressing MouseButton3 and drag up/down.   

 Press 'q' in the graphics window for hardcopy in file

 felog.jpg and  next logged HT-point screen.

 Press Alt-F4 to exit

ENDE



# main loop

for (;;) {

    plot_array_2D($file,1);

}



sub plot_array_2D {

    my ($info,$loadset) = @_;

    my ($data,$ramp,$Limits,$header) = get_detector_data_2D($info,$loadset);

    my ($x0,$x1,$dx,$y0,$y1,$dy,$z0,$z1,$dz,$c0,$c1,$dc) = @{$Limits};



# set display window size

$PDL::Graphics::TriD::GL::xsize=700;

$PDL::Graphics::TriD::GL::ysize=700;



# create new graph

$g = new PDL::Graphics::TriD::Graph();

$g->default_axes();

# "pts" ist just a name of the dataseries

$g->add_dataseries(new PDL::Graphics::TriD::Points($data,$ramp),"pts");

$g->bind_default("pts");

$g->scalethings();



# set window

my $win = PDL::Graphics::TriD::get_current_window();



#old

#$lb = PDL::Graphics::OpenGL::glpRasterFont("8x16",0,255);

#new

$lb = $win->glpRasterFont("8x16",0,255);

     print "dd\n";

$win->clear_objects();

$win->add_object($g);

$win->add_object(new TOBJ());



# twiddling loop util user presses q

$win->twiddle();



# save graphic to file   

   my $pic=grabpic3d();

   wpic $pic,"./felog.jpg";

# open next window

 $win =new PDL::Graphics::TriD::GL::Window();

$lb = $win->glpRasterFont("8x16",0,255);


}





sub get_detector_data_2D {

    my ($info,$loadset) = @_;

     if ($loadset)

    {($info->{'Numeric Data'},$info->{'header'}) = read_data_file_2D($info);

     # calculate limits     

     $xcol=0;

     $ycol=1;

     $zcol=2;

     $ccol=3;

     

     $t=$info->{'Numeric Data'};@dim=$t->dims;

     $x0=min(abs($info->{'Numeric Data'})->slice("$xcol,:"));         

     $x1=max(abs($info->{'Numeric Data'})->slice("$xcol,:"));         

     $y0=min(abs($info->{'Numeric Data'})->slice("$ycol,:"));         

     $y1=max(abs($info->{'Numeric Data'})->slice("$ycol,:"));         

     $z0=min(abs($info->{'Numeric Data'})->slice("$zcol,:"));         

     $z1=max(abs($info->{'Numeric Data'})->slice("$zcol,:"));         

     $c0=min(($info->{'Numeric Data'})->slice("$ccol,:"));         

     $c1=max(($info->{'Numeric Data'})->slice("$ccol,:"));         

     

     # calculate  step intervals

     $dx=0;for ($i=0;$dx==0&&$i<$dim[1]-1;++$i)

              {$dx=abs(-$t->at($xcol,$i)+$t->at($xcol,$i+1));}

     if ($dx==0) {$dx=1;$x1=$x0+1;}

     $dy=0;for ($i=0;$dy==0&&$i<$dim[1]-1;++$i)

              {$dy=abs(-$t->at($ycol,$i)+$t->at($ycol,$i+1));}

     if ($dy==0) {$dy=1;$y1=$y0+1;}

     $dz=0;for ($i=0;$dz==0&&$i<$dim[1]-1;++$i)

              {$dz=abs(-$t->at($zcol,$i)+$t->at($zcol,$i+1));}

     if ($dz==0) {$dz=1;$z1=$z0+1;}

     $dc=0;for ($i=0;$dc==0&&$i<$dim[1]-1;++$i)

              {$dc=abs(-$t->at($ccol,$i)+$t->at($ccol,$i+1));}

     if ($dc==0) {$dc=1;$c1=$c0+1;}

     

     # set location piddle $b and color piddle $t

    $b=abs($t->slice("0:2"))->xchg(0,1);

    $cz=sqrt((($t->slice("3,:")-$c0)/($c1-$c0))->clump(2));

     $rc = sin(-$cz*5.3+1.8)/2 + 0.5;

     $gc = sin(-$cz*5.3+3.9)/2 + 0.5;

     $bc = sin(-$cz*5.3-0.3)/2 + 0.5;

     $t=new PDL (cat $rc,$gc,$bc);

     $t=$t->reorder(1,0)->xchg(0,1);

     

     # add additional points to make them better visible

     $bb=$b;

     $tt=$t;

   for($ddx=-1;$ddx<1.0001 ;$ddx+=0.5 ){

    $ym=sqrt(-$ddx*$ddx+1.0001);

    for($ddy=-$ym ;$ddy<$ym ;$ddy+=0.5){

     $zm=sqrt(-$ddy*$ddy-$ddx*$ddx+1.0002);

     for($ddz=-$zm ;$ddz<$zm ;$ddz+=0.5 ){

     $rr= new PDL $b;$rr*=0;

     $rx=$rr->slice(":,0");

     $rx+=$ddx*$dx/50;

     $ry=$rr->slice(":,1");

     $ry+=$ddy*$dy/50;

     $rz=$rr->slice(":,2");

     $rz+=$ddz*$dz/50;

    

     $bb=$bb->append($b+$rr);

     $tt=$tt->append($t);

      }}}

      # add 2 points in order to get scales at display right

     $bb=$bb->append(new PDL [[$x1],[$y1],[$z1]]);

     $bb=$bb->append(new PDL [[$x0],[$y0],[$z0]]);

     $tt=$tt->append(new PDL [[0],[0],[0]]);

     $tt=$tt->append(new PDL [[0],[0],[0]]);

     $info->{'Numeric Data'}=$bb->xchg(0,1);

     ($info->{'Colors'})=$tt->xchg(0,1);



     ($info->{'Limits'})=[$x0,$x1,$dx,$y0,$y1,$dy,$z0,$z1,$dz,$c0,$c1,$dc];

     } 

    return ($info->{'Numeric Data'},$info->{'Colors'},$info->{'Limits'},$info->{'header'});

}





# Read 2D numeric data, skipping comment lines.

sub read_data_file_2D {

    my ($file) = @_;

    my @list = ();

    $header=$hsave;

    

    unless ($h)  # if file is not opened yet ...

    {$h = new FileHandle;

     unless (open($h, $file))

     {print STDOUT "Warning: failed to read data file \"$file\"\n";

     return undef;}

     #input comment line

     $header=<$h>; 

    }

    

     print $header;

     while (<$h>) {   # read data into buffer and return

            if (/^\s*#/) {$hsave=$_;return ((cat @list),($hsave));}

	    unless (/^.*2000/) {push(@list, new PDL (split " "));}

	}

     close $h;    

     print STDOUT "End of file \"$file\"\n";    

     $hsave=$_;return ((cat @list),($hsave));;

}





# the following package is to add text on the graph

package TOBJ;

BEGIN{@ARGV=map{glob($_)}@ARGV}{@TOBJ::ISA = qw/PDL::Graphics::TriD::Object/;}

use PDL::Graphics::OpenGLQ;

use OpenGL;



sub new {

	bless {},$_[0];

}



sub togl {

	glDisable(&GL_LIGHTING);

	line_3x_3c(

		$::nx,$::nc

	);

	glColor3f(1,0,1);

	glRasterPos3f(0,0.1,1.1);

	OpenGL::glpPrintString($::lb,$header);

	glRasterPos3f(0,0.1,0.9);

	OpenGL::glpPrintString($::lb,"fe=".$c0);

	glRasterPos3f(0,0.1,0.9);

	OpenGL::glpPrintString($::lb," "x(23)."fe=".$c1."meV");

        for ($i=0;$i<1;$i+=0.1)

	{

	glColor3f(sin(-$i*5.3+1.8)/2 + 0.5,sin(-$i*5.3+3.9)/2 + 0.5,sin(-$i*5.3-0.3)/2 + 0.5);

	glRasterPos3f(0,0.1,0.9);

	OpenGL::glpPrintString($::lb," "x(12+10*$i)."*");



	}



	glEnable(&GL_LIGHTING);

}

