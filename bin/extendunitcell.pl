#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use FileHandle;
use PDL;
use Getopt::Long;
#use PDL::Slatec;



unless ($#ARGV>=2) 
{print " program to extend crystallographic unit cell n times in r1 (or r2,r3) direction\n\n";
print " usage:  extendunitcell 3 1 4 [options]\n\n";
print " meaning take mcphas.j and generate an extended description of the unit cell 3xr1,1xr2,4xr3\n";
print " put result into results/extend.j\n
options for creating quantum dots, quantum chains and quantum planes: \n
     -r1    remove interactions beyond the extended unit cell in r1 direction\n
     -r2    remove interactions beyond the extended unit cell in r2 direction\n
     -r3    remove interactions beyond the extended unit cell in r3 direction\n
";
 exit 0;}else{print STDERR "#* $0 *\n";}


$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;my ($n1) = eval $ARGV[0];shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;my ($n2) = eval $ARGV[0];shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;my ($n3) = eval $ARGV[0];shift @ARGV;

# -------------------------------------------------------------------------------------- #
# Options and declarations sets mri to 1 if option is present in command line
# -------------------------------------------------------------------------------------- #
GetOptions("r2"=>\$mr2,"r1"=>\$mr1,
	   "r3"=>\$mr3);

 print "reading mcphas.j ....\n";
 my ($latt,$p) = getlattice("./mcphas.j");
 my ($a,$b,$c,$nofa,$nofcomponents) = @{$latt};
 print "a=".$a." b=".$b." c=".$c."\n";
 print "original primitive lattice[abc]:".$p."\n";

 #initialize output file extend.j
 my ($l)=printlattice("./mcphas.j",$n1,$n2,$n3,">./results/extend.j");
 printneighbourlist("./mcphas.j",$l,$n1,$n2,$n3,$p,$nofa);   
 close $l;
   $nofatoms=$nofa*$n1*$n2*$n3;
 # extend lattice
 $p->slice(",0")*=$n1;
 $p->slice(",1")*=$n2;
 $p->slice(",2")*=$n3;

 print "new primitive lattice[abc]:".$p."\n";

 $exit;
  
#-----------------------------------------------------------------------
# Get lattic data, reading it from file 
sub getlattice {
    my ($file) = @_;
    my $h = new FileHandle;
    my $n = 0;
    $nofcomponents=0;
  # input data int piddle
  if(open($h,$file))
  {      while(<$h>)
     {#next if /^\s*#/;
      # detect a= b= c= ...
         
      if ($a==0){($a)=extract("a",$_);}
      if ($b==0){($b)=extract("b",$_);}
      if ($c==0){($c)=extract("c",$_);}

      if ($r1x==0){($r1x)=extract("r1a",$_);}
      if ($r1y==0){($r1y)=extract("r1b",$_);}
      if ($r1z==0){($r1z)=extract("r1c",$_);}
      if ($r2x==0){($r2x)=extract("r2a",$_);}
      if ($r2y==0){($r2y)=extract("r2b",$_);}
      if ($r2z==0){($r2z)=extract("r2c",$_);}
      if ($r3x==0){($r3x)=extract("r3a",$_);}
      if ($r3y==0){($r3y)=extract("r3b",$_);}
      if ($r3z==0){($r3z)=extract("r3c",$_);}

      if ($r1x==0){($r1x)=extract("r1x",$_);}
      if ($r1y==0){($r1y)=extract("r1y",$_);}
      if ($r1z==0){($r1z)=extract("r1z",$_);}
      if ($r2x==0){($r2x)=extract("r2x",$_);}
      if ($r2y==0){($r2y)=extract("r2y",$_);}
      if ($r2z==0){($r2z)=extract("r2z",$_);}
      if ($r3x==0){($r3x)=extract("r3x",$_);}
      if ($r3y==0){($r3y)=extract("r3y",$_);}
      if ($r3z==0){($r3z)=extract("r3z",$_);}

      if ($nofatoms==0){($nofatoms)=extract("nofatoms",$_);}
      if ($nofcomponents==0){($nofcomponents)=extract("nofcomponents",$_);}

      if (/^(#!|[^#])*nofneighbours\s*=\s*/){++$n;

                                   ($nofneighbours[$n])=extract("nofneighbours",$_);
                                   ($x[$n])=extract("da",$_);
                                   ($y[$n])=extract("db",$_);
                                   ($z[$n])=extract("dc",$_);
                                   if (/^.*\Qx=\E/){($x[$n])=extract("x",$_);}
				   if (/^.*\Qy=\E/){($y[$n])=extract("y",$_);}
				   if (/^.*\Qz=\E/){($z[$n])=extract("z",$_);}
				     ($sipffilename)=extractstring("sipffilename",$_);
                                     ($charge[$n])=extractfromfile("CHARGE",$sipffilename);
                                     if($charge[$n]==""){$charge[$n]=$sipffilename;}                                 
                                             #               print "$sipffilename  charge=".$charge[$n]."\n";
                                                           

				  }

     }
     close $h; 
     if ($n!=$nofatoms) {print STDOUT "Failed to read data file \"$file\": wrong number of atoms\n";
                         return undef;}

     if($nofcomponents==0) {$nofcomponents=3;}
     return ([$a,$b,$c,$nofatoms,$nofcomponents],pdl [[$r1x,$r1y,$r1z],[$r2x,$r2y,$r2z],[$r3x,$r3y,$r3z]]);
    } else {
	print STDOUT "Warning: failed to read data file \"$file\"\n";
	return undef;
    }
}


sub printlattice {
  my ($filein,$n1,$n2,$n3,$fileout)=@_;   
#   print $rn->index($n)."\n";
    my $h = new FileHandle;
    my $l = new FileHandle;
    open($l,$fileout);
    open($h,$filein);
     while(<$h>)
     {$text=$_;
      
      ($r1x)=extract("r1x",$text);if($r1x!=""){$r1x*=$n1;$text=~s/\Qr1x=\E\s*[\-\+\d.]+/r1a=$r1x /;}
      ($r1y)=extract("r1y",$text);if($r1y!=""){$r1y*=$n1;$text=~s/\Qr1y=\E\s*[\-\+\d.]+/r1b=$r1y /;}
      ($r1z)=extract("r1z",$text);if($r1z!=""){$r1z*=$n1;$text=~s/\Qr1z=\E\s*[\-\+\d.]+/r1c=$r1z /;}

      ($r2x)=extract("r2x",$text);if($r2x!=""){$r2x*=$n2;$text=~s/\Qr2x=\E\s*[\-\+\d.]+/r2a=$r2x /;}
      ($r2y)=extract("r2y",$text);if($r2y!=""){$r2y*=$n2;$text=~s/\Qr2y=\E\s*[\-\+\d.]+/r2b=$r2y /;}
      ($r2z)=extract("r2z",$text);if($r2z!=""){$r2z*=$n2;$text=~s/\Qr2z=\E\s*[\-\+\d.]+/r2c=$r2z /;}

      ($r3x)=extract("r3x",$text);if($r3x!=""){$r3x*=$n3;$text=~s/\Qr3x=\E\s*[\-\+\d.]+/r3a=$r3x /;}
      ($r3y)=extract("r3y",$text);if($r3y!=""){$r3y*=$n3;$text=~s/\Qr3y=\E\s*[\-\+\d.]+/r3b=$r3y /;}
      ($r3z)=extract("r3z",$text);if($r3z!=""){$r3z*=$n3;$text=~s/\Qr3z=\E\s*[\-\+\d.]+/r3c=$r3z /;}

      ($r1a)=extract("r1a",$text);if($r1a!=""){$r1a*=$n1;$text=~s/\Qr1a=\E\s*[\-\+\d.]+/r1a=$r1a /;}
      ($r1b)=extract("r1b",$text);if($r1b!=""){$r1b*=$n1;$text=~s/\Qr1b=\E\s*[\-\+\d.]+/r1b=$r1b /;}
      ($r1c)=extract("r1c",$text);if($r1c!=""){$r1c*=$n1;$text=~s/\Qr1c=\E\s*[\-\+\d.]+/r1c=$r1c /;}

      ($r2a)=extract("r2a",$text);if($r2a!=""){$r2a*=$n2;$text=~s/\Qr2a=\E\s*[\-\+\d.]+/r2a=$r2a /;}
      ($r2b)=extract("r2b",$text);if($r2b!=""){$r2b*=$n2;$text=~s/\Qr2b=\E\s*[\-\+\d.]+/r2b=$r2b /;}
      ($r2c)=extract("r2c",$text);if($r2c!=""){$r2c*=$n2;$text=~s/\Qr2c=\E\s*[\-\+\d.]+/r2c=$r2c /;}

      ($r3a)=extract("r3a",$text);if($r3a!=""){$r3a*=$n3;$text=~s/\Qr3a=\E\s*[\-\+\d.]+/r3a=$r3a /;}
      ($r3b)=extract("r3b",$text);if($r3b!=""){$r3b*=$n3;$text=~s/\Qr3b=\E\s*[\-\+\d.]+/r3b=$r3b /;}
      ($r3c)=extract("r3c",$text);if($r3c!=""){$r3c*=$n3;$text=~s/\Qr3c=\E\s*[\-\+\d.]+/r3c=$r3c /;}

      ($nofatoms)=extract("nofatoms",$text);if($nofatoms!=""){$nofatoms*=$n1*$n2*$n3;$text=~s/\Qnofatoms=\E\s*[\-\+\d.]+/nofatoms=$nofatoms /;}
      last if /^(#!|[^#])*nofneighbours\s*=\s*/;
      print $l ($text);      
     }
 close $h;
 return ($l);
}



sub printneighbourlist {
  my ($filein,$l,$n1,$n2,$n3,$p,$nofatoms)=@_;   
  $inv=inv(transpose($p));
  open($h,$filein);
     while(<$h>)
     {last if /^(#!|[^#])*nofatoms\s*=\s*/;}
     while(<$h>)
     {last if /^#.*\Q**********\E/;}
  for($n=1;$n<=$nofatoms;++$n)
  { 
   # get a neighbour
   $nn=0;
     while(<$h>)
     {last if /^#.*\Q**********\E/;
      ++$nn;
      $tt[$nn]=$_;
      ($dd)=extract("nofneighbours",$tt[$nn]);
      if($dd!=""){$nofnorig=$dd;}
      }
    
   for ($i1=0;$i1<=$n1-1;++$i1){
   for ($i2=0;$i2<=$n2-1;++$i2){
   for ($i3=0;$i3<=$n3-1;++$i3){
     $da=$x[$n]+$i1*$p->at(0,0)+$i2*$p->at(0,1)+$i3*$p->at(0,2);
     $db=$y[$n]+$i1*$p->at(1,0)+$i2*$p->at(1,1)+$i3*$p->at(1,2);
     $dc=$z[$n]+$i1*$p->at(2,0)+$i2*$p->at(2,1)+$i3*$p->at(2,2);
     # check options -ri if interactions go outside extended unit cell
     $nofn=$nofnorig;
     for($i=1;$i<=$nn;++$i)
     {$text[$i]=$tt[$i];
      unless($text[$i]=~/^#/)
      {@numm=split(" ",$text[$i]);
          $dabc=pdl [$numm[0],$numm[1],$numm[2]];
          $dr123=($dabc x $inv)->slice(":,(0)");
       if($mr1){if($dr123->at(0)<-$i1){$text[$i]="";--$nofn;}
               if($dr123->at(0)>$n1-1-$i1){$text[$i]="";--$nofn;}
               }
       if($mr2&&$text[$i]){if($dr123->at(1)<-$i2){$text[$i]="";--$nofn;}
               if($dr123->at(1)>$n2-1-$i2){$text[$i]="";--$nofn;}
               }
     #  print $dr123->at(2)." ".(-$i3)." ".($n3-1-$i3)."\n";
       if($mr3&&$text[$i]){if($dr123->at(2)<-$i3){$text[$i]="";--$nofn;}
               if($dr123->at(2)>$n3-1-$i3){$text[$i]="";--$nofn;}
               }
     }
     }
      # print out atom and neighbor list
      for($i=1;$i<=$nn;++$i)
       {
        $text[$i]=~s!\Qnofneighbours=\E\s*[\-\+\d.]+!nofneighbours=$nofn!;
        $text[$i]=~s!\Qda=\E\s*[\-\+\d.]+!da=$da!;
        $text[$i]=~s!\Qdb=\E\s*[\-\+\d.]+!db=$db!;
        $text[$i]=~s!\Qdc=\E\s*[\-\+\d.]+!dc=$dc!;
        print $l ($text[$i]);
       }
     print $l ("#*************************************************************************\n");
   }}}
  } 
}

# **********************************************************************************************
# extracts number from string
# 
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
# 
# ... it stores 0.3 in the variable $standarddeviation
#
sub extract { 
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************

# **********************************************************************************************
# extracts string from string
# 
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
# 
# ... it stores 0.3 in the variable $standarddeviation
#
sub extractstring { 
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*\b([^\n\s]+)[\s\n]/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************


# **********************************************************************************************
# extracts number from file
# 
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:         
#
# ($standarddeviation)=extractfromfile("sta","data.dat");
# 
# ... it stores 0.24 in the variable $standarddeviation
#
sub extractfromfile {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";$value="";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);}}
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
