#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use FileHandle;
use PDL;
use PDL::Slatec;



unless ($#ARGV>=2) 
{print " program to extend crystallographic unit cell n times in r1 (or r2,r3) direction\n\n";
print " usage: extendunitcell 3 1 4\n\n";
print " meaning take mcphas.j, mcphas.tst and mcdiff.in and generate an extended description of the unit cell 3xr1,1xr2,4xr3\n";
print " put result into results/extend.j, results/extend.tst and results/extend.head (only header without magnetic atoms is generated from mcdiff.in)\n";
 exit 0;}


$ARGV[0]=~s/x/*/g;my ($n1) = eval $ARGV[0];shift @ARGV; 
$ARGV[0]=~s/x/*/g;my ($n2) = eval $ARGV[0];shift @ARGV; 
$ARGV[0]=~s/x/*/g;my ($n3) = eval $ARGV[0]; 

 print "reading mcphas.j ....\n";
 my ($latt,$p) = getlattice("./mcphas.j");
 my ($a,$b,$c,$nofa,$nofcomponents) = @{$latt};
 print "a=".$a." b=".$b." c=".$c."\n";
 print "original primitive lattice[abc]:".$p."\n";

 #initialize output file extendj.j
 my ($l)=printlattice("./mcphas.j",$n1,$n2,$n3,">./results/extend.j");
 printneighbourlist("./mcphas.j",$l,$n1,$n2,$n3,$p,$nofa);   
 close $l;
   
 # extend lattice
 $p->slice(",0")*=$n1;
 $p->slice(",1")*=$n2;
 $p->slice(",2")*=$n3;

 print "new primitive lattice[abc]:".$p."\n";
if (open ($h, "mcphas.tst")){
 print "reading mcphas.tst ...\n";
 $nofcomponents=0;   
 open ($l,">./results/extend.tst"); 
 while(<$h>)
 {if (/^\s*#/)
  {$text=$_;
  if($nofcomponents==0){($nofcomponents)=extract("nofcomponents",$text);}
  $text=~s/\Qnofatoms=\E\s*[\d.]+/nofatoms=$nofatoms /;
  print $l ($text);
  }
  else
  { # read configuration
   my @list =();
   push(@list,new PDL([split " "]));
   @dims=$list[0]->dims; $nr1=$dims[0];
   print @list;
   print @dims;
   print "  $nr1\n";
    for($i=2;$i<=$nofcomponents;++$i)
     {$_=<$h>;push(@list,new PDL(split " "));}
   $nr2=1;$nr3=1;
      while(<$h>)
   {last if !/^\s*\d/;
    push(@list,new PDL(split " "));
    for($i=2;$i<=$nofcomponents;++$i)
     {$_=<$h>;push(@list,new PDL(split " "));}
    $nr2++;
   } 
   if (!/^\s*#/)
   {while(<$h>)
    {last if /^\s*#/;
     $nr3++;
     for ($j=1;$j<=$nr2;++$j)
     {    for($i=1;$i<=$nofcomponents;++$i)
          {push(@list,new PDL(split " "));$_=<$h>;}
     }
     last if /^\s*#/;
    }
   }

   $conf=new PDL(cat @list);
   # extend configuration
   # write configuration
   for($ii3=1;$ii3<=$n3;++$ii3){for($i3=1;$i3<=$nr3;++$i3){
   for($ii2=1;$ii2<=$n2;++$ii2){for($i2=1;$i2<=$nr2;++$i2){
    for($i=1;$i<=$nofcomponents;++$i){
     for($ii1=1;$ii1<=$n1;++$ii1){for($i1=1;$i1<=$nr1;++$i1){
     $j=$i1;$offset=$nofcomponents*(($i3-1)*$nr2+$i2-1);
     print $l $conf->at($j-1,$offset+$i-1)." ";
     }} print $l "\n";
   }}}if($i3<$nr3||$ii3<$n3) {print $l "\n";}}}
   print $l $_;
  }
 } 
 close $h,$l;
}
else {print "\n unable to open mcphas.tst\n";}   

unless (open ($h,"mcdiff.in")) {print "\nunable to read mcdiff.in\n";exit(0);}
 print "reading mcdiff.in ...\n";
 open ($l,">./results/extend.head"); 
 while(<$h>)
 {$text=$_;
  if (/^\s*#/)
  {
   if (/^.*\Qnat=\E/){($nat)=(/^.*\Qnat=\E\s*([^\s]*)\s/);
                      $nofatoms=$nat*$n1*$n2*$n3;$text=~s/\Qnat=\E\s*[\d.]+/nat=$nofatoms /;}
   print $l ($text);
  }
  else
  {@numbers=split(" ");
   $x=$numbers[2];
   $y=$numbers[3];
   $z=$numbers[4];
   $dr1=$numbers[5];
   $dr2=$numbers[6];
   $dr3=$numbers[7];

   for ($i1=0;$i1<=$n1-1;++$i1){
   for ($i2=0;$i2<=$n2-1;++$i2){
   for ($i3=0;$i3<=$n3-1;++$i3){
     $da=$x+$i1*$p->at(0,0)+$i2*$p->at(0,1)+$i3*$p->at(0,2);
     $db=$y+$i1*$p->at(1,0)+$i2*$p->at(1,1)+$i3*$p->at(1,2);
     $dc=$z+$i1*$p->at(2,0)+$i2*$p->at(2,1)+$i3*$p->at(2,2);
     $dr1n=($dr1+$i1)/$n1;
     $dr2n=($dr2+$i2)/$n2;
     $dr3n=($dr3+$i3)/$n3;
   print $l $numbers[0]." ".$numbers[1]." ".$da." ".$db." ".$dc." ".$dr1n." ".$dr2n." ".$dr3n." ".$numbers[8]."\n";
   }}} 
  }
  last if /.*SECTION\s*3/;
 }
 
 close $h,$l;
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
                                   ($gJ[$n])=extract("gJ",$_);
				   if (/^.*\Qx=\E/){($x[$n])=extract("x",$_);}
				   if (/^.*\Qy=\E/){($y[$n])=extract("y",$_);}
				   if (/^.*\Qz=\E/){($z[$n])=extract("z",$_);}
				     ($cffilename)=extractstring("cffilename",$_);
                                     ($charge[$n])=extractfromfile("CHARGE",$cffilename);
                                     if($charge[$n]==""){$charge[$n]=$cffilename;}                                 
                                             #               print "$cffilename  charge=".$charge[$n]."\n";
                                                           

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
      
      ($r1x)=extract("r1x",$text);if($r1x!=""){$r1x*=$n1;$text=~s/\Qr1x=\E\s*[\d.]+/r1a=$r1x /;}
      ($r1y)=extract("r1y",$text);if($r1y!=""){$r1y*=$n1;$text=~s/\Qr1y=\E\s*[\d.]+/r1b=$r1y /;}
      ($r1z)=extract("r1z",$text);if($r1z!=""){$r1z*=$n1;$text=~s/\Qr1z=\E\s*[\d.]+/r1c=$r1z /;}

      ($r2x)=extract("r2x",$text);if($r2x!=""){$r2x*=$n2;$text=~s/\Qr2x=\E\s*[\d.]+/r2a=$r2x /;}
      ($r2y)=extract("r2y",$text);if($r2y!=""){$r2y*=$n2;$text=~s/\Qr2y=\E\s*[\d.]+/r2b=$r2y /;}
      ($r2z)=extract("r2z",$text);if($r2z!=""){$r2z*=$n2;$text=~s/\Qr2z=\E\s*[\d.]+/r2c=$r2z /;}

      ($r3x)=extract("r3x",$text);if($r3x!=""){$r3x*=$n3;$text=~s/\Qr3x=\E\s*[\d.]+/r3a=$r3x /;}
      ($r3y)=extract("r3y",$text);if($r3y!=""){$r3y*=$n3;$text=~s/\Qr3y=\E\s*[\d.]+/r3b=$r3y /;}
      ($r3z)=extract("r3z",$text);if($r3z!=""){$r3z*=$n3;$text=~s/\Qr3z=\E\s*[\d.]+/r3c=$r3z /;}

      ($r1a)=extract("r1a",$text);if($r1a!=""){$r1a*=$n1;$text=~s/\Qr1a=\E\s*[\d.]+/r1a=$r1a /;}
      ($r1b)=extract("r1b",$text);if($r1b!=""){$r1b*=$n1;$text=~s/\Qr1b=\E\s*[\d.]+/r1b=$r1b /;}
      ($r1c)=extract("r1c",$text);if($r1c!=""){$r1c*=$n1;$text=~s/\Qr1c=\E\s*[\d.]+/r1c=$r1c /;}

      ($r2a)=extract("r2a",$text);if($r2a!=""){$r2a*=$n2;$text=~s/\Qr2a=\E\s*[\d.]+/r2a=$r2a /;}
      ($r2b)=extract("r2b",$text);if($r2b!=""){$r2b*=$n2;$text=~s/\Qr2b=\E\s*[\d.]+/r2b=$r2b /;}
      ($r2c)=extract("r2c",$text);if($r2c!=""){$r2c*=$n2;$text=~s/\Qr2c=\E\s*[\d.]+/r2c=$r2c /;}

      ($r3a)=extract("r3a",$text);if($r3a!=""){$r3a*=$n3;$text=~s/\Qr3a=\E\s*[\d.]+/r3a=$r3a /;}
      ($r3b)=extract("r3b",$text);if($r3b!=""){$r3b*=$n3;$text=~s/\Qr3b=\E\s*[\d.]+/r3b=$r3b /;}
      ($r3c)=extract("r3c",$text);if($r3c!=""){$r3c*=$n3;$text=~s/\Qr3c=\E\s*[\d.]+/r3c=$r3c /;}

      ($nofatoms)=extract("nofatoms",$text);if($nofatoms!=""){$nofatoms*=$n1*$n2*$n3;$text=~s/\Qnofatoms=\E\s*[\d.]+/nofatoms=$nofatoms /;}
      last if /^(#!|[^#])*nofneighbours\s*=\s*/;
      print $l ($text);      
     }
 close $h;
 return ($l);
}



sub printneighbourlist {
  my ($filein,$l,$n1,$n2,$n3,$p,$nofatoms)=@_;   
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
      $text[$nn]=$_;
      }
    
   for ($i1=0;$i1<=$n1-1;++$i1){
   for ($i2=0;$i2<=$n2-1;++$i2){
   for ($i3=0;$i3<=$n3-1;++$i3){
     $da=$x[$n]+$i1*$p->at(0,0)+$i2*$p->at(0,1)+$i3*$p->at(0,2);
     $db=$y[$n]+$i1*$p->at(1,0)+$i2*$p->at(1,1)+$i3*$p->at(1,2);
     $dc=$z[$n]+$i1*$p->at(2,0)+$i2*$p->at(2,1)+$i3*$p->at(2,2);
     # print out neighbor
     for($i=1;$i<=$nn;++$i)
       {
        $text[$i]=~s!\Qda=\E\s*[\d.]+!da=$da!;
        $text[$i]=~s!\Qdb=\E\s*[\d.]+!db=$db!;
        $text[$i]=~s!\Qdc=\E\s*[\d.]+!dc=$dc!;
        print $l ($text[$i]) 
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
