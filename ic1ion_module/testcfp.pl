#!/usr/bin/perl

unless (-e "ING11.CFP")   { die "Reference cfp tables in file ING11.CFP from RCG by Rober Cowan does not exist"; } 
unless (-e "testcfp.exe") { system('g++ cfpout.cpp cfp.cpp states.cpp -o testcfp.exe')==0 or die "Fail to compile cfpout program."; }

# Gets the block of cfps for the current number of f-electrons
if ($#ARGV<0) { $n=2 } else { $n=$ARGV[0]; } 
if ($#ARGV>0) { $l=2 } else { $l=3; }

open(CFP,"ING11.CFP"); $startflag=0; $stopflag=0;
while(<CFP>)
{
   chomp;
   if($l==3)
   {
      if($_=~/^\s*F/)
      {
         @Line = split;
         if($n==$Line[1]) { $startflag=1; @headerline=@Line; }
         if($n==$Line[1]-1) { $stopflag=1; }
      }
   }
   else
   {
      if($_=~/^\s*D/)
      {
         @Line = split;
         if($n==$Line[1]) { $startflag=1; @headerline=@Line; }
         if($n==$Line[1]-1) { $stopflag=1; }
      }
   }
   if($startflag==1 && $stopflag==0)
   {
      push @cfprefblock, $_;
   }
   if(stopflag==1) { last; }
}
close(CFP);

$nchild = $headerline[2]; $nparent = $headerline[3]; @cfpref = (); @mycfpref = ();
for $i (1..$nchild) { for $j (1..$nparent) { $cfpref[$i][$j] = 0; $mycfpref[$i][$j] = 0; } }
$nparlines=int($nparent/20)+1; if($n<10) { $endcol=12; } else { $endcol=11; }

$count=1; @childstates = ();
for(@cfprefblock)
{
   @Line = split;
   if($count>1 && $count<=($nparlines+1))
   {
      push @parentstates, @Line; 
      if($n==8 && $l==3) { for $j (1..$#parentstates) 
      { 
         if($parentstates[$j] eq "2F0") { $parentstates[$j]="2F10"; } 
         if($parentstates[$j] eq "2G0") { $parentstates[$j]="2G10"; last; } 
      } }
   }
   elsif($count>$nparlines)
   {
      if(abs($Line[1])>999)
      {
         $n2=abs($Line[1])-abs(int($Line[1]/1000)*1000); $Line[1]=int($Line[1]/1000); splice(@Line,2,0,$n2);
      }
      if($Line[4]>$nparent || ($Line[4]<1 && $Line[4]!=0)) { $n2=int($Line[3]*1e10)/1e10; $n1=int(abs((abs($n2)-abs($Line[3]))*1e13)+0.5); splice(@Line,3,1,$n2,$n1); }
      if($Line[6]>$nparent || ($Line[6]<1 && $Line[6]!=0)) { $n2=int($Line[5]*1e10)/1e10; $n1=int(abs((abs($n2)-abs($Line[5]))*1e13)+0.5); splice(@Line,5,1,$n2,$n1); }
      if($Line[8]>$nparent || ($Line[8]<1 && $Line[8]!=0)) { $n2=int($Line[7]*1e10)/1e10; $n1=int(abs((abs($n2)-abs($Line[7]))*1e13)+0.5); splice(@Line,7,1,$n2,$n1); }
      if(abs($Line[1])==1) 
      { 
         if($Line[0] eq "2F0") { $Line[0]="2F10"; } 
         if($Line[0] eq "2G0") { $Line[0]="2G10"; } 
         $childstates[$Line[$endcol]-1] = $Line[0]; 
      }
      if($Line[2]!=0) { $cfpref[$Line[$endcol]-1][$Line[2]-1] = $Line[3]; }
      if($Line[4]!=0) { $cfpref[$Line[$endcol]-1][$Line[4]-1] = $Line[5]; }
      if($Line[6]!=0) { $cfpref[$Line[$endcol]-1][$Line[6]-1] = $Line[7]; }
      if($Line[8]!=0) { $cfpref[$Line[$endcol]-1][$Line[8]-1] = $Line[9]; }
   }
   $count++;
}

# Determines the cfps calculated by ic1ion
@mycfp = `./testcfp.exe $n $l`; $count=1; $ci=0;

for(@mycfp)
{
   chomp; @Line = split;
   if($count==1)
   {
      @myparentstates = @Line;
      for $iL (0..$#parentstates) { for $jL (0..$#parentstates) 
      {
         if($Line[$jL] eq $parentstates[$iL]) { push @myparind, $jL; last; }
      } }
   }
   else
   {
      if(abs($Line[1])==1) 
      { 
         for $iL (0..$#childstates) { if($Line[0] eq $childstates[$iL]) { $mychildstates[$iL] = $Line[0]; $ci=$iL; last; } }
      }
      if($Line[2]!=0) { $mycfpref[$ci][$Line[2]-1] = $Line[3]; }
      if($Line[4]!=0) { $mycfpref[$ci][$Line[4]-1] = $Line[5]; }
      if($Line[6]!=0) { $mycfpref[$ci][$Line[6]-1] = $Line[7]; }
      if($Line[8]!=0) { $mycfpref[$ci][$Line[8]-1] = $Line[9]; }
   }
   $count++
}

print "\t\t\tChecking cfp's for configuration "; if($l==2) { print "d"; } else { print "f"; } print "$n\n";
print "--------------------------------------------------------------------------------------------------------\n";
$sum=0; print "Child-Parent\tRCG_Index\tic1ion_Index\tRCG_Value\tic1ion_Value\tAbs(Difference)\n";
print "--------------------------------------------------------------------------------------------------------\n";
for $i (0..$nchild-1) { for $j (0..$nparent-1)
{
   $sum += abs($mycfpref[$i][$myparind[$j]]-$cfpref[$i][$j]);
   if(abs($mycfpref[$i][$myparind[$j]]-$cfpref[$i][$j])>1e-6)
   {
      print "$childstates[$i]-$parentstates[$j]     \t(",$i+1,",",$j+1,")   \t(",$i+1,",",$myparind[$j]+1,")    \t";
      if($cfpref[$i][$j]==0) { print "0\t\t"; } else { print "$cfpref[$i][$j]\t"; }
      if($mycfpref[$i][$myparind[$j]]==0) { print "0\t\t"; } else { print "$mycfpref[$i][$myparind[$j]]\t"; }
      print abs($mycfpref[$i][$myparind[$j]]-$cfpref[$i][$j]),"\n";
   }
} }
print "--------------------------------------------------------------------------------------------------------\n";
print "Sum of Abs(Difference) = $sum\n";
print "--------------------------------------------------------------------------------------------------------\n";
