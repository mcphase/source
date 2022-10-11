#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

# batch to remove every nth line (n>2) in output.txt

# usage: nth_rm n < input.txt > output.txt 


unless ($#ARGV >0) 

{print " program to remove every nth line (n>2) in file
\n";

 print " usage: delnthline n *.*\n";

 exit 0;}else{print STDERR "#* $0 *\n";}

# first~step

#               Match  every step'th line starting with line first.

#                 For example, ``sed -n 1~2p''  will  print  all  the

#               odd-numbered  lines  in  the  input stream, and the

#                 address 2~5 will match every fifth  line,  starting

#               with the second.  (This is a GNU extension.)

								                                   



$ARGV[0]=~s/x/*/g;$n=eval $ARGV[0];shift @ARGV;


 foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}  
   open (Fout, ">range.out");
   $i=0;
   while($line=<Fin>){
   ++$i;unless($i==$n){ print Fout $line;}
   else {$i=0;}
                     }

      close Fin;

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

 }
