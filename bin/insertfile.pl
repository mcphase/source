#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

unless ($#ARGV >1)

{print " program insertfile  used to insert file2 [and file 3 ..] to file1\n";
 print " usage: insertfile linenumber file1 file2 [file3 ...] \n";
 exit 0;}else{print STDERR "#* $0 *\n";}

$ARGV[0]=~s/x/*/g;$linenr=eval $ARGV[0];shift @ARGV;
$file1=$ARGV[0];shift @ARGV;
print "inserting to $file1 at line $linenr:\n";

   unless (open (Fin, $file1)){die "\n error:unable to open $file1\n";}
  $i=1;open (Fout, ">range.out");
   while(($line=<Fin>)&&$i<$linenr)
     {print Fout $line;++$i}

  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin1, $file)){die "\n error:unable to open $file\n";}
   print "<".$file;
   while($line=<Fin1>)
     {print Fout $line;}
      close Fin1;
   print ">\n";
 }

   while($line=<Fin>)
     {print Fout $line;++$i}

   unless (rename "range.out",$file1)

          {unless(open (Fout, ">$file1"))

      {die "\n error:unable to write to $file1\n";}

      open (Fin, "range.out");

      while($line=<Fin>){ print Fout $line;}

      close Fin;

      close Fout;

      system "del range.out";

     }