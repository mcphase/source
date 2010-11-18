#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

# batch to remove every comment line 





unless ($#ARGV >=0) 

{print " program delcomments  used to remove all comment lines (starting with #) in a file, removed comments are printed to screen (stdout)\n";

 print " usage: delcomments *.*  \n *.* .. filenname\n";

 exit 0;}



#$command="#addc";foreach $d(@ARGV){$command.= " ".$d;}; $command.="\n";



  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;


   open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print $line;}

       else{print Fout $line;}

       

      }

      close Fin;

#      print Fout $command;

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

print ">\n";

   }

