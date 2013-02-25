#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

# batch to remove every comment line 





unless ($#ARGV >=0) 

{print " program delcomments  used to remove all comment lines (starting with #) in a file, removed comments are printed to screen (stdout)\n";

 print " usage: delcomments [-s \"token\"] *.*  \n *.* .. filenname\n";
 print " -s ... option to remove comments beginning with \"token\" instead.\n";
 print "        e.g. delcomments -s \"#:\" to remove lines commented by the range command.\n";
 print " -c ... option to remove comments and replace by empty line instead of deleteing the line\n";
 exit 0;}



$delline=1;
if ($ARGV[0]=~/^-c$/) { $delline=0;shift @ARGV;}
if ($ARGV[0]=~/^-s$/) { $commentstring = $ARGV[1]; shift @ARGV; shift @ARGV; } else { $commentstring = "#"; }
if ($ARGV[0]=~/^-c$/) { $delline=0;shift @ARGV;}

  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;


   open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*$commentstring/) {print $line;
                                         if($delline==0){print Fout "\n";}
                                         }

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

