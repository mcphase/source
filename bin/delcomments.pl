#!/usr/bin/perl

# batch to remove every comment line 





unless ($#ARGV >=0) 

{print " program delcomments  used to remove all comment lines (starting with #) in a file, removed comments are printed to screen (stdout)\n";

 print " usage: delcomments [-s \"token\"] *.*  \n *.* .. filenname\n";
 print " -s ... option to remove comments beginning with \"token\" instead.\n";
 print "        e.g. delcomments -s \"#:\" to remove lines commented by the range command.\n";
 print " -k ... option to keep comments beginning with \"token\" .\n";
 print "        e.g. delcomments -k \"#!\" to keep comments starting with '#!' \n";
 print " -c ... option to remove comments and replace by empty line instead of deleteing the line\n";
 print " -fromline 10 ... option to remove comments only after line number 10\n";
 print " -toline 100  ... option to remove comments only until line number 100\n";
 exit 0;}else{print STDERR "#* $0 *";}



$delline=1;$commentstring = "#"; $fromline=1;$toline=1e100;$keep="1";
while($ARGV[0]=~/^-/)
{if ($ARGV[0]=~/^-c$/) { $delline=0;shift @ARGV;}
 if ($ARGV[0]=~/^-s$/) { shift @ARGV; $commentstring = $ARGV[0]; shift @ARGV; } 
 if ($ARGV[0]=~/^-k$/) { shift @ARGV; $keep = $ARGV[0]; shift @ARGV; } 
 if ($ARGV[0]=~/^-fromline$/) { $fromline = $ARGV[1]; shift @ARGV; shift @ARGV; } 
 if ($ARGV[0]=~/^-toline$/) { $toline = $ARGV[1]; shift @ARGV; shift @ARGV; } 
}
#print "s".$commentstring."s\n";
@ARGV=map{glob($_)}@ARGV;

  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print STDERR "<".$file;

   open (Fout, ">range.out");
   $i=0;
   while($line=<Fin>)
     {++$i;
       if ($fromline<=$i&&$i<=$toline&&$line=~/^\s*\Q$commentstring\E/
           &&!($line=~/^\s*\Q$keep\E/)) {print $line;
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
print STDERR ">";
   }

print STDERR "\n";