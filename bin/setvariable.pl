#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

unless ($#ARGV >1) {usage();}

sub usage() {

  print STDERR << "EOF";

    $0: program to set a variable in a data file, e.g. replace T=15 with T=20
        in some file

    usage: $0 variablename value file[s]

     variablename: name of the variable (e.g. "T")
     value       : value to be set (e.g. "20")
     files       : one or more filenames

    example: $0 T 20 mcdisp.mf
             replaces T=something with T=20 in file mcdisp.mf
EOF

  exit;

}
 

$varnam=$ARGV[0];shift @ARGV;
$value=$ARGV[0];shift @ARGV;



  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
   open (Fout, ">range.out");$ri=0;
   while($line=<Fin>)
     { if ($line=~/^(#!|[^#])*?\b$varnam/) {
                                            #here write modified parameter set to line
                                             $line=~s!$varnam\s*=\s*[^\s\;\n\t\*]+!$varnam=$value!;
                                           }
       print Fout $line;
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
      print ">\n";
 }   

 