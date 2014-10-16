#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

unless ($#ARGV >2) {usage();}

sub usage() {

  print STDERR << "EOF";

    $0: program to set a value in a data file

    usage: $0 row column text file[s]

     row         : row number
     column      : column number
     text        : text to be placed in file at this position
     files       : one or more filenames

    example: $0 4 5 3.142 data.dat
             replaces the number in row 4 and column 5 by 3.142 in file data.dat
EOF

  exit;

}
 

$ARGV[0]=~s/x/*/g;$row=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$column=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$const=eval $ARGV[0];shift @ARGV;



  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
   open (Fout, ">range.out");$ri=0;
   while($line=<Fin>)
     { if ($line=~/^\s*#/) {print Fout $line;}
       else{$ri+=1;
            if($ri==$row)
               {$line=~s/D/E/g;@numbers=split(" ",$line);
                $numbers[$column-1]=$const;
           	  $i=0;
		  foreach (@numbers)
		  {++$i;print Fout $numbers[$i-1]." ";}
                  print Fout "\n";
               }
             else {print Fout $line;}
           }
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

 