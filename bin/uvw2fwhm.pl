#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



use Math::Trig;
unless ($#ARGV >0) 
{print STDOUT << "EOF";

 program uvw2fwhm  used to calculate the full width half maximum from
                   2theta in degree according to fullprofs u,v,w parameters ...

 usage: uvw2fwhm u v w col *.*

 u v w  ... uvw parameters as in fullprof
 col    ... column
 *.*    ... filenname

 The column col in file *.* must contain 2theta scattering angle
 values and is overwritten with the fwhm as calculated by

 fwhm = sqr(u tan^2(theta) + v tan(theta) + w)

EOF
 exit 0;}
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$u=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$v=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$w=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$column=eval $ARGV[0];shift @ARGV;


  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

           	  $i=0;++$j;

		  foreach (@numbers)

		  {++$i;

		  if ($i==$column) {
  $zt = $numbers[$i-1] * 3.14159 / 180; #zwei theta in rad
  if ($u * tan($zt / 2) * tan($zt / 2) + $v * tan($zt / 2) + $w > 0)
  {  $fwhm = sqrt($u * tan($zt / 2) * tan($zt / 2) + $v * tan($zt / 2) + $w);
  }else
  { print "warning: at x=".$numbers[$i-1]."degree the fwhm given by parameters uvw gets negative ... continuing\n";
  }                                   $numbers[$i-1]=$fwhm;}

		  print Fout $numbers[$i-1]." ";}     

            print Fout "\n";

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



#\end{verbatim} 