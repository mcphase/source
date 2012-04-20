#!/usr/bin/perl

use Cwd;

unless ($#ARGV >=0)
{print STDOUT << "EOF";
 program script2html used to create html documentation from McPhase scripts

 usage: script2html calc1.bat calc2.bat ...

 This program creates a html file from scripts containing just the text
 in the scripts.

 input: calc1.bat, calc2.bat ...    scripts (bat files)
                                   (must be located in the current directory)
 output: stdout      ...........   html file created from the scripts
                                   (use ">" to pipe into file)
 example:

 script2html calc.bat notes.txt > report.html

 [creates the html file report.html from files calc.bat and notes.txt]

EOF
exit 0;}
 $date=localtime( time);$dir=getcwd;
print STDOUT << "EOF";
<html>
<head>
 <title>$date</title>
<style type="text/css">
.r { font-family:'Times';font-style=italic; }
.c { font-family:'Courier',monospace; }
</style>

</head><body>
 ... this document was created $date <br>
...  in directory $dir<br>
...  by the command: script2html @ARGV <br><br>
EOF
@ARGV=map{glob($_)}@ARGV;$i=0;
foreach (@ARGV) 
{$file=$_; ++$i;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<hr>Source File $i:<h1>".$file."</h1>\n";
   print '<p class="c">';
   while($line=<Fin>)
   {if ($line=~/^\s*#/||$line=~/^\s*rem/){
    $line='<span class="r">'.$line.'</span><br>'
    }else{
   $line=~s/\n/<br>\n/g;
   }
   print stdout $line;
   } print "</p>";
close Fin; 
} 
close Fout;
print "</body></html>\n";
 



