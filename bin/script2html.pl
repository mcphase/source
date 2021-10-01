#!/usr/bin/perl
use Cwd;
use File::Basename;

unless ($#ARGV >=0)
{print STDOUT << "EOF";
 program script2html used to create html documentation from McPhase scripts

 usage: script2html [options] calc1.bat [options] calc2.bat ...

 This program creates a html file from scripts containing just the text
 in the scripts.

 input: calc1.bat, calc2.bat ...    scripts (bat files)
                                   (must be located in the current directory)
 output: stdout      ...........   html file created from the scripts
                                   (use ">" to pipe into file)

 options: -fromline 3 ..........  only part of the file is html coded starting at line 3
          -toline   10 .........  only part of the file is html coded (until line 10)
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
<style type="text/css" >
.r { font-family:'Courier',monospace; }
body { font-family:'Times',monospace;font-style=italic; }
</style>


</head><body>
 ...this document was created $date <br>
 ...in directory $dir<br>
 ...by the command: script2html @ARGV <br><br>

EOF
@ARGV=map{glob($_)}@ARGV;$i=0;
@BB=@ARGV;while(@BB){if($BB[0]=~/-fromline/){shift @BB; shift @BB;}
                     if($BB[0]=~/-toline/){shift @BB;shift @BB;}
                     print '<a href="#'.$BB[0].'">'.$BB[0].'</a><br>';shift @BB;
                    } 
print "<!--This is a comment. Comments are not displayed in the browser END OF LINKS-->";
$br="<br>";
while (@ARGV)
{$linetext="";
 $fromline=1;if($ARGV[0]=~/-fromline/){shift @ARGV;$ARGV[0]=~s/x/*/g; $fromline=eval $ARGV[0];shift @ARGV;$linetext=" from line $fromline";}
 $toline=1e10;if($ARGV[0]=~/-toline/){shift @ARGV; $ARGV[0]=~s/x/*/g;$toline=eval $ARGV[0];shift @ARGV;$linetext=$linetext." up to line $toline";}
 $file=$ARGV[0];shift @ARGV; ++$i;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   # get path from filename
   $dir=dirname($file);print "\n";
   print '<a name="'.$file.'"><hr>Source File '.$i.$linetext.':<h1>'.$file.'</h1></a>';print "\n";
   #print '<p class="c">';
   $lnr=0;
   while(($line=<Fin>)&&$lnr<$toline)
   {++$lnr;if($lnr>=$fromline){
    if ($line=~/^\s*#/||$line=~/^\s*[rR][eE][mM]/)
     { # if the line starts with a comment

       # take care for verbatim \ pre commands
       if ($line=~/.*\<pre\>/&&$line!=~/.*\<pre\>.*\<\/pre\>/){$br="";}
       if ($line=~/.*\<\/pre\>/&&$line!=~/.*\<\/pre\>.*\<pre\>/){$br="<br>";} 
 
      if($line=~/.*\<\s*script2html.*\>/)
       { # look if another file should be included
         # if yes run script2htlm on this file
         ($arguments)=($line=~m/.*\<\s*script2html(.*)\>/);
          @arg=split(" ",$arguments);
          foreach(@arg){$a=$_;
                        next if($a=~/-fromline/);
                        next if($a=~/-toline/);
                        next if($aa=~/-fromline/);
                        next if($aa=~/-toline/);
                        # put directory name in front of filenames
                        $_=$dir."/".$_;
                       } continue {
                        $aa=$a;
                       } $arguments=join(' ',@arg);
        # print "script2html $arguments > ".$arg[$#arg].".htm\n";
         system("script2html $arguments > ".$arg[$#arg].".htm");
        if(-e  $arg[$#arg].".htm") {open(Fin1,$arg[$#arg].".htm");$line1=<Fin1>;
            until($line1=~/.*<!--This is a comment. Comments are not displayed in the browser END OF LINKS-->/){$line1=<Fin1>;}
            while($line1=<Fin1>){unless($line1=~/.*\<\/body\>\<\/html\>/){
                                 $line1=~s/\<hr\>Source File\s*/\<hr\>Source File $i\./;
                                 print $line1;}}
        close Fin1;unlink($arg[$#arg].".htm");
                                    }
                       else { print stderr "Error script2html: unable to open ".$arg[$#arg].".htm\n";}
       }
        else
       {
        # take care about <img src=""> commands and insert path
        $line=~s!(\s*#?\s*)\<img(.*)src\s*="(.*)"!<p style="width:50%;word-wrap: break-word; "> \1 &lt img\2src="$dir/\3"&gt </p> \<img\2src="$dir/\3"!;
        # replace html commands <...> by &aaa& ... &bbb& 
        $line=~s/\<(\/?)(a|b|q|caption|center|cite|code|col|
                         |dd|del|dfn|div|dl|dt|em|fieldset|form|frame|
                         |h1|h2|h3|h4|h5|h6|head|hr|html|img|iframe|input|ins|label|legend|li|
                         |map|meta|noframes|noscript|object|ol|optgroup|option|
                         |p|pre|small|span|sub|sup|table|tbody|textarea|tfoot|th|title|td|tr|tt|u|ul|var)([^\>]*?)\>/&aaa&\1\2\3&bbb&/g; 
        $line=~s/\<(\/?)([i])(\s*?)\>/&aaa&\1\2\3&bbb&/g;# html tag <i>

       # substitute all remaining < and > signs by the html code &gt and &lt
        $line=~s/>/&gt /g;$line=~s/</&lt /g; 
       # replace back &aaa& ... &bbb& to < ... > so that html commands are interpreted properly
        $line=~s/&aaa&/\</g;$line=~s/&bbb&/\>/g;
       $line=~s/\n/$br\n/g;  # print comments in style "c" (default)
       print  $line;
       }
    }else{ # line did not start with a comment - thus it is a command and should be printed as it is
    # substitute all  < and > signs by the html code &gt and &lt
        $line=~s/>/&gt /g;$line=~s/</&lt /g;
   $line='<span class="r">'.$line.'</span>'.$br; #print commands in style "r"
   print  $line;
   }
  
   }} # print "</p>\n";
close Fin; 
} 
close Fout;
print "<hr>\n";
print "</body></html>\n";
 



