#!/usr/bin/perl
use Cwd;
use File::Basename;

# implemented automatic creation of documentation output as (links in html, footnotes in tex)
# implemented is a hash of program names and labels in mcphase manual.tex
%implemented = qw(
acoscol             acoscol
add                 add
addcol              addcol
addj               addj
anisotropy          anisotropy
asincol             asincol
atancol             atancol
average             average
cfsplit             cfsplit
chi2                chi2
cif2mcphas          cif2mcphas
comment             comment
compare             compare
convolute           convolute
convolute2d         convolute2d
coscol              coscol
cpsingleion         cpsingleion
delcol              delcol
delcols             delcols
delcomments         delcomments
delnthline         delnthline
dif                dif
display            display
displaydensities   displaydensities
displaydensity     displaydensity
displayhtml         displayhtml
displaytext         displaytext
displaycontour      displaycontour
epsdebye            epsdebye
expcol              expcol
extendunitcell      extendunitcell
factcol             factcol
fermicol            fermicol
fform               fform
fillcol             fillcol
fitcol              fitcol
fitfermi            fitfermi
gauss               gauss
gauss2d             gauss2d
gausscol            gausscol
getvalue           getvalue
getvariable        getvariable
histcol             histcol
hkl                hkl
icsdread            icsdread
int                 int
linreg              linreg
lorentz             lorentz
lorentzcol          lorentzcol
makenn             makenn
mctest              mctest
mcphasit           runmcphas
-doeps             zeroepsex
mcdispit           mcdisp
mcdiff             mcdiff
newcol              newcol
newcols              newcols
newline             newline
pointc              pointc
potcol              potcol
radwavfunc          radwavfunc
range               range
reduce_unitcell     reduce_unitcell
rotate              rotate
rotateBlm           rotateBlm
rpvalue             rpvalue
script2html         script2html
setup_mcdiff_in     setup_mcdiff_in
setup_mcdisp_mf     setup_mcdisp_mf
setup_jqfit         setup_jqfit
setvalue           setvalue
setvariable        setvariable
shiftcol            shiftcol
sincol              sincol
singleion          Smagpoly
substitute          substitute
sumcol              sumcol
swapcol             swapcol
tancol              tancol
tanhcol             tanhcol
uvw2fwhm            uvw2fwhm
zshift              zshift
                 );
$URL="$ENV{'MCPHASE_DIR'}/doc/manual";
$DOC="$ENV{'MCPHASE_DIR'}/doc";
$DIR="$ENV{'MCPHASE_DIR'}";
require "$URL/labels.pl";

unless ($#ARGV >=0)
{print STDOUT << "EOF";
 program script2html used to create html documentation from McPhase scripts

 usage: script2html [-latex][options] calc1.bat [options] calc2.bat ...

 This program creates a html file from scripts containing just the text
 in the scripts. 

 input: calc1.bat, calc2.bat ...    scripts (bat files)
                                   (must be located in the current directory)
 output: stdout      ...........   html file created from the scripts
                                   (use ">" to pipe into file)

 In order to create a latex output use script2html -latex.

 options: -fromline 3 ..........  only part of the file is html coded starting at line 3
          -toline   10 .........  only part of the file is html coded (until line 10)

 Many html commands such as <h1> HEADER </h1> can be used in
 script commments to structure the text. Latex formulas are accepted by enclosing it 
 in brackets  \\( \\) for inline or \\[ \\] for equations  , e.g.  \\(e^{i\pi}+1=0\\) 
 Some abbreviations are acceptable: 
  \\( ... \\)      --> € ... €
  <ol><li>         --> /§
  </li><li>        --> §
  <li></ol>        --> §/

 Other useful html commands are:
 REM <h1> A Header </h1>
 REM <h2> A smaller header </2>
 REM <img src="figure.jpg">
 REM <figure> <figcaption> This graph shows ... 
 REM </figcaption> <img src="figure.jpg"> </figure>
 REM <pre> Text to be formatted exactly as written,  verbatim in latex</pre>
 

 example:

 script2html calc.bat notes.txt > calc.bat.html

 [creates the html file calc.bat.html from files calc.bat and notes.txt]

EOF
exit 0;}else{print STDERR "#* $0 *\n";}
if($ARGV[0]=~/-latex/){$l="-latex";shift @ARGV;}
if($l){print STDERR "#* Latex Output *\n";}else{print STDERR "#* Html Output *\n";}
 $date=localtime( time);$dir=getcwd;

if($l){$dir=~s/_/\\_/g;@DD=@ARGV;foreach(@DD){$_=~s/_/\\_/g;}
open (Fout, ">spverbatim.sty");
print Fout << "EOF";
%%
%% This is file `spverbatim.sty',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% spverbatim.dtx  (with options: `package')
%% 
%% This is a generated file.
%% 
%% Copyright (C) 2009 by Scott Pakin <scott+spverb\@pakin.org>
%% 
%% This file may be distributed and/or modified under the conditions of
%% the LaTeX Project Public License, either version 1.3c of this license
%% or (at your option) any later version.  The latest version of this
%% license is in:
%% 
%%    http://www.latex-project.org/lppl.txt
%% 
%% and version 1.3c or later is part of all distributions of LaTeX version
%% 2008/05/04 or later.
%% 
\\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\\ProvidesPackage{spverbatim}
    [2009/08/10 v1.0 Verbatim with breakable spaces]
\\gdef\\spverb{%
  \\bgroup
  \\let\\spverb\@ve=\\verb\@egroup
  \\def\\verb\@egroup{\\spverb\@ve\\egroup}%
  \\def\\\@xobeysp{\\mbox{}\\space}%
  \\verb
}
\\begingroup
  \\catcode`|=0
  \\catcode`[=1
  \\catcode`]=2
  \\catcode`\\{=12
  \\catcode`\\}=12
  \\catcode`\\\\=12
  |gdef|spv\@xverbatim#1\\end{spverbatim}[#1|end[spverbatim]]
|endgroup
\\newenvironment{spverbatim}{%
  \\def\\\@xobeysp{\\mbox{}\\space}%
  \\let\\\@xverbatim=\\spv\@xverbatim
  \\verbatim
}{%
}
\\endinput
%%
%% End of file `spverbatim.sty'.
EOF
close Fout;
open (Fout, ">placeins.sty");
print Fout << "EOF";
%  P L A C E I N S . S T Y          ver 2.2  April 18, 2005
%  Donald Arseneau                  asnd\@triumf.ca
%  Keep floats `in their place'; don't let them float into another section.
%  Instructions are below.
%
%  placeins.sty is freely released to the public domain.


\\def\\\@fb\@botlist{\\\@botlist}
\\def\\\@fb\@topbarrier{\\suppressfloats[t]}

\\catcode`\\V=14 % `V' is a comment character unless [verbose]

\\\@ifundefined{DeclareOption}{}%
{\\DeclareOption{below}{\\def\\\@fb\@botlist{}}
 \\DeclareOption{above}{\\def\\\@fb\@topbarrier{}}
 \\DeclareOption{section}{\\AtBeginDocument{%
     \\expandafter\\renewcommand\\expandafter\\section\\expandafter
       {\\expandafter\\\@fb\@secFB\\section}%
     \\newcommand\\\@fb\@secFB{\\FloatBarrier
     \\gdef\\\@fb\@afterHHook{\\\@fb\@topbarrier \\gdef\\\@fb\@afterHHook{}}}
     \\g\@addto\@macro\\\@afterheading{\\\@fb\@afterHHook}
     \\gdef\\\@fb\@afterHHook{}
  }}
 \\DeclareOption{verbose}{\\catcode`\\V=9 }% Activate things after `V'
 \\ProvidesPackage{placeins}[2005/04/18 \\space  v 2.2]
 \\ProcessOptions 
} % end of \\\@ifundefined

\\def\\FloatBarrier{\\par\\begingroup \\let\\\@elt\\relax
V\\edef\\\@tempa{\\write\\m\@ne{Package placeins Info: Float barrier, from
V  input line \\the\\inputlineno, processed on page \\thepage, lands on
V  page \\noexpand\\thepage. }}\\\@tempa
 \\edef\\\@tempa{\\\@fb\@botlist\\\@deferlist\\\@dbldeferlist}%
 \\ifx\\\@tempa\\\@empty V\\PackageInfo{placeins}{No floats held,}%
 \\else
    \\ifx\\\@fltovf\\relax % my indicator of recursion
       \\if\@firstcolumn V\\PackageWarning{placeins}{Some floats are stuck,}%
         \\clearpage 
       \\else V\\PackageInfo{placeins}{Eject a column and check again:}%
         \\null\\newpage\\FloatBarrier 
       \\fi
    \\else V\\PackageInfo{placeins}{Must dump some floats}%
       \\newpage \\let\\\@fltovf\\relax V\\PackageInfo{placeins}{Check again:}%
       \\FloatBarrier % recurse once only
 \\fi\\fi \\endgroup
 \\\@fb\@topbarrier }

\\catcode`\\V=11
\\endinput

%====================== BEGIN INSTRUCTIONS ===========================

  p l a c e i n s . s t y          ver 2.2  April 18, 2005
  Donald Arseneau                  asnd\@triumf.ca


Placeins.sty keeps floats `in their place', preventing them from floating
past a "\\FloatBarrier" command into another section.  To use it, declare
"\\usepackage{placeins}" and insert "\\FloatBarrier" at places that floats 
should not move past, perhaps at every "\\section".  

Option:  [section]

A more convenient way to stop floats at section boundaries is to change 
the definition of "\\section" to include "\\FloatBarrier", either at the
beginning, before "\\\@startsection", or in the `style' specification (see 
The LaTeX Companion, section 2.2.2; or 2.3 in the 1st ed).  If you specify 
"\\usepackage[section]{placeins}", then the "\\section" command will be 
redefined with "\\FloatBarrier" inserted at the beginning.

Options:  [above]  [below]

Something you may not like is that, by default, "\\FloatBarrier" is very 
strict, and will (try to) prevent a float from appearing above the start 
of the current section or below the start of the next section, even 
though the float is still on the same page as its intended section.  
Each restriction can be relaxed separately by using the "[above]" and 
"[below]" package options: "[above]" allows floats to appear above their 
section, if on the same page; "[below]" allows below.

NOTE!  The original version of placeins.sty acted like it was loaded
with the option "[above]" specified.

There is a problem with LaTeX's "\\suppressfloats" being out of step with 
the page breaking (see usenet msg <yfi656pbsn0.fsf\@triumf.ca> and thread)
which sometimes allows a float to go above a "\\FloatBarrier" placed near
the top of a page. Maybe placeins will fix it sometime later.

Option: [verbose]

There is a package option "[verbose]" that causes many messages to be
written in the log file.  It might be used to answer the question:
`How did *that* get *there*?!?'

%====================== END INSTRUCTIONS ========================

Test file integrity:  ASCII 32-57, 58-126:  !"#$%&'()*+,-./0123456789
:;<=>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~
EOF
close Fout;
print STDOUT << "EOF";
\\documentclass[twoside]{article}
\\hoffset-2.5cm
\\voffset-2.5cm
\\textwidth16cm
\\textheight25cm
\\oddsidemargin2.4cm
\\evensidemargin2.4cm
\\usepackage[pdftex]{graphicx}
% \\usepackage{makeidx}
\\usepackage{lscape}
\\usepackage{spverbatim}
\\usepackage{amssymb}
\\usepackage{amsbsy}
\\usepackage{afterpage}
\\usepackage{xcolor}
%\\usepackage{lipsum}
\\usepackage{placeins}
%\\usepackage[utf8]{inputenc}
%\\usepackage{amsmath}
%\\usepackage{empheq}
%\\usepackage{bbm}
%\\usepackage{dsfont}
\\usepackage{hyperref} 
%\\makeindex
\\newcommand\\mybar[1][black]{\\begingroup\\color{#1}\\kern1pt\\rule[-\\dp\\strutbox]{1pt}{\\baselineskip}\\kern1pt\\endgroup}
\\newcommand{\\highlight}[1]{\\colorbox{red!10}{\$\\displaystyle#1\$}}
\\newcommand{\\m}[1]{\\overline{#1}}
\\newcommand{\\M}[1]{\\underline{#1}}
\\newcommand{\\mbf}[1]{\\mathbf #1}
\\newcommand{\\V}[1]{ \\stackrel{=}{\\mathbf #1}}
\\newcommand{\\B}[1]{#1}
\\newcommand{\\prg}{\\sl}
\\newcommand{\\use}[1]{\\vspace{0.5cm} Usage: {\\prg{ #1}} \\vspace{0.5cm}}
\\newcommand{\\bra}[1]{\\langle #1|}
\\newcommand{\\ket}[1]{|#1\\rangle}
\\newcommand{\\threej}[2]{\\left( \\begin{array}{ccc} #1 \\\\ #2 \\end{array} \\right)}
\\newcommand{\\sixj}[2]{\\left\\{ \\begin{array}{ccc} #1 \\\\ #2 \\end{array} \\right\\}}
\\newcommand{\\hili}[1]{{#1}}
\\newcommand{\\hl}[1]{{#1}}
\\newcommand{\\Bell}{\\ensuremath{\\boldsymbol\\ell}}
\\newcommand{\\bm}[1]{\\boldsymbol #1}
\\newcommand{\\Trace}[1]{\\rm Tr \\{ #1 \\} }

\\begin{document}

 \\title{
\\includegraphics[bb=0 0 413 289,angle=0,height=90pt]{$DOC/figsrc/headerL.jpg}
\\includegraphics[bb=0 0 389 282,angle=20,height=140pt]{$DOC/figsrc/mcphase_logo.jpg}
\\includegraphics[bb=0 0 413 289,angle=0,height=90pt]{$DOC/figsrc/headerR.jpg}

Output of \\\\
 script2html -latex @DD \\\\
\\vspace{0.2cm}
{\\normalsize  ...in directory \\\\
 $dir}}
\\date{ $date , latex run on \\today }
\\author{McPhase Project\\thanks{mcphase@icloud.com}}
\\maketitle
\\pagestyle{myheadings}
\\markboth{\\hfill \\includegraphics[angle=20,height=20pt]{$DIR/mcphas_logo1.jpg} }
          {\\includegraphics[angle=20,height=20pt]{$DIR/mcphas_logo1.jpg} \\hfill}

\\tableofcontents
EOF

}
else{

print STDOUT << "EOF";
<!DOCTYPE html>
<html>
<head>
 <title>$date</title>
<style type="text/css" >
.r { font-family:'Courier',monospace; }
body { font-family:'Times',monospace;font-style=italic; }
</style>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex\@0.16.23/dist/katex.min.css" integrity="sha384-//SZkxyB7axjCAopkAL1E1rve+ZSPKapD89Lo/lLhcsXR+zOYl5z6zJZEFXil+q0" crossorigin="anonymous">

    <!-- The loading of KaTeX is deferred to speed up page rendering -->
    <script defer src="https://cdn.jsdelivr.net/npm/katex\@0.16.23/dist/katex.min.js" integrity="sha384-cpAIxua0Xbyc+XrpHQpCtJzGSZ6U2kS/FeyoKjnS+BgAYNV6uVUetVs/LC9+l3rs" crossorigin="anonymous"></script>

    <!-- To automatically render math in text elements, include the auto-render extension: -->
    <script defer src="https://cdn.jsdelivr.net/npm/katex\@0.16.23/dist/contrib/auto-render.min.js" integrity="sha384-hCXGrW6PitJEwbkoStFjeJxv+fSOOQKOPbJxSfM6G5sWZjAyWhXiTIIAmQqnlLlh" crossorigin="anonymous"
        onload="renderMathInElement(document.body);"></script>

</head><body>
<img src="$DOC/figsrc/headerL.jpg">
<img src="$DOC/figsrc/mcphase_logo.jpg">
<img src="$DOC/figsrc/headerR.jpg">
<center><h1> McPhase </h1></center>
 ...this document was created $date <br>
 ...in directory $dir<br>
 ...by the command: script2html @ARGV <br><br>

EOF
}
@ARGV=map{glob($_)}@ARGV;$i=0;$ii=0;
@BB=@ARGV;while(@BB){if($BB[0]=~/-fromline/){shift @BB; shift @BB;}
                     if($BB[0]=~/-toline/){shift @BB;shift @BB;}
                     if(!defined $l){print '<a href="#'.$BB[0].'">'.$BB[0].'</a><br>';}
                      shift @BB;
                    } 
if($l){print "\% ";$fignr=1;}
$br="<br>";
while (@ARGV)
{$linetext="";
 $fromline=1;if($ARGV[0]=~/-fromline/){shift @ARGV;$ARGV[0]=~s/x/*/g; $fromline=eval $ARGV[0];shift @ARGV;$linetext=" from line $fromline";}
 $toline=1e10;if($ARGV[0]=~/-toline/){shift @ARGV; $ARGV[0]=~s/x/*/g;$toline=eval $ARGV[0];shift @ARGV;$linetext=$linetext." up to line $toline";}
 $file=$ARGV[0];shift @ARGV; $i=$ii+1;$ii=$i;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   # get path from filename
   $dir=dirname($file);print "\n";
  if($l){$file=~s/_/\\_/g; print"\\markboth{\\hfill Source File ".$i.$linetext.":".$file." \\includegraphics[angle=20,height=20pt]{$DIR/mcphas_logo1.jpg} }
          {\\includegraphics[angle=20,height=20pt]{$DIR/mcphas_logo1.jpg} Source File ".$i.$linetext.":".$file."\\hfill}\n";
         print "\% ";
        }
  print "<!--This is a comment. Comments are not displayed in the browser END OF LINKS-->\n";
  if($l) {   print "\\subsection*{Source File ".$i.$linetext.":".$file."}\n";
       
         }
  else
  {   print '<a name="'.$file.'"><hr>Source File '.$i.$linetext.':<h1>'.$file.'</h1></a>';print "\n";}
   #print '<p class="c">';
   $lnr=0;$verbatim=0;
# *****************************************************************************
# Process lines 
# *****************************************************************************

   while(($line=<Fin>)&&$lnr<$toline)
   {++$lnr;if($lnr>=$fromline){
    if ($line=~/^\s*#/||$line=~/^\s*[rR][eE][mM]/)
# *****************************************************************************
     { # if the line starts with a comment
      if($l){    
 
           #remove comment for latex
           if($br){$line=~s/^\s*#//;
                   $line=~s/^\s*[rR][eE][mM]//;
                  }
            
           # substitute **** with \n****\n for latex
             $line=~s/(\*\*+)/\n\1\n/;
           # substitute --- with \n---\n for latex
             $line=~s/(\-\-\-\-\-\-+)/\n\1\n/;
             if($verbatim==1){$verbatim=0;$line="\\end\{spverbatim\}".$line;}
            }
# *****************************************************************************
       # take care for verbatim \ pre commands
       if ($line=~/.*\<pre\>/&&$line!=~/.*\<pre\>.*\<\/pre\>/){$br="";}
       if ($line=~/.*\<\/pre\>/&&$line!=~/.*\<\/pre\>.*\<pre\>/){$br="<br>";} 
       if($l){$line=~s/\<pre\>/\\begin{spverbatim}/g;$line=~s/\<\/pre\>/\\end{spverbatim}/g;} 
# *****************************************************************************
    if($br){
         if($l){# substitute # with \# for latex
             $line=~s/\043/\\\043/g;
                # substitute $ with \$ for latex
             $line=~s/\$/\\\$/g;

               }
        #  brackets opern close --- look if line should be continued in next line  ( not closed \( \[  € ) then
          # load next line, check if it is a comment, remove comment sign (at least for latex) and append 
           # line 
           req2("€","€","\\\(","\\\)");
           req2("\\\\\\(","\\\\\\)","\\\(","\\\)");
           req2("\\\\\\[","\\\\\\]","\\\[","\\\]");
# *****************************************************************************
         # substitute enumeration § by html commands
         $line=~s/\/§/\<ol\>\<li\>/g;
         $line=~s/\§\//\<\/li\>\<\/ol\>/g;
         $line=~s/\§/\<\/li\>\<li\>/g;
# *****************************************************************************

      if($line=~/.*\<\s*script2html.*\>/)
       { my ($bef)=($line=~m/(.*)\<\s*script2html.*\>/);print $bef;
         # look if another file should be included
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
         system("script2html $l $arguments > ".$arg[$#arg].".htm");
        if(-e  $arg[$#arg].".htm") {open(Fin1,$arg[$#arg].".htm");$line1=<Fin1>;++$ii;
                  
            until($line1=~/.*<!--This is a comment. Comments are not displayed in the browser END OF LINKS-->/){$line1=<Fin1>;}
            while($line1=<Fin1>){unless($line1=~/.*\<\/body\>\<\/html\>/||
                                        $line1=~/[^\%]*\\bibliographystyle\{/   ||
                                        $line1=~/[^\%]*\\bibliography/  
                                       ){
                                 if($l){$line1=~s/subsection\*\{Source File\s*(.*)\}/subsection\*{Source File $ii\.\1\}
             \\markboth\{\\hfill Source File $ii\.\1 \\includegraphics\[angle=20,height=20pt\]\{$DIR\/mcphas_logo1.jpg\} \}
            \{\\includegraphics\[angle=20,height=20pt\]\{$DIR\/mcphas_logo1.jpg\} Source File $ii\.\1 \\hfill\}\n/; # 
                                         $line1=~s/\\end\{document\}//;
                                       }
                                 else
                                 {$line1=~s/\<hr\>Source File\s*/\<hr\>Source File $ii\./;}
                                 print $line1;}}
        close Fin1;unlink($arg[$#arg].".htm");
                                 if($l){print "\\subsection*{continuing Source File $i $linetext $file}\n";}
                                  else {print '<hr>Continuing Source File '.$i.' '.$linetext.':<h1>'.$file.'</h1>';print "\n";} 
                                    }
                       else { print stderr "Error script2html: unable to open ".$arg[$#arg].".htm\n";}
       }
# *****************************************************************************
        else
       { # take care about <img src=""> commands and insert path
        if($l){if($line=~/\<figure\>/){$figure=1;$line=~s/\<figure\>/see fig.\\ref\{fig$fignr\}\n\\begin\{figure\}[ht]\\begin\{center\}/;}
               if($line=~/\<figcaption\>/){$line=~s!\<figcaption\>!\\caption\{\\label\{fig$fignr\}\n!;++$fignr;}
               if($line=~/.*\<\/figcaption\>/){$line=~s!\<\/figcaption\>!\}\n!;}
               if($line=~/\<img(.*)src\s*=/){($filename)=($line=~m|\<img.*src\s*=\s*"([^\s^>^<^=]+)"|);
                if($filename=~/\.gif$/){# notneeded ..$filen=$filename; $filen=~s/\\_/_/g; ($heightref, $widthref) = gifdim($filen);
                                       # convert gif files to jpg so they can be processed by pdflatex
                                       system("giftopnm $filename | pnmtojpeg > $filename.jpg");$filename.=".jpg";
                                      } $filename=~s/\\_/_/g;
                if($figure==1){ 
              $line=~s!(\s*#?\s*)\<img(.*)src\s*="([^"]*)"[^\>]*\>!\\includegraphics[angle=0,width=0.6\\columnwidth]\{$filename\}!;
                }else
                { $line=~s!(\s*#?\s*)\<img(.*)src\s*="([^"]*)"[^\>]*\>!see fig.\\ref\{fig$fignr\}
           \\begin\{figure\}[ht]\\begin\{center\}
           \\includegraphics[angle=0,width=0.6\\columnwidth]\{$dir/\3\}
           \\caption{\\label\{fig$fignr\}
            $filename}
           \\end\{center\}
           \\end\{figure\}
           \\afterpage\{\\FloatBarrier\}!; ++$fignr;  
                 }
             }
          if($line=~/.*\<\/figure\>/){$figure=0;$line=~s!\<\/figure\>!\\end\{center\}\\end\{figure\}\\afterpage\{\\FloatBarrier\}!;}
               
        }
        else
        {$line=~s!(\s*#?\s*)\<img(.*)src\s*="(.*)"!<p style="width:50%;word-wrap: break-word; "> \1 &lt img\2src="$dir/\3"&gt </p> \<img\2src="$dir/\3"!;
         }
# *****************************************************************************
# replace html commands
       if($l){# replace html commands <...> by nothing
        $line=~s/\<h1\>/\\section\{/g; $line=~s/\<\/h1\>/\}/g;
        $line=~s/\<h2\>/\\subsection\{/g; $line=~s/\<\/h2\>/\}/g;
        $line=~s/\<h3\>/\\subsubsection\{/g; $line=~s/\<\/h3\>/\}/g;
        $line=~s/\<h4\>/\\paragraph\{/g; $line=~s/\<\/h4\>/\}/g;
        $line=~s/\<h5\>/\\subparagraph\{/g; $line=~s/\<\/h5\>/\}/g;
        $line=~s/\<sub\>/\\(_\{/g;$line=~s/\<\/sub\>/\}\\)/g;
        $line=~s/\<ol\>/\\begin\{itemize\}/g;$line=~s/\<\/ol\>/\\end\{itemize\}/g;
        $line=~s/\<li\>/\\item /g;$line=~s/\<\/li\>//g;
        
        $line=~s/\<(\/?)(a|b|q|caption|center|cite|code|col|
                         |dd|del|dfn|div|dl|dt|em|fieldset|figure|figcaption|form|frame|
                         |h1|h2|h3|h4|h5|h6|head|hr|html|img|iframe|input|ins|label|legend|li|
                         |map|meta|noframes|noscript|object|ol|optgroup|option|
                         |p|pre|small|span|sub|sup|table|tbody|textarea|tfoot|th|title|td|tr|tt|u|ul|var)([^\>]*?)\>//g; 
        $line=~s/\<(\/?)([i])(\s*?)\>//g;# html tag <i>
        }
        else
        { # replace html commands <...> by &aaa& ... &bbb& 
        $line=~s/\<(\/?)(a|b|q|caption|center|cite|code|col|
                         |dd|del|dfn|div|dl|dt|em|fieldset|figure|figcaption|form|frame|
                         |h1|h2|h3|h4|h5|h6|head|hr|html|img|iframe|input|ins|label|legend|li|
                         |map|meta|noframes|noscript|object|ol|optgroup|option|
                         |p|pre|small|span|sub|sup|table|tbody|textarea|tfoot|th|title|td|tr|tt|u|ul|var)([^\>]*?)\>/&aaa&\1\2\3&bbb&/g; 
        $line=~s/\<(\/?)([i])(\s*?)\>/&aaa&\1\2\3&bbb&/g;# html tag <i>
        }
       if($l)
       { 
# *****************************************************************************
 #unless we are in an equation treat _ & ^ | < > # $ symbols
# ...  first remove all normal brackets except \( \) \[ \]
$line=~s/(?<!\\)\[/myleftrectangularbracket/g;
$line=~s/(?<!\\)\]/myrightrectangularbracket/g;
$line=~s/(?<!\\)\(/myleftangularbracket/g;
$line=~s/(?<!\\)\)/myrightangularbracket/g;
             # substitute underscore with \_ for latex - use lookahead to exclude being between () []
             #  brackets \( \) or \[ \]  math mode of latex
             $line=~s/_(?![^\[\]\(\)]*\\[\]\)])/\\_/g; # s/_/\\_/g;
             # substitute & with \& for latex
             $line=~s/\&(?![^\[\]\(\)]*\\[\]\)])/\\\&/g; # ~s/\&/\\\&/g;
             # substitute ^ with \^ for latex
             $line=~s/\^(?![^\[\]\(\)]*\\[\]\)])/\\\^/g; # ~s/\^/\\\^/g;
             # substitute | with $|$ for latex
             $line=~s/\|(?![^\[\]\(\)]*\\[\]\)])/\$\|\$/g; # ~s/\|/\$\|\$/g;
             # substitute <> with $<$ $>$ for latex
             $line=~s/>(?![^\[\]\(\)]*\\[\]\)])/\$\>\$ /g;
             $line=~s/<(?![^\[\]\(\)]*\\[\]\)])/\$\<\$ /g;
$line=~s/myleftrectangularbracket/\[/g; # substitute back all brackets
$line=~s/myrightrectangularbracket/\]/g;
$line=~s/myleftangularbracket/\(/g;
$line=~s/myrightangularbracket/\)/g;

           #do substituion to get in latex an equation\( \) \[ \] should become $ and begin equation ...
             $line=~s/\\\(/\$/g; # inline math
             $line=~s/\\\)/\$/g;
             $line=~s/\\\[/\n\\begin\{equation\}\n/g; # equation
             $line=~s/\\\]/\n\\end\{equation\}\n/g; 
             

        }else{
       # substitute all remaining < and > signs by the html code &gt and &lt
        $line=~s/>/&gt /g;$line=~s/</&lt /g; 
       # replace back &aaa& ... &bbb& to < ... > so that html commands are interpreted properly
        $line=~s/&aaa&/\</g;$line=~s/&bbb&/\>/g;
       $line=~s/\n/$br\n/g;  # print comments in style "c" (default)
       }
       print  $line;
       }
# *****************************************************************************
    } # fi $br
    else
    { print  $line;  }
    }else{ 
# line did not start with a comment - thus it is a command and should be printed as it is
if($l) {if($line=~/\S/&&$verbatim==0){$verbatim=1;$line="\\begin\{spverbatim\} ".$line;}
      my @to_delete;
      foreach(keys %implemented)
        {my $com=$_;  
                my  $comr=$com;$comr=~s/_/\\_/g;
         if($line=~/.*\s$_\s/){# $_  matches a command ? --> insert a footnote with exlanation of the command
                          # and delete command from hash %implemented so footnotes do not double on next use
                   # scan doc/*.tex for %script2html_begin{singleion} some text to be processed
                   #                    %script2html_end{singleion} some text to be processed
                   my $ftexfile="results/".$_.".tex";
                   open(FOUT, '>', $ftexfile);
                    opendir my $dir, $DOC; my @files = readdir $dir;
                    foreach(@files){if($_=~/\.tex$/){$store=0;
                                      open(FH, '<',$DOC."/".$_) or die $!;while(<FH>){
                                      if($_=~/\%script2html_begin\{$com\}/){$store=1;$_=~s/\%script2html_begin\{$com\}/$comr:/;}
                                      if($_=~/\%script2html_end\{$com\}/){$store=0;$_=~s/\%script2html_end\{$com\}//;print FOUT $_;}
                                      if($store==1){ print FOUT $_;}
                                     } close FH;
                                   }
                                   }
                   close FOUT;
                  $line=$line."\\end\{spverbatim\}\n \\dots for details on $comr see footnote \\footnote\{\\input\{".$ftexfile."\}\}\n\\begin\{spverbatim\}";                  
                  push @to_delete, $com;
                   }
         }  
       
       # delete keys which occured already
       foreach (@to_delete) {delete($implemented{$_}); }
       }else{
   

    # substitute all  < and > signs by the html code &gt and &lt
        $line=~s/>/&gt /g;$line=~s/</&lt /g;
   $line='<span class="r"> '.$line.' </span>'.$br; #print commands in style "r"
   foreach(keys %implemented)
   {
    if($line=~/.*\s$_\s/){# $_  matches a command ? --> insert a link to the formula in 
                    $label=$implemented{$_}; # this is the label of an equation etc
                    $link=$external_labels{$label}."#$label";  # this is the link to manual/node...html#label
                    # insert the link here
                    $line=~s/$_/\<A HREF="$link"\>$_\<\/A\>/;           
                   }
    }   

  }
   print  $line;
   }
  
   }} # print "</p>\n";
close Fin; 
if($verbatim==1){$verbatim=0;print "\\end\{spverbatim\}";}
       
} 
close Fout;
if($l){

print "\\bibliographystyle\{".$DOC."\/physrev\}   \% here you should update any list by\n";
print "\\bibliography\{".$DOC."\/li120914\}   \% bibtex - ing the database\n";

print "\\end\{document\}\n";}
else{
print "<hr>\n";
print "</body></html>\n";
 }

sub gifdim ($) {
    my $filename = $_[0];

    open(GIF, $filename) || return (undef, undef);
    my $buf = '';
    my $n = read GIF, $buf, 10;
    close GIF;

    return (undef, undef) if $n < 10;
    my ($head, $width, $height) = unpack("A6vv", $buf);
    return (undef, undef) unless $head =~ /^GIF8[79]a/;
    return \($width, $height);
}


sub req2   # check if all $S are closed by $Z symbols and substitute $S by $SS and $Z by $ZZ in $line
{my ($S,$Z,$SS,$ZZ)=@_;
          while(cu($S,$Z)){if (!($more=<Fin>)||!($more=~/^\s*#/||$more=~/^\s*[rR][eE][mM]/)) 
                                       { die "Error reading  line $lnr unclosed $S \n$line\n";}
                               ++$lnr;
                                  $more=~s/^\s*#//;$more=~s/^\s*[rR][eE][mM]//; # remove comment
                                  $line=$line.$more; # attach
                                  }

          if($line=~/$S/){#print STDERR $line."1\n";
 $line=~s/$S(.*?)$Z/$SS\1$ZZ/sg;
#print STDERR $line."2\n";
                          }
                
}
    
sub cu   # check if all $S are closed by $Z symbols in $line
{my ($S,$Z)=@_;my $check=$line;
 while($check=~/$S(.*?)$Z/s) {#print STDERR $check."3\n";
                              $check=~s/$S(.*?)$Z/\1/s;} 
if($check=~/$S/){#print STDERR "cu  $S $Z true\n";
return true;}else{#print STDERR "cu $S $Z false\n";
return undef;}
}