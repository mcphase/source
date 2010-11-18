#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}





unless ($#ARGV >0) 

{print " program plotbook to arrange many (single side) ps-images to one booklet\n";

 print " needs: latex, pstops, perl\n"; 

 print " usage: plotbook *.ps\n";

 exit 0;}



open (Fout,"> plotbookfig.tex");

$i=0;

  foreach (@ARGV)

  {++$i;

   $file=$_;

   $fileeps=$file;

   $fileeps=~s/\Q.ps\E/.eps/;

   if (!-e $fileeps){system("ps2epsi ".$file." ".$fileeps);}

   print Fout "\n".'\\begin{figure}[btp]%h=here, t=top, b=bottom, p=separate figure page'."\n";

   print Fout '\\begin{center}\\leavevmode'."\n";

   print Fout '\\includegraphics[angle=0, width=1.0\\textwidth]{'.$fileeps.'}'."\n";

#   print Fout '\\caption{'.$dirtex.'/setup'.$i.'.ps}'."\n".'\\label{'.$dir.'/setup'.$i.'}'."\n";

   print Fout '\\end{center}\\end{figure}'."\n\n";

   if($i==6){print Fout '\\clearpage'."\n";$i=0;}

  }



close Fout;



open (Fout,"> plotbook.tex");

print Fout 

'%\documentclass[twoside]{article}

\documentclass[a4paper,10pt,twoside]{article}

\hoffset -2.5cm    %-2.5

\voffset -1cm    %-4

\textwidth16cm

\textheight26cm

\oddsidemargin2.4cm

\evensidemargin2.4cm

\usepackage{graphicx,rotating}

\usepackage{epsfig}

\newcommand{\mbf}[1]{\mbox{\bf $#1$}}

\newcommand{\prg}{\textit}

\newcommand{\use}[1]{\vspace{0.5cm} Usage: {\prg{ #1}} \vspace{0.5cm}}



\begin{document}





\input{plotbookfig.tex}



\end{document}';

close Fout;

 

system ("latex plotbook");

system ("dvips plotbook");

system ("pstops '4:0L@.7(21cm,0)+1L@.7(21cm,14.85cm),2L@.7(21cm,0)+3L@.7(21cm,14.85cm)' plotbook.ps plotbook2.ps");

system ("pstops '4:0L@.7(21cm,0)+1L@.7(21cm,14.85cm),2L@.7(21cm,0)+3L@.7(21cm,14.85cm)' plotbook2.ps plotbook4.ps");

system ("pstops '4:0L@.7(21cm,0)+1L@.7(21cm,14.85cm),2L@.7(21cm,0)+3L@.7(21cm,14.85cm)' plotbook4.ps plotbook8.ps");

system ("rm plotbook.tex plotbook.aux plotbook.dvi plotbook.log plotbookfig.tex");

#pstops '4:0L@.7(21cm,0)+1L@.7(21cm,14.85cm),2L@.7(21cm,0)+3L@.7(21cm,14.85cm)' many1.ps many2.ps

#pstops '4:0L@.7(21cm,0)+1L@.7(21cm,14.85cm),2L@.7(21cm,0)+3L@.7(21cm,14.85cm)' many2.ps many3.ps



#pstops  "2:0L@.65(21cm,0)+1L@.65(21cm,14.85cm)" many1 > dp.ps

#pstops  "2:-1" plot.ps > dev.ps   

#pstops  "2:0" many1.ps > dodd.ps







  

#\end{verbatim} 