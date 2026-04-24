/* ---------------------
 * modified BubbleChartDemo1.java
 * to implement display program
 * M . Rotter 2010
 * ---------------------
 * (C) Copyright 2003-2008, by Object Refinery Limited.
 */


//package demo; 
import java.awt.geom.Point2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;

import java.awt.Color;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import javax.imageio.ImageIO;

//import javax.swing.JButton;
//import javax.swing.SwingConstants;
import java.io.*;
import java.util.EventListener;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.LegendItemSource;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.axis.Axis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.DefaultIntervalXYDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.chart.ui.UIUtils;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.HorizontalAlignment;

import expr.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;

//import com.sun.image.codec.jpeg.JPEGCodec;
//import com.sun.image.codec.jpeg.JPEGImageEncoder;

/**
 * A bubble chart demo.
 */
public class display extends ApplicationFrame implements KeyListener,WindowListener  {

static final int MAX_NOF_FILES = 20;
static myStringfunc SF=new myStringfunc();
static int xy[]={0,0,0,0};
 static String gVlines = "";
 static String gHlines = "";
//static Frame frame;
//static Frame popup;
//static ToolTipManager ToolTipManager;

  static class MyChartMouseListener implements ChartMouseListener {
 ChartPanel panel;

        /**
         * Creates a new mouse listener.
         *
         * @param panel  the panel.
         */
        public MyChartMouseListener(ChartPanel panel) {
            this.panel = panel;
        }
        /**
         * Callback method for receiving notification of a mouse click on a
         * chart.
         *
         * @param event  information about the event.
         */
        public void chartMouseClicked(ChartMouseEvent e) {
     xy[0]=e.getTrigger().getX();
     xy[1]=e.getTrigger().getY();
      XYPlot plot = (XYPlot) chart.getPlot();
        double ya= plot.getDomainAxis().getLowerBound() ;
        double xa= plot.getRangeAxis().getLowerBound();
        double ye= plot.getDomainAxis().getUpperBound() ;
        double xe= plot.getRangeAxis().getUpperBound();

Rectangle2D plotArea = panel.getScreenDataArea();
//    System.out.println("x="+ plotArea.getMaxX()+" y="+ plotArea.getMaxY());
//    System.out.println("x="+ plotArea.getMinX()+" y="+ plotArea.getMinY());
//    System.out.println("x="+ xa +" y="+ ya);
//    System.out.println("x="+ xe +" y="+ ye);
//    System.out.println("x="+ xy[0] +" y="+ xy[1]);
double chartX=(xy[0]-plotArea.getMinX())/(plotArea.getMaxX()-plotArea.getMinX())*(xe-xa)+xa;
double chartY=-(xy[1]-plotArea.getMinY())/(plotArea.getMaxY()-plotArea.getMinY())*(ye-ya)+ye;
    System.out.println("x="+ chartX +" y="+ chartY);
    }
        /**
         * Callback method for receiving notification of a mouse movement on a
         * chart.
         *
         * @param event  information about the event.
         */
        public void chartMouseMoved(ChartMouseEvent event) {
            // ignore
        }
        

    }

/* static class MyMouseListener implements MouseListener {
 Frame frame;
 public MyMouseListener(Frame frame) {
            this.frame = frame;
        }
public void mouseEntered(MouseEvent e){
if (popup == null) {popup = new Frame();TextArea textArea = new TextArea("Some text to display like a tooltip.");
    frame.add(popup);frame.pack();}else{popup.show();}//System.out.println("Key pressed ");
}
public void mouseReleased(MouseEvent e){}
public void mousePressed(MouseEvent e){}
public void mouseClicked(MouseEvent e){}

public void mouseExited(MouseEvent e){
      if (popup != null) popup.hide();
      }

}
*/

public void windowClosing(WindowEvent e) {
         windowclose();
        }
static   DataOutputStream GnuOutputStream;

static public void windowclose(){
        if(jpgfilename.length()!=0)
         {  BufferedImage image= chart.createBufferedImage(panel.getWidth(),panel.getHeight(),BufferedImage.TYPE_INT_RGB,null);
           try {
                // write the image as a jpg
                ImageIO.write(image,"jpg",new File(jpgfilename));
              } catch(Exception f) {
                f.printStackTrace();
              }
          
                //dispose();
// here we create  gnuplot file
    try{ //Oeffnen der gnu Datei
      System.out.println("Writing "+jpgfilename+".gnu");
     GnuOutputStream = new DataOutputStream(new FileOutputStream(new File(jpgfilename+".gnu")));
     }
   catch (FileNotFoundException e)
      {
      System.out.println("Error opening " + e.getLocalizedMessage());System.exit(0);
       }
         
    
   try{ int w=panel.getWidth();
        int h=panel.getHeight();
       GnuOutputStream.writeBytes("""
set term jpeg enhanced size """+" "+w+" , "+ h + "\n"+ """ 
#set term png enhanced size """+" "+w+" , "+ h + "\n"+ """
#set terminal postscript eps  enhanced color "Arial" 22

set style line 1 lt 1 lw 7 lc rgb "blue" ps 0.3
set style line 2 lt 1 lw 7 lc rgb "red" ps 0.3
set style line 3 lt 1 lw 7 lc rgb "forest-green" ps 0.3
set style line 4 lt 1 lw 7 lc rgb "black" ps 0.3
set style line 5 lt 1 lw 7 lc rgb "magenta" ps 0.3
set style line 6 lt 1 lw 7 lc rgb "orange" ps 0.3
set style line 7 lt 1 lw 7 lc rgb "cyan" ps 0.3
set style line 8 lt 1 lw 7 lc rgb "brown" ps 0.3
set style line 9 lt 2 lw 7 lc rgb "blue" ps 0.3
set key right center
#unset key
#set yr [0.01:20]
#set xr [0.5:1.5]\n """);
GnuOutputStream.writeBytes("set xr ["+chart.getXYPlot().getRangeAxis().getRange().getLowerBound()+":"+
         +chart.getXYPlot().getRangeAxis().getRange().getUpperBound()+"]\n");
GnuOutputStream.writeBytes("set yr ["+chart.getXYPlot().getDomainAxis().getRange().getLowerBound()+":"+
         +chart.getXYPlot().getDomainAxis().getRange().getUpperBound()+"]\n");
GnuOutputStream.writeBytes(gHlines);
GnuOutputStream.writeBytes(gVlines);

 GnuOutputStream.writeBytes("set xlabel '"+chart.getXYPlot().getRangeAxis().getLabel()+"'\n");
 GnuOutputStream.writeBytes("set ylabel '"+chart.getXYPlot().getDomainAxis().getLabel()+"'\n");
 GnuOutputStream.writeBytes("set title '"+chart.getTitle()+"'\n");
 GnuOutputStream.writeBytes("set out '"+jpgfilename+"'\n");
GnuOutputStream.writeBytes("""
#set size ratio 2
#set origin 0, 0
 plot """);
 for(int i=0;i<noffiles;i+=1)
       {
GnuOutputStream.writeBytes(" \""+file[i]+"\" using ");
if(colx[i].contains("c")){GnuOutputStream.writeBytes("("+colx[i].replace("c","$").replace("$os","cos").replace("x","*").replace("e*p","exp")+")");
}else{GnuOutputStream.writeBytes(colx[i]);}
GnuOutputStream.writeBytes(":");
if(coly[i].contains("c")){GnuOutputStream.writeBytes("("+coly[i].replace("c","$").replace("$os","cos").replace("x","*").replace("e*p","exp")+")");
}else{GnuOutputStream.writeBytes(coly[i]);}
GnuOutputStream.writeBytes(" with points ls "+(i+1));
if(i<noffiles-1)GnuOutputStream.writeBytes(", \\\n");
       }
// "results/001mcdisp.qei" using 7:9:(sqrt($10)*2) with points pt 4 pointsize variable , \\
//      "results/002mcdisp.qei" using 7:9:(sqrt($10)*2) with points pt 6 pointsize variable 
GnuOutputStream.writeBytes("""

replot
      """);
       GnuOutputStream.close();


    }
    //Sonstiger Dateifehler
    catch (IOException e)
    { System.out.println("File Error: " + e.getLocalizedMessage());
    }

   } // fi jpgfilename
                System.exit(0);
}
        public void windowOpened(WindowEvent e) {}
        public void windowActivated(WindowEvent e) {}
        public void windowIconified(WindowEvent e) {}
        public void windowDeiconified(WindowEvent e) {}
        public void windowDeactivated(WindowEvent e) {}
        public void windowClosed(WindowEvent e) {}




  public void keyPressed(KeyEvent e) {}
  public void keyReleased(KeyEvent e) {}
 public void keyTyped(KeyEvent e) {
                                    if (e.getKeyChar()=='-'||e.getKeyChar()=='-'){
                                               XYPlot plot = (XYPlot) chart.getPlot();
                                              for (int i=0;i<noffiles;++i)
                                             { if(colyerr[i]=="0"&&colxerr[i]=="0")
                                               { XYErrorRenderer renderer = (XYErrorRenderer) plot.getRenderer(i);
                                                boolean lines=renderer.getSeriesLinesVisible(i);
                                                boolean sym=renderer.getSeriesShapesVisible(i);
                                                if(lines&&sym){lines=false;sym=true;}
                                                else if(lines&&!sym){lines=true;sym=true;}
                                                else if(!lines&&sym){lines=true;sym=false;}

                                                renderer.setSeriesLinesVisible(i,lines);
                                                renderer.setSeriesShapesVisible(i,sym);
                                                update_legend();
                                               }
                                            }}

                                   if (e.getKeyChar()=='S'||e.getKeyChar()=='s'){scale=0.5*scale;for (int i=0;i<noffiles;++i){reload_data(i);};update_legend();}
                                   if (e.getKeyChar()=='B'||e.getKeyChar()=='b'){scale=2*scale;for (int i=0;i<noffiles;++i){reload_data(i);};update_legend();}

//                                    if (e.getKeyChar()=='_'||e.getKeyChar()=='_'){chart.setXAxisVisible(!chart.isXAxisVisible());}
//                                    if (e.getKeyChar()=='|'||e.getKeyChar()=='|'){chart.setYAxisVisible(!chart.isYAxisVisible());}
//                                    if (e.getKeyChar()=='l'||e.getKeyChar()=='L'){chart.setLegendVisible(!chart.isLegendVisible());}
//                                    if (e.getKeyChar()=='s'||e.getKeyChar()=='S'){bRot.setVisible(!bRot.isVisible());}
//                                    if (e.getKeyChar()=='g'||e.getKeyChar()=='G'){chart.getXAxis().setGridVis(!chart.getXAxis().getGridVis());
//                                                                                  chart.getYAxis().setGridVis(!chart.getYAxis().getGridVis());}
//                                    if (e.getKeyChar()=='x'){chart.getXAxis().setLabelPrecision(chart.getXAxis().getLabelPrecision()+1);}
//                                    if (e.getKeyChar()=='X'){chart.getXAxis().setLabelPrecision(chart.getXAxis().getLabelPrecision()-1);}
//                                    if (e.getKeyChar()=='y'){chart.getYAxis().setLabelPrecision(chart.getYAxis().getLabelPrecision()+1);}
//                                    if (e.getKeyChar()=='Y'){chart.getYAxis().setLabelPrecision(chart.getYAxis().getLabelPrecision()-1);}
//                                    if (e.getKeyChar()=='t'){chart.getXAxis().setNumMinTicks(chart.getXAxis().getNumMinTicks()+1);
//                                                             chart.getYAxis().setNumMinTicks(chart.getYAxis().getNumMinTicks()+1);}
//                                    if (e.getKeyChar()=='T'){chart.getXAxis().setNumMinTicks(chart.getXAxis().getNumMinTicks()-1);
//                                                             chart.getYAxis().setNumMinTicks(chart.getYAxis().getNumMinTicks()-1);}
//                                    //if (e.getKeyChar()=='p'){ chart.getXAxis().setLogScaling(!chart.getXAxis().getLogScaling());}
//                                    //if (e.getKeyChar()=='q'){ chart.getYAxis().setLogScaling(!chart.getYAxis().getLogScaling());}


                                  // repaint();
//System.out.println("Key pressed ");
                                   }

      
                                   

  public static void main(String[] args) {
      xmin=1e30;xmax=-1e30;detymin=true;detymax=true;doexit=false;
      ymin=1e30;ymax=-1e30;detxmin=true;detxmax=true;
      detxText=true;detyText=true;detTitle=true;detsTitle=true;detdim=true;logx=false;logy=false;
      prefxsize=500;prefysize=270;
           String ss; String s;
      if (args.length<1)
      {System.out.println("- too few arguments...");
       System.out.println("  program display - show and watch data file by viewing a xy graphic on screen\n");
       System.out.println("use as:  display [-options] xcol[excolerr] ycol[eycolerr][bcolbubble] filename [xcol1[] ycol1[] filename1 ...]\n");
       System.out.println("         xcol,ycol ... column to be taken as x-, y- axis in a lineplot, expressions such as 'c1xc2+1*(c3==2)*(c5<7)' are allowed\n");
       System.out.println("         to plot sum/ product of columns, including math function such as abs,acos,asin,atan,... see complete list in manual ");
       System.out.println("	 filename ..... filename of datafile");
       System.out.println("	 Data files may contain lines to tune the display output, such as");
       System.out.println("	 # displaytitle=My new Graph");
       System.out.println("	 # displayytext=intensity");
       System.out.println("	 # displayxtext=meV ");
       System.out.println("       if optional errorcolumns are added then instead of lines symbols and errorbars are shown");
       System.out.println("	  if optional bubblecolumns are added then instead of lines bubbles with area corresponding to");
       System.out.println("	  bubblecolumn are shown (toggle bubblesize with 's' and 'b' by factor 2)");
       System.out.println("	  (toggle lines also with '-' key))");
//    System.out.println("	 # displaylegend=false (toggle also with 'L' key)\n");
       System.out.println("       options:  -o file.jpg  create a jpg file on exiting, also create file.jpg.gnu to be used in gnuplot");
       System.out.println("                 -c file.jpg  only creates a jpg file and exit immediately, also create file.jpg.gnu ");
       System.out.println("                 -logx -logy  make x(y) a logarithmic axis");
       System.out.println("                 -xmin 23.3 the application sets the minimum of the display xaxis to 23.3");
       System.out.println("                 -xmax -ymin -ymax -xtext -ytext -title -stitle...similar");
       System.out.println("                 -s -l -sl shows symbols/lines/both");
       System.out.println("                 -vlines 2|(201),3.4,12.3 shows vertical lines at specified x values");
       System.out.println("                          a text to be written as line label can be added by inserting | and adding the text");
       System.out.println("                 -hlines 2,3.4,12.3 shows horizontal lines at specified y values");
       System.out.println("                 -g shows gridlines");
       System.out.println("                 -dim 400 200  set dimension of plot (in pixels width 400 height 200)\n");
       System.out.println("                 Press Enter to Continue");

       System.exit(0);
      } scale=0.01;
       file = new String[MAX_NOF_FILES];
       lastmod = new long[MAX_NOF_FILES];
       colx = new String[MAX_NOF_FILES];
       coly = new String[MAX_NOF_FILES];
       colxerr = new String[MAX_NOF_FILES];
       colyerr = new String[MAX_NOF_FILES];
       Double p = Double.valueOf(0.0);
       //      System.out.println(sx+" "+sy);
       //      p.valueOf(strLine);
       //    double[] myDatax = {};
       int j=0;int k=0; jpgfilename="";showgrid=false;showlines=false;showsymbols=true;
       String title="display";
       s=args[0];s=SF.TrimString(s); // command line arguments are treated here
       //look if options are present
       while(SF.TrimString(s).substring(0, 1).equalsIgnoreCase("-"))
          {// yes there are options

           if(SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-o")) // option "-o file.jpg"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             jpgfilename=SF.FirstWord(s);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
           else if(SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-c")) // option "-c file.jpg"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             jpgfilename=SF.FirstWord(s);doexit=true;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-g")) // option "-g"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             showgrid=true;
            }
            else if(SF.TrimString(s).equalsIgnoreCase("-s")) // option "-s"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             showsymbols=true;showlines=false;
            }
            else if(SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-l")) // option "-l"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             showlines=true;showsymbols=false;
            }
            else if(SF.TrimString(s).substring(0, 3).equalsIgnoreCase("-sl")) // option "-sl"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             showlines=true;showsymbols=true;
            }
            else if(SF.TrimString(s).substring(0, 4).equalsIgnoreCase("-dim")) // option "-dim 500 223"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detdim=false;ss=SF.FirstWord(s);prefxsize=p.valueOf(ss).intValue();
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
                          ss=SF.FirstWord(s);prefysize=p.valueOf(ss).intValue();
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-logx")) // option "-logx"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             logx=true;
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-logy")) // option "-logy"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             logy=true;
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-xmin")) // option "-xmin 23"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detxmin=false;ss=SF.FirstWord(s);xmin=p.parseDouble(ss);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-xmax")) // option "-xmax 23"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detxmax=false;ss=SF.FirstWord(s);xmax=p.parseDouble(ss);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-ymin")) // option "-ymin 23"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detymin=false;ss=SF.FirstWord(s);ymin=p.parseDouble(ss);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-ymax")) // option "-ymax 23"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detymax=false;ss=SF.FirstWord(s);ymax=p.parseDouble(ss);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-ytext")) // option "-ytext meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detyText=false;ss=SF.FirstWord(s);yText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-xtext")) // option "-xtext "
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detxText=false;ss=SF.FirstWord(s);xText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-title")) // option "-title text"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detTitle=false;ss=SF.FirstWord(s);Title=ss;Title=Title.replace("_", " ");
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 7).equalsIgnoreCase("-stitle")) // option "-stitle smalltext"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detsTitle=false;ss=SF.FirstWord(s);sTitle=ss;sTitle=sTitle.replace("_", " ");
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 7).equalsIgnoreCase("-hlines")) // option "-hlines 3,2,4"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             ss=SF.FirstWord(s);Hlines=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 7).equalsIgnoreCase("-vlines")) // option "-vlines 3,2,4"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             ss=SF.FirstWord(s);Vlines=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else {break;}
          }

       for(int i=k;s.length()>0;	i+=0)
       {Integer pp;
       ss=SF.FirstWord(s);
       colx[j]=SF.DataCol(ss);       title=title+" "+ss;
       colxerr[j]=SF.ErrorCol(ss);
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       coly[j]=SF.DataCol(ss);       title=title+" "+ss;
       colyerr[j]=SF.ErrorCol(ss);
       if (colyerr[j]=="0") {colyerr[j]=SF.BubbleCol(ss);}
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       file[j]=ss;lastmod[j]=0; title=title+" "+ss;++j;if(j>=MAX_NOF_FILES){System.out.println("ERROR: maximum number of files"+j+" exceeded, recompile with larger MAX_NOF_FILES\n\n");System.exit(0);}
       s=SF.DropWord(s); if (s.length()==0&&i<args.length-1){++i;s=args[i];s=SF.TrimString(s);}
       }noffiles=j;
        display demo = new display(title,prefxsize,prefysize);
       
        
        demo.pack();
       // RefineryUtilities.centerFrameOnScreen(demo);
        demo.setVisible(true);
        final Thread updater = demo.new UpdaterThread();
        //updater.setDaemon(true);
        updater.start();
        Runtime.getRuntime().addShutdownHook(new Thread()
                    {    @Override
                         public void run() { 
                       //  System.out.println("Exiting");
                                           }
                    });
        
      if(doexit==true){windowclose();
                       }
     } //main
// JButton bRot=new JButton("save display.jpg");                       //erstellt einen Button
// Box.Filler bRot1=new Box.Filler (new Dimension(350,10),new Dimension(350,10),new Dimension(370,10));                       //erstellt einen Button
// AbstractButton bRot= new AbstractButton();
 static Expr pclx;
 static Expr pcly;
 static Expr pclxerr;
 static Expr pclyerr;
 static int noffiles;
 static String[] file;
 static String jpgfilename;
 static long[] lastmod;
 static String[] colx;
 static String[] coly;
 static String[] colxerr;
 static String[] colyerr;
 static double scale;
 static double xmin,xmax,ymin,ymax;
 static Integer prefxsize,prefysize;
 static boolean detxmin,detymin,detxmax,detymax,detxText,detyText,detTitle,detsTitle,detdim,doexit;
 static boolean showgrid,showlines,showsymbols,logx,logy;
 static String [] legend; 
 static String yText = "";
 static String xText = "";
 static String Title = "";
 static String sTitle = "";
 static String Vlines = "";
 static String Hlines = "";
 static LegendTitle Legendt;
 static DefaultXYZDataset bdataset;
 static DefaultIntervalXYDataset dataset;
 static JFreeChart chart;
 public ChartPanel chartPanel;
 static ChartPanel panel;
// static JFrame displayFrame;
  
   /**
     * A demonstration application showing a bubble chart.
     *
     * @param title  the frame title.
     */
          
 public display(String title,Integer prefxsize,Integer prefysize) {
        super(title);
        addKeyListener(this); 
        
        //addWindowListener(new MyWindowListener(this,chartPanel));
        //addWindowListener(this);
        //displayFrame= new JFrame();
        dataset = new DefaultIntervalXYDataset();
        JFreeChart chart = createChart(dataset);
        
        ChartPanel chartPanel = new ChartPanel(chart,true,true,true,true,true);
        panel= chartPanel;
       // chartPanel.addChartMouseListener(this);
        chartPanel.addChartMouseListener(new MyChartMouseListener(chartPanel));
       // chartPanel.addMouseListener(new MyMouseListener(this));        
       ToolTipManager.setToolTipText(chartPanel,"Press b/s for bigger/smaller bubbles, - toggles lines/points, use mouse to zoom");
       // chartPanel.addMouseListener(ToolTipManager);
 
        chartPanel.setDomainZoomable(true);
        chartPanel.setRangeZoomable(true);
        //bRot.setHorizontalAlignment(SwingConstants.LEFT);
        //chartPanel.add(bRot);
        chartPanel.setPreferredSize(new java.awt.Dimension(prefxsize, prefysize));
        //chartPanel.setAlignmentX(Component.RIGHT_ALIGNMENT);
        setContentPane(chartPanel);
          // get the top-level container in the Frame (= Window)
        // bRot.setAlignmentY(Component.RIGHT_ALIGNMENT);
        //bRot.setLocation(10,10);
        //bRot.doLayout();
         //bRot.setSize(100,100);
        //bRot.setBounds(10,10,40,40);
        //bRot.setOpaque(true);
        //setLayout(new BorderLayout());
        //add(bRot,BorderLayout.NORTH);
        //add(displayFrame,BorderLayout.SOUTH);
        //setLayout(new FlowLayout(0));
        //setLayout(new CardLayout());
         //bRot.list();
        //add(bRot1);
        //add(bRot);
//   bRot.addActionListener(new ActionListener(){
//    public void actionPerformed(ActionEvent ed){
//    try{
//         FileOutputStream fos=new FileOutputStream("display.jpg");
//         BufferedImage image= chart.createBufferedImage(chartPanel.getWidth(),chartPanel.getHeight(),BufferedImage.TYPE_INT_RGB,null);
//         JPEGImageEncoder encoder= JPEGCodec.createJPEGEncoder(fos);
//         encoder.encode(image);
//         fos.close();
//    }    catch (FileNotFoundException e)
//    {    System.out.println("File not found: " + e.getLocalizedMessage());
         //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
//    }
         //Sonstiger Dateifehler
//         catch (IOException e)
//    {    System.out.println("Dateifehler: " + e.getLocalizedMessage());
         //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
//    }

//      } });

        

                                 
 }// constructor

    /**
     * Creates a chart.
     *
     * @param dataset  the dataset.
     *
     * @return The chart.
     */
    private static JFreeChart createChart(IntervalXYDataset dataset) {
        chart = ChartFactory.createScatterPlot(
                Title, yText, xText, dataset,
                PlotOrientation.HORIZONTAL, true, true, false);
         chart.getTitle().setFont(new Font("SansSerif", Font.PLAIN, 13));
        TextTitle stit=new TextTitle(sTitle);
        stit.setFont(new Font("SansSerif", Font.PLAIN, 10));
        stit.setPosition(RectangleEdge.BOTTOM);
        stit.setHorizontalAlignment(HorizontalAlignment.CENTER);
        chart.addSubtitle(stit);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setForegroundAlpha(1.0f);
plot.setRangeGridlinesVisible(showgrid); 
plot.setRangeGridlinePaint(Color.BLACK);

plot.setDomainGridlinesVisible(showgrid);
plot.setDomainGridlinePaint(Color.BLACK);
       //default to not include zero
        


        bdataset = new DefaultXYZDataset();
       

//        XYBubbleRenderer renderer = ( XYBubbleRenderer)plot.getRenderer();
//    XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
     XYErrorRenderer renderer = new XYErrorRenderer();
     XYBubbleRenderer brenderer = new XYBubbleRenderer(2);
        renderer.setCapLength(0.0);
        renderer.setSeriesPaint(0, Color.blue);
        renderer.setSeriesPaint(1, Color.red);
        renderer.setSeriesPaint(2, Color.green);
        renderer.setSeriesPaint(3, Color.black); // dark
        renderer.setSeriesPaint(4, Color.orange);
        renderer.setSeriesPaint(5, Color.pink);
        
        brenderer.setSeriesPaint(1, Color.blue);
        brenderer.setSeriesPaint(0, Color.red);
        brenderer.setSeriesPaint(3, Color.green);
        brenderer.setSeriesPaint(2, Color.black); // dark
        brenderer.setSeriesPaint(5, Color.orange);
        brenderer.setSeriesPaint(4, Color.pink);
       renderer.setDefaultToolTipGenerator(new StandardXYToolTipGenerator());
       brenderer.setDefaultToolTipGenerator(new StandardXYToolTipGenerator());
    
               
   for(int i=6;i<=MAX_NOF_FILES;++i){ renderer.setSeriesPaint(i, new Color(70*i%256,140*i % 256,210*i % 256));}
   for(int i=6;i<=MAX_NOF_FILES;++i){ brenderer.setSeriesPaint(i, new Color(70*i%256,140*i % 256,210*i % 256));}
           //renderer.setPlotShapes(true);
           //renderer.setShapesFilled(true);
          //renderer.setSeriesShapesVisible(0, true);
          //renderer.setSeriesShapesVisible(1, true);
          //renderer.setSeriesShapesVisible(2, true);
            renderer.setSeriesShape(0, new Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0));
            renderer.setSeriesShape(1, new Rectangle2D.Double(-3.0, -3.0, 6.0, 6.0));
           brenderer.setSeriesShape(0, new Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0));
           brenderer.setSeriesShape(1, new Rectangle2D.Double(-3.0, -3.0, 6.0, 6.0));
           // renderer.setSeriesShape(1, ShapeUtilities.createDiamond(4.0f));

        // increase the margins to account for the fact that the auto-range
        // doesn't take into account the bubble size...
        NumberAxis yAxis = (NumberAxis) plot.getDomainAxis();
        yAxis.setLowerMargin(0.15);
        yAxis.setUpperMargin(0.15);
        NumberAxis xAxis = (NumberAxis) plot.getRangeAxis();
        xAxis.setLowerMargin(0.15);
        xAxis.setUpperMargin(0.15);
        xAxis.setAutoRangeIncludesZero(false);
        yAxis.setAutoRangeIncludesZero(false);
        Double p = Double.valueOf(0.0);
    for(int i=0;i<noffiles;++i){
             if(colyerr[i].startsWith("b")) {    plot.setRenderer(i,brenderer);
                      plot.setDataset(i,bdataset);
                           //            legendItemsNew.add(brenderer.getLegendItem(i,i));
                      }else{plot.setRenderer(i,renderer);
                           plot.setDataset(i,dataset);
                           renderer.setSeriesLinesVisible(i,showlines);
                           renderer.setSeriesShapesVisible(i,showsymbols);
                           }
        

        reload_data(i);
                               }
     if(xmax<xmin||ymax<ymin){System.out.println("No data to plot");System.exit(1);}
     xAxis.setRange(xmin-(xmax-xmin)*0.04,xmax+(xmax-xmin)*0.04);
     yAxis.setRange(ymin-(ymax-ymin)*0.04,ymax+(ymax-ymin)*0.04);
  if(logy){
     LogAxis ylogAxis = new LogAxis(yText);
      if(ymin-(ymax-ymin)*0.04>0)
     {ylogAxis.setRange(ymin-(ymax-ymin)*0.04,ymax+(ymax-ymin)*0.04);}
     plot.setDomainAxis(0,ylogAxis);
        }
  if(logx){
     LogAxis xlogAxis = new LogAxis(xText);
     xlogAxis.setLowerMargin(0.15);
     xlogAxis.setUpperMargin(0.15);
     if(xmin-(xmax-xmin)*0.04>0)
     {xlogAxis.setRange(xmin-(xmax-xmin)*0.04,xmax+(xmax-xmin)*0.04);}
     plot.setRangeAxis(0,xlogAxis);
           }
// this is for plotting a line 
//     XYLineAnnotation axy = new  XYLineAnnotation(0.0, 0.0, 1.0, 0.0);
//     plot.addAnnotation(axy);
// we plot vertical lines at the positions specified in the numbers of string Vlines

    String hl [] = Hlines.split(",");
for (String s : hl) {
if(!s.isEmpty()){
    String sn [] = s.split("\\|"); 
   double y =p.parseDouble(sn[0]); 
 XYLineAnnotation axy = new  XYLineAnnotation(y, xmin, y, xmax);
 gHlines=gHlines+"set object polygon from "+xmin+","+y+" to "+xmax+","+y +" to "+xmin+","+y+ "\n";

plot.addAnnotation(axy);
   if(sn.length>1){
XYTextAnnotation t = new XYTextAnnotation(sn[1],y,xmax+0.02*(xmax-xmin));
gHlines=gHlines+"set label  \""+sn[1]+"\" at "+(xmax+0.02*(xmax-xmin))+","+y+"\n";
plot.addAnnotation(t);
    }
 }
}
    String vl [] = Vlines.split(",");
for (String s : vl) {
if(!s.isEmpty()){ String sn [] = s.split("\\|"); 
   double x =p.parseDouble(sn[0]); 
 XYLineAnnotation axy = new  XYLineAnnotation(ymin, x, ymax, x);
 gVlines=gVlines+"set object polygon from "+x+","+ymin+" to "+x+","+ymax +" to "+x+","+ymin+" \n";

plot.addAnnotation(axy);
 if(sn.length>1){
XYTextAnnotation t = new XYTextAnnotation(sn[1],ymax+0.02*(ymax-ymin),x);
gVlines=gVlines+"set label  \""+sn[1]+"\" at "+x+","+(ymax+0.08*(ymax-ymin))+" center\n";

plot.addAnnotation(t);
    }
 }
}


     update_legend();
     return chart;
    }


      /**
     * A thread for updating the dataset.
     */
    private class UpdaterThread extends Thread {
        /**
         * @see java.lang.Runnable#run()
         */
        public void run() {
            setPriority(MIN_PRIORITY); // be nice
          while(true){
                try {
                    sleep(500);                
      File fileIni;
      for (int i=0;i<noffiles;++i)
           {fileIni = new File(file[i]);
            if(fileIni.lastModified()!=lastmod[i]){lastmod[i]=fileIni.lastModified(); reload_data(i);
            }
           }}
                catch (IndexOutOfBoundsException e) {
                    // suppress
                }
                catch (InterruptedException e) {
                    // suppress
                }
 }}}

protected static void reload_data(int i){    try{
            File fileIni;
            String s="";
            //XYDataset ds = chart.getXYPlot().getDataset(i);
            //ds.getData().removeAllElements();
            int maxnofpoints=10;int j=maxnofpoints;
           while(j==maxnofpoints)
           {double [][] data=new double [6][maxnofpoints];//={{0,1},{0,1},{0,1}};
            double [][] bdata=new double [3][maxnofpoints];
            fileIni = new File(file[i]);
            //?ffnen der Datei
             DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));
             String strLine;
             String sx;
             String sy;
             String sxe;
             String sye;
             String sxerr;
             String syerr;
             String clx = colx[i];
             String cly = coly[i];
             String clxerr = colxerr[i];
             String clyerr = colyerr[i];
           Double p = Double.valueOf(0.0);
             boolean bubbles=false;
             j=0;int dxtf=0; int dytf=0;
          //                System.out.println(clx+" "+cly+" "+clxerr+" "+clyerr);
             if(clyerr.startsWith("b")){clyerr=clyerr.substring(1);bubbles=true;}
             
             if(clx.contains("c")){// parse expression
              try { pclx = Parser.parse(clx.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }
             if(cly.contains("c")){// parse expression
              try { pcly = Parser.parse(cly.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }
             if(clxerr.contains("c")){// parse expression
              try { pclxerr = Parser.parse(clxerr.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }
             if(clyerr.contains("c")){// parse expression
              try { pclyerr = Parser.parse(clyerr.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }

             //Auslesen der Datei
            while (inStream.available() > 0&&j<maxnofpoints)
            {
             strLine = inStream.readLine();
             if (strLine==null) break;
             if (strLine.length() == 0) continue;
      // replace tabs by spaces
      strLine=strLine.replaceAll("[\t\n\u000B\u0009\f]"," ");

// treat comment lines and read variables which might be there to tune plotting
             if(SF.TrimString(strLine).substring(0, 1).equalsIgnoreCase("#"))
             {// remove the comment
              strLine = strLine.substring(1, strLine.length());
      for(int i1=0;i1<=strLine.length();++i1)
       {//if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";chart.addLegend(chart.getXYPlot().Legendt);}}
        //if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";Legendt=chart.getLegend();chart.removeLegend();}}
        if(detxText==true){
            if(i1<=strLine.length()-13&&strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext="))
              {chart.getXYPlot().getRangeAxis().setLabel(strLine.substring(i1+13,strLine.length()));dxtf=1;}
                            }
        if(detyText==true){
            if(i1<=strLine.length()-13&&strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext="))
              {chart.getXYPlot().getDomainAxis().setLabel(strLine.substring(i1+13,strLine.length()));dytf=1;}
                          }
        //if(i1<=strLine.length()-17){if(strLine.substring(i1,i1+17).equalsIgnoreCase("displaylines=true")){chart.setLineVisible(true);}}
        //if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylines=false")){chart.setLineVisible(false);}}
        if(detTitle==true){
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displaytitle=")){chart.setTitle(strLine.substring(i1+13,strLine.length()));}}
                          }
        //if(detsTitle==true){
        //if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+14).equalsIgnoreCase("displaystitle=")){chart.setsubTitle(strLine.substring(i1+14,strLine.length()));}}
        //                  }
        }
        // if no data has yet been read  -go through string and try to find automatically column headers
        if(detxText==true&&dxtf==0&&j==0&&SF.NofCols(strLine)>0)
             if(clx.contains("c")){chart.getXYPlot().getRangeAxis().setLabel(clx);}
             else{chart.getXYPlot().getRangeAxis().setLabel(SF.NthWord(strLine,p.valueOf(clx).intValue()));}
        // if no data has yet been read  -go through string and try to find automatically column headers
        if(detyText==true&&dytf==0&&j==0&&SF.NofCols(strLine)>0)
              if(cly.contains("c")){chart.getXYPlot().getDomainAxis().setLabel(cly);}
             else{chart.getXYPlot().getDomainAxis().setLabel(SF.NthWord(strLine,p.valueOf(cly).intValue()));}
        continue;
             }  // fi is a comment
               // select colx and coly
try{
            if(clx.contains("c")||cly.contains("c")||clxerr.contains("c")||clyerr.contains("c")){     
           Variable [] c=new Variable[SF.NofCols(strLine)+1];
               for(int ii=0;ii<=SF.NofCols(strLine);++ii){c[ii]=Variable.make("c"+ii);
                           if(ii==0){c[ii].setValue(j);}
                          else{String cv=SF.NthWord(strLine,ii);
                                 cv=cv.replace("+-"," ");cv=SF.NthWord(cv,1);cv=cv.replace('D','E');
                                 cv=cv.replaceAll("[a-d,f-z,A-D,F-Z]"," ");cv=SF.NthWord(cv,1);
//System.out.println("c"+ii+"  "+cv);
                              c[ii].setValue(p.parseDouble(cv));
                           }
                               }
                //System.out.println(expr.value());
                if(clx.contains("c")){sx=s.valueOf(pclx.value());}else{sx=SF.NthWord(strLine,p.valueOf(clx).intValue());if(clx=="0"){sx=s.valueOf(j);}}
                if(cly.contains("c")){sy=s.valueOf(pcly.value());}else{sy=SF.NthWord(strLine,p.valueOf(cly).intValue());if(cly=="0"){sy=s.valueOf(j);}}
                if(clxerr.contains("c")){sxe=s.valueOf(pclxerr.value());}else{sxe=SF.NthWord(strLine,p.valueOf(clxerr).intValue());}
                if(clyerr.contains("c")){sye=s.valueOf(pclyerr.value());}else{sye=SF.NthWord(strLine,Math.abs(p.valueOf(clyerr).intValue()));}
               }   
           else{sx=SF.NthWord(strLine,p.valueOf(clx).intValue());if(clx=="0"){sx=s.valueOf(j);}

                 sy=SF.NthWord(strLine,p.valueOf(cly).intValue());if(cly=="0"){sy=s.valueOf(j);}
                 sxe=SF.NthWord(strLine,p.valueOf(clxerr).intValue());
                 sye=SF.NthWord(strLine,Math.abs(p.valueOf(clyerr).intValue()));
                }

            //System.out.println(sx+" "+sy+" "+sxe+" "+sye);

              
   if(sx.length()!=0&&sy.length()!=0&&sxe.length()!=0&&sye.length()!=0
      &&!sx.contains("Infinity")&&!sy.contains("Infinity")&&!sxe.contains("Infinity")&&!sye.contains("Infinity")
      ){
               sx=sx.replace("+-"," ");sx=SF.NthWord(sx,1);
                   sy=sy.replace("+-"," ");sy=SF.NthWord(sy,1);
                    sx=sx.replace('D','E');
                    sy=sy.replace('D','E');
                   sxe=sxe.replace("+-"," "); // if possible take number after +- as error bar
                   sxerr=SF.NthWord(sxe,2);if(sxerr.length()==0){sxerr=SF.NthWord(sxe,1);}
                   sye=sye.replace("+-"," ");
                   syerr=SF.NthWord(sye,2);if(syerr.length()==0){syerr=SF.NthWord(sye,1);}
                    sxerr=sxerr.replace('D','E');
                    syerr=syerr.replace('D','E');
                    if(bubbles)
                   {bdata[1][j]=p.parseDouble(sx);
                     if (detxmin&bdata[1][j]<xmin){xmin=bdata[1][j];}
                     if (detxmax&bdata[1][j]>xmax){xmax=bdata[1][j];}
                    bdata[0][j]=p.parseDouble(sy);
                     if (detymin&bdata[0][j]<ymin){ymin=bdata[0][j];}
                     if (detymax&bdata[0][j]>ymax){ymax=bdata[0][j];}
                    bdata[2][j]=p.parseDouble(syerr);
                    if (bdata[2][j]<0){bdata[2][j]=0;}
                    bdata[2][j]=scale*Math.sqrt(bdata[2][j]);
                  }else
                  {if(clxerr=="0"){sxerr="0";}
                     if(clyerr=="0"){syerr="0";}
                     data[0][j]=p.parseDouble(sy);
                      if (detymin&data[0][j]<ymin){ymin=data[0][j];}
                      if (detymax&data[0][j]>ymax){ymax=data[0][j];}
                     data[1][j]=p.parseDouble(sy)+p.parseDouble(syerr);
                      if (detymax&data[1][j]>ymax){ymax=data[1][j];}
                     data[2][j]=p.parseDouble(sy)-p.parseDouble(syerr);
                      if (detymin&data[2][j]<ymin){ymin=data[2][j];}
                     data[3][j]=p.parseDouble(sx);
                      if (detxmin&data[3][j]<xmin){xmin=data[3][j];}
                      if (detxmax&data[3][j]>xmax){xmax=data[3][j];}
                     data[4][j]=p.parseDouble(sx)+p.parseDouble(sxerr);;
                      if (detxmax&data[4][j]>xmax){xmax=data[4][j];}
                     data[5][j]=p.parseDouble(sx)-p.parseDouble(sxerr);;
                      if (detxmin&data[5][j]<xmin){xmin=data[5][j];}
                   }
                    ++j;
                   }
}
                   catch(NumberFormatException e){if(j>0){--j;}System.exit(1);
                                                  }
                                                          
               }          //         System.out.println(ymin+" "+ymax);

               if(j==maxnofpoints){maxnofpoints*=2;j=maxnofpoints;}
                 else {
               if (j>0)
               {// here fill the rest of the array with the same values
                for(int jj=j;jj<maxnofpoints;++jj)
                  {data[0][jj]=data[0][j-1];data[1][jj]=data[1][j-1];data[2][jj]=data[2][j-1];
                    if(!bubbles){data[3][jj]=data[3][j-1];data[4][jj]=data[4][j-1];data[5][jj]=data[5][j-1];
                                }
                   
                  }

               if(!bubbles)
                   {//dataset.removeSeries(file[i]+s.valueOf(i));
                    dataset.addSeries(file[i]+s.valueOf(i),data);
                    
                   }
                else
                   {//bdataset.removeSeries(file[i]+s.valueOf(i));
                      bdataset.addSeries(file[i]+s.valueOf(i),bdata);
                    }
               }
              }
             }
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
             
    }
    catch(EOFException e)
    {System.out.println("EOF: " + e.getLocalizedMessage());
    }
    catch (FileNotFoundException e)
    {System.out.println("File not found: " + e.getLocalizedMessage());
    }
    //Sonstiger Dateifehler
    catch (IOException e)
    {System.out.println("Dateifehler: " + e.getLocalizedMessage());
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }
}


private static void update_legend () {
XYPlot plot = (XYPlot) chart.getPlot();
LegendItemCollection legendItemsOld = plot.getLegendItems();
final LegendItemCollection legendItemsNew = new LegendItemCollection();

for(int i = 0; i<noffiles&&i<=legendItemsOld.getItemCount(); i++){
    legendItemsNew.add(legendItemsOld.get(i));
}
LegendItemSource source = new LegendItemSource() {
    LegendItemCollection lic = new LegendItemCollection();
    {lic.addAll(legendItemsNew);}
    public LegendItemCollection getLegendItems() {
        return lic;
    }
};
LegendItemSource [] s={source};
chart.getLegend().setSources(s);

//    repaint();

  }




} // display


