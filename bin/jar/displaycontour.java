
/* ---------------------
 * BubbleChartDemo1.java
 * ---------------------
 * (C) Copyright 2003-2008, by Object Refinery Limited.
 */

//package demo;

import java.awt.Color;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.JPanel;
import java.io.*;
import javax.imageio.ImageIO;
//import java.io.*;
//import java.util.EventListener;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.GrayPaintScale;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.data.Range;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
//import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.general.DatasetUtils;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;

import expr.*;

//import com.sun.image.codec.jpeg.JPEGCodec;
//import com.sun.image.codec.jpeg.JPEGImageEncoder;

/**
 * A bubble chart demo.
 */
public class displaycontour extends ApplicationFrame implements WindowListener {
// Button bRot=new Button("save display.jpg");                       //erstellt einen Button
static myStringfunc SF=new myStringfunc();
static final int MAX_NOF_FILES = 1;
static int noffiles;
 static Expr pclx;
 static Expr pcly;
 static Expr pclint;

 static String[] file;
 static String jpgfilename;
 static long[] lastmod;
 static String[] colx;
 static String[] coly;
 static String[] colint;
 static double scale;
 static double bw;
 static double bh;
 static double xmin,xmax,ymin,ymax,zmin,zmax;
 static Integer prefxsize,prefysize;
 static boolean detxmin,detymin,detzmin,detxmax,detymax,detzmax,detxText,detyText,detzText,detTitle,detdim,doexit,showgrid;
 static String [] legend; 
 static String xText = "";
 static String yText = "";
 static String zText = "";
 static String Title = "";
 static String Vlines = "";
 static String gVlines = "";
 static String gHlines = "";
 static String Hlines = "";
 static LegendTitle Legendt;
 static DefaultXYZDataset dataset;
 static NumberAxis zAxis;
 static JFreeChart chart;
 static ChartPanel panel;
 public ChartPanel chartPanel;
 static   DataOutputStream GnuOutputStream;

public void windowClosing(WindowEvent e) {
         windowclose();
        }
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

unset key
#set cbrange [-21:21]
#set cbtics -20, 5, 20
#set contours
#set cntrparam cubic
#set cntrparam levels incremental -20, 5, 20
#set cntrlabel onecolor
#unset colorbox
#unset hidden3d
set palette viridis positive
#set contourfill cbtics
#set pm3d scansauto border retrace
set view map
#set tics scale 0
#set key inside samplen .1 reverse
#set key title "&{----} z = x^2 + y^2(1-x)^3" \n 
#set yr [0.01:20]
#set xr [0.5:1.5]

""");

GnuOutputStream.writeBytes("set xr ["+chart.getXYPlot().getDomainAxis().getRange().getLowerBound()+":"+
         +chart.getXYPlot().getDomainAxis().getRange().getUpperBound()+"]\n");
GnuOutputStream.writeBytes("set yr ["+chart.getXYPlot().getRangeAxis().getRange().getLowerBound()+":"+
         +chart.getXYPlot().getRangeAxis().getRange().getUpperBound()+"]\n");
GnuOutputStream.writeBytes(gHlines);
GnuOutputStream.writeBytes(gVlines);

 GnuOutputStream.writeBytes("set xlabel '"+chart.getXYPlot().getDomainAxis().getLabel()+"'\n");
 GnuOutputStream.writeBytes("set ylabel '"+chart.getXYPlot().getRangeAxis().getLabel()+"'\n");
 GnuOutputStream.writeBytes("set title '"+chart.getTitle()+"'\n");
 GnuOutputStream.writeBytes("set out '"+jpgfilename+"'\n");
GnuOutputStream.writeBytes("""
#set size ratio 2
#set origin 0, 0
#splot g(x,y) with contourfill fs solid border notitle, \\
#      g(x,y) nosurface lt black title "Contour levels dz = 5"
# set dgrid3d  splines
 splot """);
 for(int i=0;i<noffiles;i+=1)
       {
GnuOutputStream.writeBytes(" \""+file[i]+"\" using ");
if(colx[i].contains("c")){GnuOutputStream.writeBytes("("+colx[i].replace("c","$").replace("$os","cos").replace("x","*").replace("e*p","exp")+")");
}else{GnuOutputStream.writeBytes(colx[i]);}
GnuOutputStream.writeBytes(":");
if(coly[i].contains("c")){GnuOutputStream.writeBytes("("+coly[i].replace("c","$").replace("$os","cos").replace("x","*").replace("e*p","exp")+")");
}else{GnuOutputStream.writeBytes(coly[i]);}
GnuOutputStream.writeBytes(":");
if(colint[i].contains("c")){GnuOutputStream.writeBytes("("+colint[i].replace("c","$").replace("$os","cos").replace("x","*").replace("e*p","exp")+")");
}else{GnuOutputStream.writeBytes(colint[i]);}
GnuOutputStream.writeBytes("  with pm3d  \n");
if(i<noffiles-1)GnuOutputStream.writeBytes(", \\\n");
       }
// "results/001mcdisp.qei" using 7:9:(sqrt($10)*2) with points pt 4 pointsize variable , \\
//      "results/002mcdisp.qei" using 7:9:(sqrt($10)*2) with points pt 6 pointsize variable 
 for(int i=0;i<noffiles;i+=1)
       {
GnuOutputStream.writeBytes("""
# if 'with pm3d' does not work you can try 'with points pointtype 5  palette'
# or to insert empty lines between changes of data column use the command
 """);
GnuOutputStream.writeBytes("# comment -cc ");
if(colx[i].contains("c")){GnuOutputStream.writeBytes("("+colx[i].replace("c","$").replace("$os","cos").replace("x","*").replace("e*p","exp")+")");
}else{GnuOutputStream.writeBytes(colx[i]);}
GnuOutputStream.writeBytes(" \" \" "+file[i]+"\n");
   }
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

    /**
     * A demonstration application showing a bubble chart.
     *
     * @param title  the frame title.
     */
    public displaycontour(String title,Integer prefxsize,Integer prefysize) {
        super(title);


        dataset = new DefaultXYZDataset();
        JFreeChart chart=createChart(dataset);

        //chartPanel = createDemoPanel();
        ChartPanel chartPanel = new ChartPanel(chart,true,true,true,true,true);
        
        panel= chartPanel;
        chartPanel.setPreferredSize(new java.awt.Dimension(prefxsize, prefysize));
        setContentPane(chartPanel);
        chartPanel.setAlignmentY(Component.LEFT_ALIGNMENT);
//        chartPanel.add(bRot);

//   bRot.addActionListener(new ActionListener(){
//    public void actionPerformed(ActionEvent ed){
//    try{
//         FileOutputStream fos=new FileOutputStream("display.jpg");
//         BufferedImage image= chart.createBufferedImage(chartPanel.getWidth(),chartPanel.getHeight(),BufferedImage.TYPE_INT_RGB,null); 
//         JPEGImageEncoder encoder= JPEGCodec.createJPEGEncoder(fos); 
//         encoder.encode(image);
//         fos.close();
//    }    catch (FileNotFoundException e)
//    {
//         System.out.println("File not found: " + e.getLocalizedMessage());
//         //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
//    }
//         //Sonstiger Dateifehler
//         catch (IOException e)
//    {
//         System.out.println("Dateifehler: " + e.getLocalizedMessage());
//         //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
//    }

//      } });
                                               
    }

    /**
     * Creates a chart.
     *
     * @param dataset  the dataset.
     *
     * @return The chart.
     */
    private  JFreeChart createChart(DefaultXYZDataset dataset) {
        NumberAxis xAxis = new NumberAxis(xText);
//         xAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
         xAxis.setLowerMargin(0.05);
         xAxis.setUpperMargin(0.05);
         NumberAxis yAxis = new NumberAxis(yText);
//         yAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
         yAxis.setLowerMargin(0.05);
         yAxis.setUpperMargin(0.05);
         zAxis = new NumberAxis(zText);
        chart = ChartFactory.createScatterPlot(
                Title, xText, yText, dataset,
                PlotOrientation.VERTICAL, false, false, false);

         XYBlockRenderer renderer = new XYBlockRenderer();
          XYPlot plot = (XYPlot) chart.getPlot();
//         XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
//         plot.setBackgroundPaint(Color.lightGray);
         plot.setDomainGridlinesVisible(false);
         plot.setRangeGridlinePaint(Color.white);
         chart = new JFreeChart(Title, plot);
         chart.removeLegend();
         chart.setBackgroundPaint(Color.white);
for(int i=0;i<noffiles;++i){
                           plot.setRenderer(i,renderer);
                           plot.setDataset(i,dataset);
         reload_data(i);
                               }
int[] red =
    {  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  7, 23, 39, 55, 71, 87,103,
       119,135,151,167,183,199,215,231,
       247,255,255,255,255,255,255,255,
       255,255,255,255,255,255,255,255,
       255,246,228,211,193,175,158,140};
    int[] green =
    {  0,  0,  0,  0,  0,  0,  0,  0,
       0, 11, 27, 43, 59, 75, 91,107,
       123,139,155,171,187,203,219,235,
       251,255,255,255,255,255,255,255,
       255,255,255,255,255,255,255,255,
       255,247,231,215,199,183,167,151,
       135,119,103, 87, 71, 55, 39, 23,
       7,  0,  0,  0,  0,  0,  0,  0};
    int[] blue =
    {  0,143,159,175,191,207,223,239,
       255,255,255,255,255,255,255,255,
       255,255,255,255,255,255,255,255,
       255,247,231,215,199,183,167,151,
       135,119,103, 87, 71, 55, 39, 23,
       7,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0};

     if(xmax<=xmin||ymax<=ymin||zmax<=zmin){System.out.println("No data to plot");System.exit(1);}
      xAxis.setRangeWithMargins(new Range(xmin-(xmax-xmin)*0.04,xmax+(xmax-xmin)*0.04),true,true);
      yAxis.setRangeWithMargins(new Range(ymin-(ymax-ymin)*0.04,ymax+(ymax-ymin)*0.04),true,true);
         LookupPaintScale scale = new LookupPaintScale(zmin, zmax,Color.red);
         for(int i=0;i<64;++i){double value=zmin+i*(zmax-zmin)/64;//System.out.println(value);
         scale.add(value,new Color(red[i],green[i],blue[i]));
                                  }
         PaintScaleLegend zscale = new PaintScaleLegend(scale,zAxis);
        zscale.setVisible(true);
          chart.addSubtitle(zscale);
         renderer.setPaintScale(scale);
        renderer.setBlockWidth(bw);
        renderer.setBlockHeight(bh);
        plot.setBackgroundPaint(Color.white);
        plot.setForegroundAlpha(1.0f);
plot.setRangeGridlinesVisible(showgrid);
plot.setRangeGridlinePaint(Color.WHITE);

plot.setDomainGridlinesVisible(showgrid);
plot.setDomainGridlinePaint(Color.WHITE);
//System.out.println(DatasetUtils.findDomainBounds(dataset, false).getLowerBound());
//System.out.println(renderer.findRangeBounds(dataset));

// this is for plotting a line 
//     XYLineAnnotation axy = new  XYLineAnnotation(0.0, 0.0, 1.0, 0.0);
//     plot.addAnnotation(axy);
// we plot vertical lines at the positions specified in the numbers of string Vlines

    String vl [] = Vlines.split(",");Double p = new Double(0.0);
for (String s : vl) {
if(!s.isEmpty()){
    String sn [] = s.split("\\|"); 
   double x =p.parseDouble(sn[0]); 
 XYLineAnnotation axy = new  XYLineAnnotation(x, ymin, x, ymax+0.01*(ymax-ymin));
 gVlines=gVlines+"set object polygon from "+x+","+ymin+" to "+x+","+(ymax+0.01*(ymax-ymin)) +" to "+x+","+ymin+"  \n";
plot.addAnnotation(axy);
   if(sn.length>1){
XYTextAnnotation t = new XYTextAnnotation(sn[1],x,ymax+0.02*(ymax-ymin));
gVlines=gVlines+"set label  \""+sn[1]+"\" at "+x+","+(ymax+0.03*(ymax-ymin))+" center\n";
plot.addAnnotation(t);
    }
 }
}
    String hl [] = Hlines.split(",");
for (String s : hl) {
if(!s.isEmpty()){ String sn [] = s.split("\\|"); 
   double y =p.parseDouble(sn[0]); 
 XYLineAnnotation axy = new  XYLineAnnotation(xmin, y, xmax+0.01*(xmax-xmin), y);
 gHlines=gHlines+"set object polygon from "+xmin+","+y+" to "+(xmax+0.01*(xmax-xmin))+","+y +" to "+xmin+","+y+" \n";
plot.addAnnotation(axy);
 if(sn.length>1){
XYTextAnnotation t = new XYTextAnnotation(sn[1],xmax+0.02*(xmax-xmin),y);
gHlines=gHlines+"set label  \""+sn[1]+"\" at "+(xmax+0.02*(xmax-xmin))+","+y+"\n";

plot.addAnnotation(t);
    }
 }
}
         return chart;
    }

    /**
     * Creates a sample dataset.
     *
     * @return A sample dataset.
     */
  /*  public static XYZDataset createDataset() {
        
         dataset = new DefaultXYZDataset();
        //double[] x = {2.1, 2.3, 2.3, 2.2, 2.2, 1.8, 1.8, 1.9, 2.3, 3.8};
        //double[] y = {14.1, 11.1, 10.0, 8.8, 8.7, 8.4, 5.4, 4.1, 4.1, 25};
        //double[] z = {2.4, 2.7, 2.7, 2.2, 2.2, 2.2, 2.1, 2.2, 1.6, 4};
        //double[][] series = new double[][] { x, y, z };
        //dataset.addSeries("Series 1", series);
        return dataset;
    }*/

    /**
     * Creates a panel for the demo (used by SuperDemo.java).
     *
     * @return A panel.
     */
/*    public  ChartPanel createDemoPanel() {
        JFreeChart chart = createChart();
        ChartPanel chartPanel = new ChartPanel(chart);
       	
        chartPanel.setDomainZoomable(true);
        chartPanel.setRangeZoomable(true);
        return chartPanel;
    }
*/
    /**
     * Starting point for the demonstration application.
     *
     * @param args  ignored.
     */
    public static void main(String[] args) {
xmin=1e30;xmax=-1e30;detymin=true;detymax=true;doexit=false;
ymin=1e30;ymax=-1e30;detxmin=true;detxmax=true;
zmin=1e30;zmax=-1e30;detzmin=true;detzmax=true;
      detxText=true;detyText=true;detzText=true;detTitle=true;detdim=true;
      prefxsize=500;prefysize=270;
     
          String ss,s;
      if (args.length<3)
      {System.out.println("- too few arguments...\n");
       System.out.println("  program displaycontour - show and watch data file by viewing a xy graphic on screen\n");
       System.out.println("use as:  displaycontour [options] xcol ycol intcol filename \n");
       System.out.println("         xcol,ycol,intcol ... column to be taken as x-, y- and intensity-axis");
       System.out.println("                              x-spacings have to be equal, also y-spacings may not vary ");
       System.out.println("	 filename ..... filename of datafile");
       System.out.println("	 Data files may contain lines to tune the display output, such as");
       System.out.println("	 # displaytitle=My new Graph");
       System.out.println("	 # displayytext=intensity");
       System.out.println("	 # displayxtext=meV \n");
       System.out.println("        options:   -o file.jpg create a jpg file on exiting, also create file.jpg.gnu to be used in gnuplot");
       System.out.println("                   -c file.jpg create a jpg file and exit immediately, also create file.jpg.gnu ");
       System.out.println("                   -xmin 23.3 the application sets the minimum of the display xaxis to 23.3");
       System.out.println("                   -xmax -ymin -ymax -xtext -ytext -title  similar");
       System.out.println("                   -vlines 2|(201),3.4,12.3 shows vertical lines at specified x values");
       System.out.println("                          a text to be written as line label can be added by inserting | and adding the text");
       System.out.println("                   -hlines 2,3.4,12.3 shows horizontal lines at specified y values");
       System.out.println("                   -g shows gridlines");
       System.out.println("                   -dim 400 200  set dimension of plot (in pixels width 400 height 200)\n");
       System.out.println("                 Press Enter to Continue");
       System.exit(0);
      }
       file = new String[args.length/3];
       lastmod = new long[args.length/3];
       colx = new String[args.length/3];
       coly = new String[args.length/3];
       colint = new String[args.length/3];
       Double p = new Double(0.0);
       //      System.out.println(sx+" "+sy);
       //      p.valueOf(strLine);
       //    double[] myDatax = {};
  
int j=0;int k=0; jpgfilename="";showgrid=false;
       String title="displaycontour";
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
            else if(SF.TrimString(s).substring(0, 4).equalsIgnoreCase("-dim")) // option "-dim 500 223"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detdim=false;ss=SF.FirstWord(s);prefxsize=p.valueOf(ss).intValue();
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
                          ss=SF.FirstWord(s);prefysize=p.valueOf(ss).intValue();
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-zmax")) // option "-zmax 23"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detzmax=false;ss=SF.FirstWord(s);zmax=p.parseDouble(ss);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 5).equalsIgnoreCase("-zmin")) // option "-zmin 23"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detzmin=false;ss=SF.FirstWord(s);zmin=p.parseDouble(ss);
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
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
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-ztext")) // option "-ztext meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detzText=false;ss=SF.FirstWord(s);zText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-ytext")) // option "-ytext meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detyText=false;ss=SF.FirstWord(s);yText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-xtext")) // option "-xtext meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detxText=false;ss=SF.FirstWord(s);xText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-title")) // option "-title meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detTitle=false;ss=SF.FirstWord(s);Title=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 7).equalsIgnoreCase("-hlines")) // option "-hlines 3,2,4"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detTitle=false;ss=SF.FirstWord(s);Hlines=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 7).equalsIgnoreCase("-vlines")) // option "-vlines 3,2,4"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detTitle=false;ss=SF.FirstWord(s);Vlines=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else {System.out.println("ERROR: option,"+SF.TrimString(s)+" not implemented !\n\n");System.exit(0);}
          }
       for(int i=k;s.length()>0;	i+=0)
       {Integer pp;
       ss=SF.FirstWord(s);
       colx[j]=SF.DataCol(ss);       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       coly[j]=SF.DataCol(ss);       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       colint[j]=SF.DataCol(ss);       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       file[j]=ss;lastmod[j]=0; title=title+" "+ss;++j;if(j>MAX_NOF_FILES){System.out.println("ERROR: maximum number of files"+j+" exceeded, recompile with larger MAX_NOF_FILES\n\n");System.exit(0);}
       s=SF.DropWord(s); if (s.length()==0&&i<args.length-1){++i;s=args[i];s=SF.TrimString(s);}
       }noffiles=j;
//       System.out.println(xText);
        displaycontour demo = new displaycontour(title,prefxsize,prefysize);
        demo.pack();
        //RefineryUtilities.centerFrameOnScreen(demo);
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

    } // main




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
           for (int i=0;i<1;++i)
                  {fileIni = new File(file[i]);
                   if(fileIni.lastModified()!=lastmod[i]){lastmod[i]=fileIni.lastModified(); reload_data(i);}
                  }
           }
 catch (IndexOutOfBoundsException e) {
                    // suppress
                }
                catch (InterruptedException e) {
                    // suppress
                }
}}}

protected void reload_data(int i)
      {  File fileIni;
            String s="";
       try{
            //XYDataset ds = chart.getXYPlot().getDataset(i);
            //ds.getData().removeAllElements();
            int maxnofpoints=1000;int j=maxnofpoints;double xold=0,yold=0;bh=0;bw=0;
           while(j==maxnofpoints)           
           {double [][] data=new double [3][maxnofpoints];//={{0,1},{0,1},{0,1}};             

            fileIni = new File(file[i]);
            //?ffnen der Datei
             DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));
             String strLine;
             String sx;
             String sy;
             String sint;
             String clx = colx[i];
             String cly = coly[i];   
             String clint = colint[i];
            Double p = Double.valueOf(0.0);
             if(clx.contains("c")){// parse expression
              try { pclx = Parser.parse(clx.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }
             if(cly.contains("c")){// parse expression
              try { pcly = Parser.parse(cly.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }
            if(clint.contains("c")){// parse expression
              try { pclint = Parser.parse(clint.replace("x","*").replace("e*p","exp")); } catch (SyntaxException e) { System.err.println(e.explain()); System.exit(1); }
                                 }
            
             j=0;int dxtf=0; int dytf=0;int dztf=0;
             //Auslesen der Datei
            while (inStream.available() > 0&&j<maxnofpoints)
            {
             strLine = inStream.readLine();
             if (strLine==null) break;
             if (strLine.length() == 0) continue;
      // replace tabs by spaces
      strLine=strLine.replaceAll("[\t\n\u000B\u0009\f]"," ");

             if(SF.TrimString(strLine).substring(0, 1).equalsIgnoreCase("#"))
             {// remove the comment
              strLine = strLine.substring(1, strLine.length());
      for(int i1=0;i1<=strLine.length();++i1)
       {//if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";chart.addLegend(chart.getXYPlot().Legendt);}}
        //if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";Legendt=chart.getLegend();chart.removeLegend();}}
        if(detxText==true){
            if(i1<=strLine.length()-13&&strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext="))
              {chart.getXYPlot().getDomainAxis().setLabel(strLine.substring(i1+13,strLine.length()));dxtf=1;}
                           }
        if(detyText==true){
            if(i1<=strLine.length()-13&&strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext="))
              {chart.getXYPlot().getRangeAxis().setLabel(strLine.substring(i1+13,strLine.length()));dytf=1;}
                           }
        if(detzText==true){
            if(i1<=strLine.length()-13&&strLine.substring(i1,i1+13).equalsIgnoreCase("displayztext="))
              {zAxis.setLabel(strLine.substring(i1+13,strLine.length()));dztf=1;}
                           }
        //if(i1<=strLine.length()-17){if(strLine.substring(i1,i1+17).equalsIgnoreCase("displaylines=true")){chart.setLineVisible(true);}}
        //if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylines=false")){chart.setLineVisible(false);}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displaytitle=")){chart.setTitle(strLine.substring(i1+13,strLine.length()));}}
        }
         // if no data has yet been read  -go through string and try to find automatically column headers
          if(detxText==true&&dxtf==0&&j==0&&SF.NofCols(strLine)>0)
                  if(clx.contains("c")){chart.getXYPlot().getDomainAxis().setLabel(clx);}
                  else{if(SF.NofCols(strLine)>p.valueOf(clx).intValue())chart.getXYPlot().getDomainAxis().setLabel(SF.NthWord(strLine,p.valueOf(clx).intValue()));}
          if(detyText==true&&dytf==0&&j==0&&SF.NofCols(strLine)>0)
                 if(cly.contains("c")){chart.getXYPlot().getRangeAxis().setLabel(cly);}
                 else{if(SF.NofCols(strLine)>p.valueOf(cly).intValue())chart.getXYPlot().getRangeAxis().setLabel(SF.NthWord(strLine,p.valueOf(cly).intValue()));}
          if(detzText==true&&dztf==0&&j==0&&SF.NofCols(strLine)>0)
                 if(clint.contains("c")){zAxis.setLabel(clint);}
                 else{if(SF.NofCols(strLine)>p.valueOf(clint).intValue())zAxis.setLabel(SF.NthWord(strLine,p.valueOf(clint).intValue()));}
//System.out.println(SF.NthWord(strLine,p.valueOf(clint).intValue())+"ddd"+zAxis.getLabel());
        continue;
             }  // fi is a comment
            
          // select colx and coly
            if(clx.contains("c")||cly.contains("c")||clint.contains("c")){     
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
                if(clx.contains("c")){sx=s.valueOf(pclx.value());}else{sx=SF.NthWord(strLine,p.valueOf(clx).intValue());if(clx=="0"){sx=s.valueOf(j);}}
                if(cly.contains("c")){sy=s.valueOf(pcly.value());}else{sy=SF.NthWord(strLine,p.valueOf(cly).intValue());if(cly=="0"){sy=s.valueOf(j);}}
                if(clint.contains("c")){sint=s.valueOf(pclint.value());}else{sint=SF.NthWord(strLine,p.valueOf(clint).intValue());if(clint=="0"){sint=s.valueOf(j);}}
                  }   
           else{sx=SF.NthWord(strLine,p.valueOf(clx).intValue());if(clx=="0"){sx=s.valueOf(j);}
                 sy=SF.NthWord(strLine,p.valueOf(cly).intValue());if(cly=="0"){sy=s.valueOf(j);}
                 sint=SF.NthWord(strLine,p.valueOf(clint).intValue());if(clint=="0"){sint=s.valueOf(j);}
                  }
//             System.out.println(sx+" "+sy+" "+clx+" "+cly);

              
   if(sx.length()!=0&&sy.length()!=0&&sint.length()!=0){
               try{sx=sx.replace("+-"," ");sx=SF.NthWord(sx,1);
                   sy=sy.replace("+-"," ");sy=SF.NthWord(sy,1);
                   sint=sint.replace("+-"," ");sint=SF.NthWord(sint,1);
                    sx=sx.replace('D','E');
                    sy=sy.replace('D','E');
                    sint=sint.replace('D','E');
                     data[0][j]=p.parseDouble(sx);if(j>0){double bwg=Math.abs(data[0][j]-xold);if(bwg>0&&(bwg<bw||bw==0)){bw=bwg;}}
                      xold=data[0][j];
                      if (detxmin&data[0][j]<xmin){xmin=data[0][j];}
                      if (detxmax&data[0][j]>xmax){xmax=data[0][j];}
                     data[1][j]=p.parseDouble(sy);if(j>0){double bhg=Math.abs(data[1][j]-yold);if(bhg>0&&(bhg<bh||bh==0)){bh=bhg;}}
                     yold=data[1][j];
                      if (detymin&data[1][j]<ymin){ymin=data[1][j];}
                      if (detymax&data[1][j]>ymax){ymax=data[1][j];}
                     data[2][j]=p.parseDouble(sint);
                      if (detzmin&data[2][j]<zmin){zmin=data[2][j];}
                      if (detzmax&data[2][j]>zmax){zmax=data[2][j];}
                    ++j;
                   }
                   catch(NumberFormatException e){if(j>0){--j;}//System.exit(1);
                                                  }
                                                          }
               }                
      // System.out.println("x:"+xmin+" "+xmax);
      // System.out.println("y:"+ymin+" "+ymax);
      // System.out.println("z:"+zmin+" "+zmax);
      chart.getXYPlot().getDomainAxis().setRangeWithMargins(new Range(xmin-(xmax-xmin)*0.04,xmax+(xmax-xmin)*0.04),true,true);
      chart.getXYPlot().getRangeAxis().setRangeWithMargins(new Range(ymin-(ymax-ymin)*0.04,ymax+(ymax-ymin)*0.04),true,true);
      
               if(j==maxnofpoints){maxnofpoints*=2;j=maxnofpoints;}
                 else {
               if (j>0)
               {double [][] series=new double [3][j];
                // here fill the rest of the array with the same values
                for(int jj=0;jj<j;++jj)
                  {series[0][jj]=data[0][jj];series[1][jj]=data[1][jj];series[2][jj]=data[2][jj];
                  }

                    dataset.addSeries(file[i]+s.valueOf(i),series);
                   
               }
              }
             
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
   }}
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
    }

    catch (FileNotFoundException e)
    {
      System.out.println("File not found: " + e.getLocalizedMessage());
    }

    //Sonstiger Dateifehler
    catch (IOException e)
    {
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }
//    repaint();
//System.out.println("Displaycontour: Data reloaded ymax:"+ymax);
  }
    
 

}