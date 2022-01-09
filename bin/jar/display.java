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
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.DefaultIntervalXYDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.ui.RectangleEdge;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.labels.StandardXYToolTipGenerator;

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
static Frame frame;
static Frame popup;
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
static public void windowclose(){
        if(jpgfilename.length()!=0)
         {  BufferedImage image= chart.createBufferedImage(panel.getWidth(),panel.getHeight(),BufferedImage.TYPE_INT_RGB,null);
           try {
                // write the image as a jpg
                ImageIO.write(image,"jpg",new File(jpgfilename));
              } catch(Exception f) {
                f.printStackTrace();
              }
          }
                //dispose();
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
                                             { if(colyerr[i]==0&&colxerr[i]==0)
                                               { XYErrorRenderer renderer = (XYErrorRenderer) plot.getRenderer(i);
                                                renderer.setSeriesLinesVisible(i,!renderer.getSeriesLinesVisible(i));
                                                renderer.setSeriesShapesVisible(i,!renderer.getSeriesShapesVisible(i));
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
      detxText=true;detyText=true;detTitle=true;
           String ss; String s;
      if (args.length<1)
      {System.out.println("- too few arguments...\n");
       System.out.println("  program display - show and watch data file by viewing a xy graphic on screen\n\n");
       System.out.println("use as:  display [-options] xcol[excolerr] ycol[eycolerr][bcolbubble] filename [xcol1[] ycol1[] filename1 ...]\n\n");
       System.out.println("         xcol,ycol ... column to be taken as x-, y- axis\n in a lineplot");
       System.out.println("        when using option -o file.jpg the application creates a jpg file on exiting\n");
       System.out.println("        when using option -xmin 23.3 the application sets the minimum of the display xaxis to 23.3\n");
       System.out.println("        similar are options -xmax -ymin -ymax -xtext -ytext -title....\n");
       System.out.println("        when using option -c file.jpg the application only creates a jpg file and exits immediatly\n");
       System.out.println("        if optional errorcolumns are added then instead of lines symbols and errorbars are shown\n");
       System.out.println("	  if optional bubblecolumns are added then instead of lines bubbles with area corresponding to\n");
       System.out.println("	  bubblecolumn are shown (toggle bubblesize with 's' and 'b')\n");
       System.out.println("	  (toggle lines also with '-' key))\n");
       System.out.println("	 filename ..... filename of datafile\n\n");
    System.out.println("	 Data files may contain lines to tune the display output, such as\n");
    System.out.println("	 # displaytitle=My new Graph\n");
    System.out.println("	 # displayytext=intensity\n");
    System.out.println("	 # displayxtext=meV \n");
//    System.out.println("	 # displaylegend=false (toggle also with 'L' key)\n\n");
       System.exit(0);
      } scale=0.01;
       file = new String[MAX_NOF_FILES];
       lastmod = new long[MAX_NOF_FILES];
       colx = new int[MAX_NOF_FILES];
       coly = new int[MAX_NOF_FILES];
       colxerr = new int[MAX_NOF_FILES];
       colyerr = new int[MAX_NOF_FILES];
       Double p = new Double(0.0);
       //      System.out.println(sx+" "+sy);
       //      p.valueOf(strLine);
       //    double[] myDatax = {};
       int j=0;int k=0; jpgfilename="";
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
             detyText=false;ss=SF.FirstWord(s);xText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-xtext")) // option "-xtext meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detxText=false;ss=SF.FirstWord(s);yText=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else if(SF.TrimString(s).substring(0, 6).equalsIgnoreCase("-title")) // option "-title meV"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
             detTitle=false;ss=SF.FirstWord(s);Title=ss;
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
            else {System.out.println("ERROR: option,"+SF.TrimString(s)+" not implemented !\n\n");System.exit(0);}
          }
       for(int i=k;s.length()>0;	i+=0)
       {Integer pp;
       ss=SF.FirstWord(s);
       colx[j]=p.valueOf(SF.DataCol(ss)).intValue();       title=title+" "+ss;
       colxerr[j]=p.valueOf(SF.ErrorCol(ss)).intValue();
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       coly[j]=p.valueOf(SF.DataCol(ss)).intValue();       title=title+" "+ss;
       colyerr[j]=p.valueOf(SF.ErrorCol(ss)).intValue();
       if (colyerr[j]==0) {colyerr[j]=-p.valueOf(SF.BubbleCol(ss)).intValue();}
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       file[j]=ss;lastmod[j]=0; title=title+" "+ss;++j;if(j>=MAX_NOF_FILES){System.out.println("ERROR: maximum number of files"+j+" exceeded, recompile with larger MAX_NOF_FILES\n\n");System.exit(0);}
       s=SF.DropWord(s); if (s.length()==0&&i<args.length-1){++i;s=args[i];s=SF.TrimString(s);}
       }noffiles=j;
        display demo = new display(title);
       
        
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
 static int noffiles;
 static String[] file;
 static String jpgfilename;
 static long[] lastmod;
 static int[] colx;
 static int[] coly;
 static int[] colxerr;
 static int[] colyerr;
 static double scale;
 static double xmin,xmax,ymin,ymax;
 static boolean detxmin,detymin,detxmax,detymax,detxText,detyText,detTitle,doexit;
 static String [] legend; 
 static String xText = "";
 static String yText = "";
 static String Title = "";
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
          
 public display(String title) {
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
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
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
                Title, xText, yText, dataset,
                PlotOrientation.HORIZONTAL, true, true, false);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setForegroundAlpha(1.0f);
      //  plot.setDomainGridlinesVisible(true);
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
        renderer.setSeriesPaint(3, Color.black);
        renderer.setSeriesPaint(4, Color.orange);
        renderer.setSeriesPaint(5, Color.pink);

        brenderer.setSeriesPaint(1, Color.blue);
        brenderer.setSeriesPaint(0, Color.red);
        brenderer.setSeriesPaint(3, Color.green);
        brenderer.setSeriesPaint(2, Color.black);
        brenderer.setSeriesPaint(5, Color.orange);
        brenderer.setSeriesPaint(4, Color.pink);
       renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
       brenderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
    

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
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setLowerMargin(0.15);
        domainAxis.setUpperMargin(0.15);
        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setLowerMargin(0.15);
        rangeAxis.setUpperMargin(0.15);
        rangeAxis.setAutoRangeIncludesZero(false);
        domainAxis.setAutoRangeIncludesZero(false);
        
    for(int i=0;i<noffiles;++i){
         if(colyerr[i]>=0){plot.setRenderer(i,renderer);
                           plot.setDataset(i,dataset);
                           renderer.setSeriesLinesVisible(i,false);
                           renderer.setSeriesShapesVisible(i, true);
                           }
            else {    plot.setRenderer(i,brenderer);
                      plot.setDataset(i,bdataset);
                           //            legendItemsNew.add(brenderer.getLegendItem(i,i));
                      }

        reload_data(i);
                               }
     rangeAxis.setRange(xmin-(xmax-xmin)*0.04,xmax+(xmax-xmin)*0.04);
     domainAxis.setRange(ymin-(ymax-ymin)*0.04,ymax+(ymax-ymin)*0.04);

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
             String sxerr;
             String syerr;
             int clx = colx[i];
             int cly = coly[i];
             int clxerr = colxerr[i];
             int clyerr = colyerr[i];

             j=0;
             //Auslesen der Datei
            while (inStream.available() > 0&&j<maxnofpoints)
            {
             strLine = inStream.readLine();
             if (strLine==null) break;
             if (strLine.length() == 0) continue;
             if(SF.TrimString(strLine).substring(0, 1).equalsIgnoreCase("#"))
             {
      for(int i1=0;i1<=strLine.length();++i1)
       {//if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";chart.addLegend(chart.getXYPlot().Legendt);}}
        //if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";Legendt=chart.getLegend();chart.removeLegend();}}
        if(detyText==true&&i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext=")){chart.getXYPlot().getRangeAxis().setLabel(strLine.substring(i1+13,strLine.length()));}}
        if(detxText==true&&i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext=")){chart.getXYPlot().getDomainAxis().setLabel(strLine.substring(i1+13,strLine.length()));}}
        //if(i1<=strLine.length()-17){if(strLine.substring(i1,i1+17).equalsIgnoreCase("displaylines=true")){chart.setLineVisible(true);}}
        //if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylines=false")){chart.setLineVisible(false);}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displaytitle=")){chart.setTitle(strLine.substring(i1+13,strLine.length()));}}
        }
        continue;
             }
             // select colx and coly
     // replace tabs by spaces
      strLine=strLine.replaceAll("[\t\n\u000B\u0009\f]"," ");
                 sx=SF.NthWord(strLine,clx);
                 sy=SF.NthWord(strLine,cly);
                 sxerr=SF.NthWord(strLine,clxerr);
                 syerr=SF.NthWord(strLine,Math.abs(clyerr));
             // System.out.println(sx+" "+sy+" "+sxerr+" "+syerr);

               Double p = new Double(0.0);
   if(sx.length()!=0&&sy.length()!=0&&sxerr.length()!=0&&syerr.length()!=0){
               try{
                    sx=sx.replace('D','E');
                    sy=sy.replace('D','E');
                    sxerr=sxerr.replace('D','E');
                    syerr=syerr.replace('D','E');
                   if(clyerr>=0)
                   { if(clxerr==0){sxerr="0";}
                     if(clyerr==0){syerr="0";}
                     data[0][j]=p.parseDouble(sy);
                     if (detymin&data[0][j]<ymin){ymin=data[0][j];}
                     if (detymax&data[0][j]>ymax){ymax=data[0][j];}
                     data[1][j]=p.parseDouble(sy)+p.parseDouble(syerr);
                     data[2][j]=p.parseDouble(sy)-p.parseDouble(syerr);
                     data[3][j]=p.parseDouble(sx);
                     if (detxmin&data[3][j]<xmin){xmin=data[3][j];}
                     if (detxmax&data[3][j]>xmax){xmax=data[3][j];}
                     data[4][j]=p.parseDouble(sx)+p.parseDouble(sxerr);;
                     data[5][j]=p.parseDouble(sx)-p.parseDouble(sxerr);;
                   //System.out.println(data[0][j]+" "+data[3][j]+" "+j);
                   }
                   else
                   {bdata[1][j]=p.parseDouble(sx);
                     if (detxmin&bdata[1][j]<xmin){xmin=bdata[1][j];}
                     if (detxmax&bdata[1][j]>xmax){xmax=bdata[1][j];}
                    bdata[0][j]=p.parseDouble(sy);
                     if (detymin&bdata[0][j]<ymin){ymin=bdata[0][j];}
                     if (detymax&bdata[0][j]>ymax){ymax=bdata[0][j];}
                    bdata[2][j]=p.parseDouble(syerr);
                    if (bdata[2][j]<0){bdata[2][j]=0;}
                    bdata[2][j]=scale*Math.sqrt(bdata[2][j]);
                  }
                    ++j;
                   }
                   catch(NumberFormatException e){--j;//System.exit(1);
                                                  }
                                                          }
               }
               if(j==maxnofpoints){maxnofpoints*=2;j=maxnofpoints;}
                 else {
               if (j>0)
               {// here fill the rest of the array with the same values
                for(int jj=j;jj<maxnofpoints;++jj)
                  {data[0][jj]=data[0][j-1];data[1][jj]=data[1][j-1];data[2][jj]=data[2][j-1];
                    if(clyerr>=0){data[3][jj]=data[3][j-1];data[4][jj]=data[4][j-1];data[5][jj]=data[5][j-1];
                                }
                   
                  }

               if(clyerr>=0)
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


