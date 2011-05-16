
/* ---------------------
 * BubbleChartDemo1.java
 * ---------------------
 * (C) Copyright 2003-2008, by Object Refinery Limited.
 */

//package demo;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;

import java.awt.Color;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
//import javax.swing.JButton;
//import javax.swing.SwingConstants;
import java.io.*;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGImageEncoder;

/**
 * A bubble chart demo.
 */
public class displaybubbles extends ApplicationFrame {

static final int MAX_NOF_FILES = 10;
static myStringfunc SF=new myStringfunc();

 public static void main(String[] args) {
          String ss; String s;
      if (args.length<1)
      {System.out.println("- too few arguments...\n");
       System.out.println("  program displaybubbles - show and watch data file by viewing a xy graphic on screen\n\n");
       System.out.println("use as:  displaybubbles xcol ycol intcol filename [xcol1 ycol1 intcol filename1 ...]\n\n");
       System.out.println("         xcol,ycol,intcol ... column to be taken as x-, y- and bubble-area axis\n");
       System.out.println("	 filename ..... filename of datafile\n\n");
       System.exit(0);
      }
       file = new String[MAX_NOF_FILES];
       lastmod = new long[MAX_NOF_FILES];
       colx = new int[MAX_NOF_FILES];
       coly = new int[MAX_NOF_FILES];
       colint = new int[MAX_NOF_FILES];
       Double p = new Double(0.0);
       //      System.out.println(sx+" "+sy);
       //      p.valueOf(strLine);
       //    double[] myDatax = {};
       int j=0;
       String title="displaybubbles";
       s=args[0];s=SF.TrimString(s); // command line arguments are treated here
       for(int i=0;s.length()>0;	i+=0)
       {Integer pp;
       ss=SF.FirstWord(s);
       colx[j]=p.valueOf(ss).intValue();       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       coly[j]=p.valueOf(ss).intValue();       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       colint[j]=p.valueOf(ss).intValue();       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       file[j]=ss;lastmod[j]=0; title=title+" "+ss;++j;if(j>=MAX_NOF_FILES){System.out.println("ERROR: maximum number of files"+j+" exceeded, recompile with larger MAX_NOF_FILES\n\n");System.exit(0);}
       s=SF.DropWord(s); if (s.length()==0&&i<args.length-1){++i;s=args[i];s=SF.TrimString(s);}
       }noffiles=j;
        displaybubbles demo = new displaybubbles(title);
        demo.pack();
        RefineryUtilities.centerFrameOnScreen(demo);
        demo.setVisible(true);
        final Thread updater = demo.new UpdaterThread();
        updater.setDaemon(true);
        updater.start();
    } //main

 JButton bRot=new JButton("save display.jpg");                       //erstellt einen Button
 Box.Filler bRot1=new Box.Filler (new Dimension(350,10),new Dimension(350,10),new Dimension(370,10));                       //erstellt einen Button
// AbstractButton bRot= new AbstractButton();
 static int noffiles;
 static String[] file;
 static long[] lastmod;
 static int[] colx;
 static int[] coly;
 static int[] colint;
 static String [] legend; 
 static String xText = "";
 static String yText = "";
 static String Title = "";
 static LegendTitle Legendt;
 static DefaultXYZDataset dataset;
 static JFreeChart chart;
 static JPanel chartPanel;
// static JFrame displayFrame;
    /**
     * A demonstration application showing a bubble chart.
     *
     * @param title  the frame title.
     */
 public displaybubbles(String title) {
        super(title);
        //displayFrame= new JFrame();
        chartPanel = createDemoPanel();
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
        add(bRot1);
        add(bRot);

   bRot.addActionListener(new ActionListener(){
    public void actionPerformed(ActionEvent ed){
    try{
         FileOutputStream fos=new FileOutputStream("display.jpg");
         BufferedImage image= chart.createBufferedImage(chartPanel.getWidth(),chartPanel.getHeight(),BufferedImage.TYPE_INT_RGB,null); 
         JPEGImageEncoder encoder= JPEGCodec.createJPEGEncoder(fos); 
         encoder.encode(image);
         fos.close();
    }    catch (FileNotFoundException e)
    {    System.out.println("File not found: " + e.getLocalizedMessage());
         //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }
         //Sonstiger Dateifehler
         catch (IOException e)
    {    System.out.println("Dateifehler: " + e.getLocalizedMessage());
         //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }

      } });
                                               
 }// constructor

    /**
     * Creates a chart.
     *
     * @param dataset  the dataset.
     *
     * @return The chart.
     */
    private static JFreeChart createChart(XYZDataset dataset) {
        chart = ChartFactory.createBubbleChart(
                Title, xText, yText, dataset,
                PlotOrientation.HORIZONTAL, true, true, false);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setForegroundAlpha(1.0f);

        XYBubbleRenderer renderer = ( XYBubbleRenderer)plot.getRenderer();
//    XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

        renderer.setSeriesPaint(0, Color.blue);
        renderer.setSeriesPaint(1, Color.red);
        renderer.setSeriesPaint(2, Color.green);
        renderer.setSeriesPaint(3, Color.black);
        renderer.setSeriesPaint(4, Color.orange);
        renderer.setSeriesPaint(5, Color.pink);
       
   for(int i=6;i<=MAX_NOF_FILES;++i){ renderer.setSeriesPaint(i, new Color(70*i%256,140*i % 256,210*i % 256));}
           //renderer.setPlotShapes(true);
           //renderer.setShapesFilled(true);
          //renderer.setSeriesShapesVisible(0, true);
          //renderer.setSeriesShapesVisible(1, true);
          //renderer.setSeriesShapesVisible(2, true);
            renderer.setSeriesShape(0, new Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0));
            renderer.setSeriesShape(1, new Rectangle2D.Double(-3.0, -3.0, 6.0, 6.0));
           // renderer.setSeriesShape(1, ShapeUtilities.createDiamond(4.0f));

        // increase the margins to account for the fact that the auto-range
        // doesn't take into account the bubble size...
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setLowerMargin(0.15);
        domainAxis.setUpperMargin(0.15);
        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setLowerMargin(0.15);
        rangeAxis.setUpperMargin(0.15);
        return chart;
    }

    /**
     * Creates a sample dataset.
     *
     * @return A sample dataset.
     */
    public static XYZDataset createDataset() {
        
         dataset = new DefaultXYZDataset();
        //double[] x = {2.1, 2.3, 2.3, 2.2, 2.2, 1.8, 1.8, 1.9, 2.3, 3.8};
        //double[] y = {14.1, 11.1, 10.0, 8.8, 8.7, 8.4, 5.4, 4.1, 4.1, 25};
        //double[] z = {2.4, 2.7, 2.7, 2.2, 2.2, 2.2, 2.1, 2.2, 1.6, 4};
        //double[][] series = new double[][] { x, y, z };
        //dataset.addSeries("Series 1", series);
        return dataset;
    }

    /**
     * Creates a panel for the demo (used by SuperDemo.java).
     *
     * @return A panel.
     */
    public static JPanel createDemoPanel() {
        JFreeChart chart = createChart(createDataset());
        ChartPanel chartPanel = new ChartPanel(chart);
       	
        chartPanel.setDomainZoomable(true);
        chartPanel.setRangeZoomable(true);
        return chartPanel;
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
                }
                catch (InterruptedException e) {
                    // suppress
                }

      File fileIni;
      int filechanged=0;
      for (int i=0;i<noffiles;++i)
           {fileIni = new File(file[i]);
            if(fileIni.lastModified()!=lastmod[i]){lastmod[i]=fileIni.lastModified();filechanged=1;}
           }
       if(filechanged==1)
      { try{
           for (int i=0;i<noffiles;++i)
           {           
            String s="";
            //XYDataset ds = chart.getXYPlot().getDataset(i);
            //ds.getData().removeAllElements();
            int maxnofpoints=1000;int j=maxnofpoints;
           while(j==maxnofpoints)           
           {double [][] data=new double [3][maxnofpoints];//={{0,1},{0,1},{0,1}};             

            fileIni = new File(file[i]);
            //?ffnen der Datei
             DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));
             String strLine;
             String sx;
             String sy;
             String sint;
             int clx = colx[i];
             int cly = coly[i];   
             int clint = colint[i];

             j=0;
             //Auslesen der Datei
            while (inStream.available() > 0&&j<maxnofpoints)
            {
             strLine = inStream.readLine();
             if ((strLine.length() == 0)
             ||(SF.TrimString(strLine).substring(0, 1).equalsIgnoreCase("#")))
             {
      for(int i1=0;i1<=strLine.length();++i1)
       {//if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";chart.addLegend(chart.getXYPlot().Legendt);}}
        //if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";Legendt=chart.getLegend();chart.removeLegend();}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext=")){chart.getXYPlot().getRangeAxis().setLabel(strLine.substring(i1+13,strLine.length()));}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext=")){chart.getXYPlot().getDomainAxis().setLabel(strLine.substring(i1+13,strLine.length()));}}
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
                 sint=SF.NthWord(strLine,clint);
              //System.out.println(sx+" "+sy+" "+sint);

               Double p = new Double(0.0);
   if(sx.length()!=0&&sy.length()!=0&&sint.length()!=0){
               try{
                    sx=sx.replace('D','E');
                    sy=sy.replace('D','E');
                    sint=sint.replace('D','E');
                   data[1][j]=p.parseDouble(sx);
                   data[0][j]=p.parseDouble(sy);
                   data[2][j]=p.parseDouble(sint);
                   if (data[2][j]<0){data[2][j]=0;}
                   data[2][j]=Math.sqrt(data[2][j]);
                    ++j;
                   }
                   catch(NumberFormatException e){System.exit(1);}
                                                          }
               }   
               if(j==maxnofpoints){maxnofpoints*=2;j=maxnofpoints;}
               else
              {dataset.removeSeries(file[i]);
              dataset.addSeries(file[i],data);}
             }
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
             }
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
    repaint();
  }
 }}}

} // displaybubbles