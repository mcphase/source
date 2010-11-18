/*
 * $Id: JGridDemo.java,v 1.8 2001/02/05 23:28:46 dwd Exp $
 *
 * This software is provided by NOAA for full, free and open release.  It is
 * understood by the recipient/user that NOAA assumes no liability for any
 * errors contained in the code.  Although this software is released without
 * conditions or restrictions in its use, it is expected that appropriate
 * credit be given to its author and to the National Oceanic and Atmospheric
 * Administration should the software be included by the recipient as an
 * element in other product development.
 */
//package gov.noaa.pmel.sgt.demo;

//import gov.noaa.pmel.sgt.demo.TestData;

import gov.noaa.pmel.sgt.swing.JPlotLayout;
import gov.noaa.pmel.sgt.swing.JClassTree;
import gov.noaa.pmel.sgt.swing.prop.GridAttributeDialog;
import gov.noaa.pmel.sgt.JPane;
import gov.noaa.pmel.sgt.AbstractPane;
import gov.noaa.pmel.sgt.GridAttribute;
import gov.noaa.pmel.sgt.ContourLevels;
import gov.noaa.pmel.sgt.CartesianRenderer;
import gov.noaa.pmel.sgt.CartesianGraph;
import gov.noaa.pmel.sgt.GridCartesianRenderer;
import gov.noaa.pmel.sgt.IndexedColorMap;
import gov.noaa.pmel.sgt.ColorMap;
import gov.noaa.pmel.sgt.LinearTransform;

import gov.noaa.pmel.sgt.dm.SGTData;
import gov.noaa.pmel.sgt.dm.SGTPoint;
import gov.noaa.pmel.sgt.dm.SGTLine;
import gov.noaa.pmel.sgt.dm.SGTGrid;
import gov.noaa.pmel.sgt.dm.SGTMetaData;
import gov.noaa.pmel.sgt.dm.SimplePoint;
import gov.noaa.pmel.sgt.dm.SimpleLine;
import gov.noaa.pmel.sgt.dm.SimpleGrid;
import gov.noaa.pmel.sgt.dm.Collection;
import gov.noaa.pmel.sgt.SGLabel;

import gov.noaa.pmel.util.GeoDate;
import gov.noaa.pmel.util.TimeRange;
import gov.noaa.pmel.util.Range2D;
import gov.noaa.pmel.util.Dimension2D;
import gov.noaa.pmel.util.Rectangle2D;
import gov.noaa.pmel.util.Point2D;
import gov.noaa.pmel.util.IllegalTimeValue;

import gov.noaa.pmel.util.GeoDateArray;

import java.awt.*;
import java.awt.print.*;
import javax.swing.*;

import java.awt.event.*;
import java.io.*;
import java.lang.*;
/**
 * Example demonstrating how to use <code>JPlotLayout</code>
 * to create a raster-contour plot.
 * 
 * @author Donald Denbo
 * @version $Revision: 1.8 $, $Date: 2001/02/05 23:28:46 $
 * @since 2.0
 */
public class displaycontour extends JApplet {
 static   double vals[]={0.,1.};
 static String[] file;
 static String title;
 static String xaxistext;
 static String yaxistext;
 static String zaxistext;
 
 static String comment;
 static String comment1;
 static int[] colx;
 static int[] coly;
 static int[] colz;


  static JPlotLayout rpl_;
  private GridAttribute gridAttr_;
  JButton edit_;
  JButton space_ = null;
  JButton tree_;
  JButton print_;
  JButton layout_;
  PageFormat pageFormat = PrinterJob.getPrinterJob().defaultPage();

  public void init() {
    double[] axis1, axis2;
    double[] values;
    double[] zero;

    axis1=new double[10];
    axis2=new double[10];
        
    int count;
        for(count=0; count < 10; count++) {
      axis1[count] = count;
      axis2[count] = count;
    }
    values= new double[10*10];
        for(count=0; count < 100; count++) {
      values[count] = 1.0;
    }
    
    /*
     * Create the demo in the JApplet environment.
     */
    getContentPane().setLayout(new BorderLayout(0,0));
    setBackground(java.awt.Color.white);
    setSize(600,600);
    JPanel main = new JPanel();
    rpl_ = makeGraph(values,axis1,axis2);
    JPanel button = makeButtonPanel(false);
    rpl_.setBatch(true);
    main.add(rpl_, BorderLayout.CENTER);
    JPane gridKeyPane = rpl_.getKeyPane();
    gridKeyPane.setSize(new Dimension(600,100));
    main.add(gridKeyPane, BorderLayout.SOUTH);
    getContentPane().add(main, "Center");
    getContentPane().add(button, "South");
    rpl_.setBatch(false);
  }

  JPanel makeButtonPanel(boolean mark) {
    JPanel button = new JPanel();
    button.setLayout(new FlowLayout());
    tree_ = new JButton("Tree View");
    MyAction myAction = new MyAction();
    tree_.addActionListener(myAction);
    button.add(tree_);
    edit_ = new JButton("Edit GridAttribute");
    edit_.addActionListener(myAction);
    button.add(edit_);

    print_ = new JButton("Print...");
    print_.addActionListener(myAction);
    button.add(print_);
    layout_ = new JButton("Page Layout...");
    layout_.addActionListener(myAction);
    button.add(layout_);
    /*
     * Optionally leave the "mark" button out of the button panel
     */
    if(mark) {
      space_ = new JButton("Add Mark");
      space_.addActionListener(myAction);
      button.add(space_);
    }
    return button;
  }
  public static void main(String[] args) {
  if (args.length<2)
  {System.out.println("- too few arguments...\n");
    System.out.println("  program displaycontour - show data file by viewing a contour/color graphic on screen\n\n");
    System.out.println("use as:  displayclontour xcol ycol zcol filename\n\n");
    System.out.println("         xcol,ycol,zcol ... column to be taken as x- and y-axis and color-axis\n");
    System.out.println("	     filename ..... filename of datafile\n\n");
  System.exit(0);
  }

  String ss;
  file = new String[args.length];
  colx = new int[args.length];
  coly = new int[args.length];
  colz = new int[args.length];
   Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {};
 int j=0;
// title= new String;
 title="displaycontour";
 for(int i=0; i<args.length-1;	i+=4)
 {file[j]=args[i+3];
  Integer pp;
  ss=args[i];xaxistext=ss;
  colx[j]=p.valueOf(ss).intValue();
  ss=args[i+1];yaxistext=ss;
  coly[j]=p.valueOf(ss).intValue();
  ss=args[i+2];zaxistext=ss;
  colz[j]=p.valueOf(ss).intValue();
  ++j; 
  title=title+" "+args[i]+" "+args[i+1]+" "+args[i+2]+" "+args[i+3];
 }
 File fileIni;

    double[] axis1, axis2,a1,a2;
    double[] values,vv;
    double[] zero;
    a1=new double[100000];
    a2=new double[100000];
    vv=new double[100000];
    int ii=0;int jj=0;
    
 try{int i=0;
 String s="";
 fileIni = new File(file[i]);

    //ffnen der Datei
    DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));
    String strLine;
    String sx;
    String sy;
    String sz;
    int clx = colx[i];
    int cly = coly[i];   
    int clz = colz[i];   
    int washere=0;
     //Auslesen der Datei
    while (inStream.available() > 0)
    {
      strLine = inStream.readLine();
      if ((strLine.length() == 0)
        ||(TrimString(strLine).substring(0, 1).equalsIgnoreCase("#")))
      {if((TrimString(strLine).substring(0, 1).equalsIgnoreCase("#"))&&washere==1){washere=2;comment1=strLine ;}
       if((TrimString(strLine).substring(0, 1).equalsIgnoreCase("#"))&&washere==0){washere=1;comment=strLine ;}

        continue;
      }
      
      // select colx and coly
      sx=TrimString(strLine);
      sy=TrimString(strLine);
      sz=TrimString(strLine);
      int cx =clx-1;
      int cy =cly-1;      
      int cz =clz-1;

      while (cx>0)
      {--cx;
       int iPos = sx.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sx=sx.substring(iPos);
       sx=TrimString(sx); 
      }

      while (cy>0)
      {--cy;
       int iPos = sy.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sy=sy.substring(iPos);
       sy=TrimString(sy); 
      }
      
      while (cz>0)
      {--cz;
       int iPos = sz.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sz=sz.substring(iPos);
       sz=TrimString(sz); 
      }
            
       cx=sx.indexOf(" ");
       cy=sy.indexOf(" ");
       cz=sz.indexOf(" ");
       if (cx>0) {sx=sx.substring(0,cx);}
       if (cy>0) {sy=sy.substring(0,cy);}
       if (cz>0) {sz=sz.substring(0,cz);}

      Double pp = new Double(0.0);
      a1[ii]= pp.parseDouble(sx);
      a2[ii]= pp.parseDouble(sy);
      vv[ii]= pp.parseDouble(sz);
   //   System.out.println(sx+" x "+sy+" x "+sz);
    
      if(ii>0){if(a2[ii]<a2[ii-1]&&jj==0){jj=ii;}}
      
     ++ii;
    }  

 }
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
      //EntSession.CWatch("Unplaniges 'End Of File' in DSN-Konfigurationsdatei!");
    }

    catch (FileNotFoundException e)
    {
      System.out.println("File not found: " + e.getLocalizedMessage());
      //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }

    //Sonstiger Dateifehler
    catch (IOException e)
    {
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }

    
    axis1=new double[ii/jj];
    axis2=new double[jj];
      System.out.println("mesh "+ii+"points in grid "+(ii/jj)+"x"+jj);
        
    int count;
        for(count=0; count < ii/jj; count++) {
      axis1[count] = a1[jj*count];}
        for(count=0; count < jj; count++) {
      axis2[count] = a2[count];
    }
    
    values= new double[ii];
        for(count=0; count < ii; count++) {
      values[count] = vv[count];
    }



    /*
     * Create the demo as an application
     */
    displaycontour gd = new displaycontour();
    /*
     * Create a new JFrame to contain the demo.
     */
    JFrame frame = new JFrame("displaycontour");
    JPanel main = new JPanel();
    main.setLayout(new BorderLayout());
    frame.setSize(600,600);
    frame.getContentPane().setLayout(new BorderLayout());
    /*
     * Listen for windowClosing events and dispose of JFrame
     */
    frame.addWindowListener(new java.awt.event.WindowAdapter() {
      public void windowClosing(java.awt.event.WindowEvent event) {
	JFrame fr = (JFrame)event.getSource();
	fr.setVisible(false);
	fr.dispose();
	System.exit(0);
      }
      public void windowOpened(java.awt.event.WindowEvent event) {
	rpl_.getKeyPane().draw();
      }
    });
    /*
     * Create button panel with "mark" button
     */
    JPanel button = gd.makeButtonPanel(true);
    /*
     * Create JPlotLayout and turn batching on.  With batching on the
     * plot will not be updated as components are modified or added to
     * the plot tree.
     */
    rpl_ = gd.makeGraph(values,axis1,axis2);
    rpl_.setBatch(true);
    /*
     * Layout the plot, key, and buttons.
     */
    main.add(rpl_, BorderLayout.CENTER);
    JPane gridKeyPane = rpl_.getKeyPane();
    gridKeyPane.setSize(new Dimension(600,100));
    rpl_.setKeyLayerSizeP(new Dimension2D(6.0, 1.0));
    rpl_.setKeyBoundsP(new Rectangle2D.Double(0.0, 1.0, 6.0, 1.0));
    main.add(gridKeyPane, BorderLayout.SOUTH);
    frame.getContentPane().add(main, BorderLayout.CENTER);
    frame.getContentPane().add(button, BorderLayout.SOUTH);
    frame.pack();
    frame.setVisible(true);
    /*
     * Turn batching off. JPlotLayout will redraw if it has been
     * modified since batching was turned on.
     */
    rpl_.setBatch(false);
  }

  void edit_actionPerformed(java.awt.event.ActionEvent e) {
    /*
     * Create a GridAttributeDialog and set the renderer.
     */
    GridAttributeDialog gad = new GridAttributeDialog();
    gad.setJPane(rpl_);
    CartesianRenderer rend = ((CartesianGraph)rpl_.getFirstLayer().getGraph()).getRenderer();
    gad.setGridCartesianRenderer((GridCartesianRenderer)rend);
    //        gad.setGridAttribute(gridAttr_);
    gad.setVisible(true);
  }

    void tree_actionPerformed(java.awt.event.ActionEvent e) {
      /*
       * Create a JClassTree for the JPlotLayout objects
       */
        JClassTree ct = new JClassTree();
        ct.setModal(false);
        ct.setJPane(rpl_);
        ct.show();
    }

  void print_actionPerformed(ActionEvent e) {
    Color saveColor;

    PrinterJob printJob = PrinterJob.getPrinterJob();
    printJob.setPrintable(rpl_, pageFormat);
    printJob.setJobName("Grid Demo");
    if(printJob.printDialog()) {
      try {
        saveColor = rpl_.getBackground();
        if(!saveColor.equals(Color.white)) {
          rpl_.setBackground(Color.white);
        }
        rpl_.setPageAlign(AbstractPane.TOP,
                          AbstractPane.CENTER);
        RepaintManager currentManager = RepaintManager.currentManager(rpl_);
        currentManager.setDoubleBufferingEnabled(false);
        printJob.print();
        currentManager.setDoubleBufferingEnabled(true);
        rpl_.setBackground(saveColor);
      } catch (PrinterException pe) {
        System.out.println("Error printing: " + pe);
      }
    }

  }

  void layout_actionPerformed(ActionEvent e) {
    PrinterJob pj = PrinterJob.getPrinterJob();
    pageFormat = pj.pageDialog(pageFormat);
  }


    
  JPlotLayout makeGraph(double[] values,double[]axis1,double[]axis2) {
    /*
     * This example uses a pre-created "Layout" for raster time
     * series to simplify the construction of a plot. The
     * JPlotLayout can plot a single grid with
     * a ColorKey, time series with a LineKey, point collection with a
     * PointCollectionKey, and general X-Y plots with a
     * LineKey. JPlotLayout supports zooming, object selection, and
     * object editing.
     */                
    SimpleGrid newData;
    SGTMetaData xMeta;
    SGTMetaData yMeta;
    SGTMetaData zMeta;
    SGLabel keyLabel = new SGLabel("Key Label", "", new Point2D.Double(0.0, 0.0));
    keyLabel.setHeightP(0.16);


      keyLabel.setText("displaycontour");
      xMeta = new SGTMetaData("col"+xaxistext, "");
      yMeta = new SGTMetaData("col"+yaxistext, "");
      zMeta = new SGTMetaData("col"+zaxistext, "");
  
  newData = new SimpleGrid(values, axis1, axis2, "Test Series");

    newData.setXMetaData(xMeta);
    newData.setYMetaData(yMeta);
    newData.setZMetaData(zMeta);
    newData.setKeyTitle(keyLabel);

    JPlotLayout rpl;
    ContourLevels clevels;
    /*
     * Create a test grid with sinasoidal-ramp data.
     */
    Range2D xr = new Range2D(190.0f, 250.0f, 1.0f);
    Range2D yr = new Range2D(0.0f, 45.0f, 1.0f);
    /*
     * Create the layout without a Logo image and with the
     * ColorKey on a separate Pane object.
     */   
    rpl = new JPlotLayout(true, false, false, "test layout", null, true);
    rpl.setEditClasses(false);
    /*
     * Create a GridAttribute for CONTOUR style.
     */
    Range2D datar = new Range2D(-20.0f, 45.0f, 5.0f);
    clevels = ContourLevels.getDefault(datar);
    gridAttr_ = new GridAttribute(clevels);
    /*
     * Create a ColorMap and change the style to RASTER_CONTOUR.
     */
    ColorMap cmap = createColorMap(datar);
    gridAttr_.setColorMap(cmap);
    gridAttr_.setStyle(GridAttribute.RASTER_CONTOUR);
    /*
     * Add the grid to the layout and give a label for
     * the ColorKey.
     */
    rpl.addData(newData, gridAttr_, "col"+zaxistext);
    /*
     * Change the layout's three title lines.
     */        
    rpl.setTitles(title,
                  comment,
                  comment1);
    /*
     * Resize the graph  and place in the "Center" of the frame.
     */
    rpl.setSize(new Dimension(600,500));
    
    /*
     * Resize the key Pane, both the device size and the physical
     * size. Set the size of the key in physical units and place
     * the key pane at the "South" of the frame.
     */
    rpl.setLayerSizeP(new Dimension2D(6.0, 6.0));
    rpl.setKeyLayerSizeP(new Dimension2D(6.0, 1.02));
    rpl.setKeyBoundsP(new Rectangle2D.Double(0.01, 1.01, 5.98, 1.0));

    return rpl;
  }

  ColorMap createColorMap(Range2D datar) {
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

    IndexedColorMap cmap = new IndexedColorMap(red, green, blue);
    cmap.setTransform(new LinearTransform(0.0, (double)red.length,
					  datar.start, datar.end));
    return cmap;
  }
    
  class MyAction implements java.awt.event.ActionListener {
        public void actionPerformed(java.awt.event.ActionEvent event) {
           Object obj = event.getSource();
           
	   if(obj == edit_) 
             edit_actionPerformed(event);
	   
	   if(obj == space_) 
	     System.out.println("  <<Mark>>");
	   
	   if(obj == tree_)
	       tree_actionPerformed(event);

           if(obj == print_) 
               print_actionPerformed(event);
	       
           if(obj == layout_) 
               layout_actionPerformed(event);

        }
    }

 static private String FirstWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(0,iPos);
       fw=TrimString(fw); 
       }
 return(fw); 
 }

 static private String DropWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(iPos);
       fw=TrimString(fw); 
       }
 return(fw); 
 }

 static private String TrimString(String strSource)
 {
    while ((strSource.startsWith(" "))
      && (strSource.length() > 0))
      {
        strSource = strSource.substring(1, strSource.length());
      }

    while ((strSource.endsWith(" "))
        && (strSource.length() > 0))
      {
        strSource = strSource.substring(0, strSource.length() - 1);
      }

    return(strSource);
 }


}


/*
 public static void main(String[] args){ 

 Frame myFrame = new Frame(title);
 display myPanel = new display();
	myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
 myPanel.initChart();
 myFrame.add(myPanel);
 myFrame.setSize(300,300);
 myFrame.setVisible(true);
 myPanel.start();	

//	frame = new JFrame("FileChooserDemo");
//	frame.addWindowListener(new WindowAdapter() {
//	    public void windowClosing(WindowEvent e) {System.exit(0);}
//	});
//	frame.getContentPane().add("Center", panel);
//	frame.pack();
//	frame.setVisible(true);

//	panel.updateState();


 }

 static private String FirstWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(0,iPos);
       fw=TrimString(fw); 
       }
 return(fw); 
 }

 static private String DropWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(iPos);
       fw=TrimString(fw); 
       }
 return(fw); 
 }

 static private String TrimString(String strSource)
 {
    while ((strSource.startsWith(" "))
      && (strSource.length() > 0))
      {
        strSource = strSource.substring(1, strSource.length());
      }

    while ((strSource.endsWith(" "))
        && (strSource.length() > 0))
      {
        strSource = strSource.substring(0, strSource.length() - 1);
      }

    return(strSource);
 }

}*/