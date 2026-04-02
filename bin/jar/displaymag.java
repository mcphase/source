import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.SwingUtilities;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class displaymag extends Panel implements Runnable {
 XYSeriesCollection collection = new XYSeriesCollection();
 JFreeChart chart;
 ChartPanel chartPanel;
 XYLineAndShapeRenderer renderer;
 Thread myThread = null;
 static String[] file;
 static int[] colx;
 static int[] coly;
 static int ctr=0;
 static FileInputStream ff;

 public void start(){ myThread = new Thread (this); myThread.start();}

 public void stop(){myThread = null;}

 public void run(){ while(myThread!=null&&ctr<=10){
    try{Thread.sleep(500);
       }catch(Exception ignored){}
 String sT="";
 File fileIni;
 try{
 String s="abcdefghikl";
 final XYSeries[] newSeries = new XYSeries[file.length];
 for (int i=0;i<file.length;++i)
 {newSeries[i] = new XYSeries(s.substring(i,i+1));

 fileIni = new File(file[i]);
ff = new FileInputStream(fileIni);
    DataInputStream inStream = new DataInputStream(ff);

    String strLine;
    String sx;
    String sy;
    int clx = colx[i];
    int cly = coly[i];
    int clT=3;

     while (inStream.available() > 0)
    { strLine = inStream.readLine();
        ctr=0;
      if ((strLine.length() == 0)
        ||(strLine.substring(0, 1).equalsIgnoreCase("#")))
      {
        continue;
      }

      // select colx and coly
      sx=TrimString(strLine);
      sy=TrimString(strLine);
      sT=TrimString(strLine);
      int cx =clx-1;
      int cy =cly-1;
      int cT =clT-1;

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

      while (cT>0)
      {--cT;
       int iPos = sT.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sT=sT.substring(iPos);
       sT=TrimString(sT);
      }

       cx=sx.indexOf(" ");
       cy=sy.indexOf(" ");
       cT=sT.indexOf(" ");
       if (cx>0) {sx=sx.substring(0,cx);}
       if (cy>0) {sy=sy.substring(0,cy);}
       if (cT>0) {sT=sT.substring(0,cT);}

      Double p = new Double(0.0);
 try{
      newSeries[i].add(p.parseDouble(sx),p.parseDouble(sy));
      }
      catch(NumberFormatException e){;}
    }
    ff.close();
    }

   // swap data on the EDT to avoid concurrent modification
   final String finalT = sT;
   final XYSeriesCollection newCollection = new XYSeriesCollection();
   for (int i=0;i<newSeries.length;++i) newCollection.addSeries(newSeries[i]);
   SwingUtilities.invokeLater(new Runnable() { public void run() {
     chart.getXYPlot().setDataset(newCollection);
     chart.setTitle("T="+finalT);
   }});

 }
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
    }

    catch (FileNotFoundException e)
    {++ctr;
    }

    catch (IOException e)
    {++ctr;
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
    }

 }}

 public void update(Graphics g){paint(g);}

 protected void initChart(){
    chart = ChartFactory.createXYLineChart(
        "Magnetisation", "Magnetic Field (T)", "M[mb/T/ion]", collection,
        PlotOrientation.VERTICAL, true, true, false);
    renderer = new XYLineAndShapeRenderer(false, true);
    chart.getXYPlot().setRenderer(renderer);
    chart.getXYPlot().getDomainAxis().setMinorTickMarksVisible(true);
    chart.getXYPlot().getDomainAxis().setMinorTickCount(5);
    chart.getXYPlot().getRangeAxis().setMinorTickMarksVisible(true);
    chart.getXYPlot().getRangeAxis().setMinorTickCount(5);

    String s="abcdefghikl";
    for (int i=0;i<file.length;++i)
    {collection.addSeries(new XYSeries(s.substring(i,i+1)));
     renderer.setSeriesLinesVisible(i, false);
     renderer.setSeriesShapesVisible(i, true);
     // diamond marker, size 7
     double sz=7;
     renderer.setSeriesShape(i, new java.awt.geom.Path2D.Double(){{moveTo(0,-sz);lineTo(sz,0);lineTo(0,sz);lineTo(-sz,0);closePath();}});
    }

    chartPanel = new ChartPanel(chart);
    setLayout(new BorderLayout());
    add(chartPanel, BorderLayout.CENTER);
 }

 public static void main(String[] args){
  String ss;
  file = new String[args.length/3];
  colx = new int[args.length/3];
  coly = new int[args.length/3];
   Double p = new Double(0.0);
 int j=0;
 String title="Magnetisation";
 for(int i=0; i<args.length-1;	i+=3)
 {file[j]=args[i+2];
  Integer pp;
  ss=args[i];
  colx[j]=p.valueOf(ss).intValue();
  ss=args[i+1];
  coly[j]=p.valueOf(ss).intValue();
  ++j;
}

 Frame myFrame = new Frame(title);
 displaymag myPanel = new displaymag();
	myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
 myPanel.initChart();
 myFrame.add(myPanel);
       myFrame.pack();
 myFrame.setSize(400,400);
 myFrame.setLocation(0,400);
 myFrame.setVisible(true);
 myPanel.start();
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
