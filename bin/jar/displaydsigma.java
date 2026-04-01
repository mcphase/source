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

public class displaydsigma extends Panel implements Runnable {
 XYSeriesCollection collection = new XYSeriesCollection();
 JFreeChart chart;
 ChartPanel chartPanel;
 XYLineAndShapeRenderer renderer;
 Button bRot=new Button("save spectrum.jpg");
 Thread myThread = null;
 static String[] file;
 static int[] colx;
 static int[] coly;
 static String [] legend;
 static String xText = "";
 static String yText = "";
 static FileInputStream ff;

 public void start(){ myThread = new Thread (this); myThread.start();}

 public void stop(){myThread = null;}

 public void run(){ while(myThread!=null){
    try{Thread.sleep(500);
       }catch(Exception ignored){}

 File fileIni;
 boolean legendVisible = false;
 boolean linesVisible = true;
 String titleText = null;
 try{
 String s="abcdefghd";
 final XYSeries[] newSeries = new XYSeries[file.length];
 for (int i=0;i<file.length;++i)
 {newSeries[i] = new XYSeries(s.substring(i,i+1));

 fileIni = new File(file[i]);
 ff = new FileInputStream(fileIni);
    DataInputStream inStream = new DataInputStream(ff);
    String strLine;
    String tit;
    String sx="1";
    String sy="1";
    int clx = colx[i];
    int cly = coly[i];

       tit =  inStream.readLine();
    int nofpoints=0;
    while (inStream.available() > 0)
    { strLine = inStream.readLine();

      if ((strLine.length() == 0)
        ||(strLine.substring(0, 1).equalsIgnoreCase("#")))
      {
      for(int i1=0;i1<=strLine.length();++i1)
       {if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";legendVisible=true;}}
        if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";legendVisible=false;}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext=")){xText=strLine.substring(i1+13,strLine.length());}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext=")){yText=strLine.substring(i1+13,strLine.length());}}
        if(i1<=strLine.length()-17){if(strLine.substring(i1,i1+17).equalsIgnoreCase("displaylines=true")){linesVisible=true;}}
        if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylines=false")){linesVisible=false;}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displaytitle=")){titleText=strLine.substring(i1+13,strLine.length());}}
        }

        continue;
      }

      // select colx and coly
      sx=TrimString(strLine);
      sy=TrimString(strLine);
      int cx =clx-1;
      int cy =cly-1;

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
       cx=sx.indexOf(" ");
       cy=sy.indexOf(" ");
       if (cx>0) {sx=sx.substring(0,cx);}
       if (cy>0) {sy=sy.substring(0,cy);}

      Double p = new Double(0.0);
 try{
      newSeries[i].add(p.parseDouble(sx),p.parseDouble(sy));
      ++nofpoints;
      }
      catch(NumberFormatException e){;}
    }
    if (nofpoints<2){
      Double p = new Double(0.0);
 try{
      newSeries[i].add(p.parseDouble(sx)*1.1,p.parseDouble(sy)*1.1);
      ++nofpoints;
      }
      catch(NumberFormatException e){;}
    }

   ff.close();
    }

   // swap data on the EDT to avoid concurrent modification
   final String finalXText = xText;
   final String finalYText = yText;
   final String finalTitle = titleText;
   final boolean finalLegend = legendVisible;
   final boolean finalLines = linesVisible;
   final XYSeriesCollection newCollection = new XYSeriesCollection();
   for (int i=0;i<newSeries.length;++i) newCollection.addSeries(newSeries[i]);
   SwingUtilities.invokeLater(new Runnable() { public void run() {
     chart.getXYPlot().setDataset(newCollection);
     chart.getXYPlot().getDomainAxis().setLabel(finalXText);
     chart.getXYPlot().getRangeAxis().setLabel(finalYText);
     if(finalTitle!=null) chart.setTitle(finalTitle);
     chart.getLegend().setVisible(finalLegend);
     for(int j=0;j<newSeries.length;++j) renderer.setSeriesLinesVisible(j,finalLines);
   }});

 }
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
    }

    catch (FileNotFoundException e)
    {
      System.out.println("File not found: " + e.getLocalizedMessage());
    }

    catch (IOException e)
    {
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
    }

 }}

 public void update(Graphics g){paint(g);}

 protected void initChart(){
    chart = ChartFactory.createXYLineChart(
        "Scattering Cross Section", xText, yText, collection,
        PlotOrientation.VERTICAL, true, true, false);
    renderer = new XYLineAndShapeRenderer(true, true);
    chart.getXYPlot().setRenderer(renderer);
    chart.getLegend().setVisible(false);

    String s="abcdefghd";
    for (int i=0;i<file.length;++i)
    {collection.addSeries(new XYSeries(s.substring(i,i+1)));
     renderer.setSeriesLinesVisible(i, true);
     renderer.setSeriesShapesVisible(i, true);
     // diamond marker, size varies per series
     final double sz=1+3*i;
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
 String title="Scattering Cross Section (barn/meV/sr/f.u.)";
 for(int i=0; i<args.length-1;	i+=3)
 {file[j]=args[i+2];
  Integer pp;
  ss=args[i];
  colx[j]=p.valueOf(ss).intValue();
  ss=args[i+1];
  coly[j]=p.valueOf(ss).intValue();
  ++j;
}
 legend = new String[file.length];

 Frame myFrame = new Frame(title);
 displaydsigma myPanel = new displaydsigma();
	myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
 myPanel.initChart();
       myFrame.add(myPanel.bRot);
 myFrame.add(myPanel);
       myFrame.pack();
 myFrame.setSize(400,400);
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
