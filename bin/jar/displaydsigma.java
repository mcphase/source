import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javachart.chart.*;
import java.io.*;
import java.lang.*;
import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGImageEncoder;

public class displaydsigma extends Panel implements Runnable {
 public LineChart chart = new LineChart("sample");
 Button bRot=new Button("save spectrum.jpg");                       //erstellt einen Button
 Thread myThread = null;
 static   double vals[]={0.,1.};
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
       // here do something

 File fileIni;
 try{
 for (int i=0;i<file.length;++i)
 {String s="abcdefghd";
      Dataset ds = chart.getDataset(s.substring(i,i+1));
      ds.getData().removeAllElements();

 fileIni = new File(file[i]);
 ff = new FileInputStream(fileIni);
    //ffnen der Datei
    DataInputStream inStream = new DataInputStream(ff);
    String strLine;
    String tit;
    String sx="1";
    String sy="1";
    int clx = colx[i];
    int cly = coly[i];   

//    displaydsigma app = new displaydsigma();
//    app.setSize(640, 640);

//    Dataset ds = chart.getDataset("xy");
    ds.getGc().setMarkerStyle(Gc.MK_DIAMOND);
    ds.getGc().setMarkerSize(1+3*i);
//    ds.getGc().setLineWidth(1);
//    ds.getGc().setFillColor(Color.blue);
       tit =  inStream.readLine();
     //Auslesen der Datei
    int nofpoints=0;
    while (inStream.available() > 0)
    { strLine = inStream.readLine();

      if ((strLine.length() == 0)
        ||(strLine.substring(0, 1).equalsIgnoreCase("#")))
      {
      for(int i1=0;i1<=strLine.length();++i1)
       {if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";chart.setLegendVisible(true);}}
        if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";chart.setLegendVisible(false);}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext=")){xText=strLine.substring(i1+13,strLine.length());}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext=")){yText=strLine.substring(i1+13,strLine.length());}}
        if(i1<=strLine.length()-17){if(strLine.substring(i1,i1+17).equalsIgnoreCase("displaylines=true")){chart.setLineVisible(true);}}
        if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylines=false")){chart.setLineVisible(false);}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displaytitle=")){chart.getBackground().setTitleString(strLine.substring(i1+13,strLine.length()));}}
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
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {p.parseDouble(sx)};
//    double[] myDatay = {p.parseDouble(sy)};
 try{
      Datum d = new Datum(p.parseDouble(sx),p.parseDouble(sy),null);
      ds.addDatum(d);++nofpoints;
      }
      catch(NumberFormatException e){;}
    }  
    if (nofpoints<2){
      Double p = new Double(0.0);
 try{
      Datum d = new Datum(p.parseDouble(sx)*1.1,p.parseDouble(sy)*1.1,null);
      ds.addDatum(d);++nofpoints;
      }
      catch(NumberFormatException e){;}


    }
     
   ff.close();
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
    }
// double[] myDatax = {1, 3, 2, 3,33};
// double[] myDatay = {123, 432, 223, 345,33};

   chart.getXAxis().setTitleString(xText); 
   chart.getYAxis().setTitleString(yText); 
 
// app.setVisible(true);
//  repaint();

 }
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
      //EntSession.CWatch("Unplanmges 'End Of File' in DSN-Konfigurationsdatei!");
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

 

       //  
        repaint();
 }}

 protected void initChart(){ 
    chart.setLineVisible(true);
    chart.setLegendVisible(false);
chart.getXAxis().setTitleString("hallo"); 
chart.getYAxis().setTitleString(yText); 
    String s="abcdefghijkl";
 for (int i=0;i<file.length;++i)
   {//char ii=i;   
    chart.addDataset(s.substring(i,i+1),vals,vals);
   }  
    bRot.addActionListener(new ActionListener(){
    public void actionPerformed(ActionEvent ed){
    try{
         FileOutputStream fos=new FileOutputStream("spectrum.jpg");
         BufferedImage image= new BufferedImage(chart.getWidth(),chart.getHeight(), BufferedImage.TYPE_INT_RGB); 
         Graphics g=image.getGraphics();
         paint(g);
         JPEGImageEncoder encoder= JPEGCodec.createJPEGEncoder(fos); 
         encoder.encode(image);
         fos.close();
    }    catch (FileNotFoundException e)
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

      }
                                                 });
 }


 
 public void update(Graphics g){paint(g);}

 public void paint(Graphics g){try{chart.paint(this,g);}catch(ArrayIndexOutOfBoundsException e){;}}

 public static void main(String[] args){ 
  String ss;
  file = new String[args.length/3];
  colx = new int[args.length/3];
  coly = new int[args.length/3];
   Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {};
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

 Frame myFrame = new Frame(title);
 displaydsigma myPanel = new displaydsigma();
	myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
 myPanel.initChart();
       myFrame.add(myPanel.bRot);                                          //fuegt dem JFrame den Button hinzu
       myFrame.pack();
 myFrame.add(myPanel);
 myFrame.setSize(400,400);
 myFrame.setVisible(true);
 myPanel.start();	
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
