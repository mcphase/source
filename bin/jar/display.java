import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javachart.chart.*;
import java.io.*;
import java.lang.*;
import java.util.StringTokenizer;
import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGImageEncoder;
public class display extends Panel implements Runnable, MouseListener,KeyListener
{
 static final int MAX_NOF_FILES = 10;

 public static void main(String[] args){
  String ss;String s;
  if (args.length<1)
  {System.out.println("- too few arguments...\n");
    System.out.println("  program display - show and watch data file by viewing a xy graphic on screen\n\n");
    System.out.println("use as:  display xcol ycol filename [xcol1 ycol1 filename1 ...]\n\n");
    System.out.println("         xcol,ycol ... column to be taken as x- and y-axis (toggle axis on/off with '_','|' key)\n");
    System.out.println("	 filename ..... filename of datafile.\n\n");
    System.out.println("	 jpg can be saved (button toggle visible 's')\n");
    System.out.println("	 grid can be toggled visible 'g'\n");
    System.out.println("	 precision of x/y axis can be changed by 'x','X' and 'y','Y'\n");
    System.out.println("	 Data files may contain lines to tune the display output, such as\n");
    System.out.println("	 # displayytext=intensity\n");
    System.out.println("	 # displayxtext=meV \n");
    System.out.println("	 # displaylines=true (toggle also with '-' key))\n");
    System.out.println("	 # displaytitle=My new Graph\n");
    System.out.println("	 # displaylegend=false (toggle also with 'L' key)\n\n");
  System.exit(0);
  }

  file = new String[MAX_NOF_FILES];
  legend = new String[MAX_NOF_FILES];
  lastmod = new long[MAX_NOF_FILES];
  colx = new int[MAX_NOF_FILES];
  coly = new int[MAX_NOF_FILES];
   Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {};
 int j=0;
 String title="display";

 s=args[0];s=SF.TrimString(s);
       for(int i=0;s.length()>0;	i+=0)
       {Integer pp;
       ss=SF.FirstWord(s);
       colx[j]=p.valueOf(ss).intValue();       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       coly[j]=p.valueOf(ss).intValue();       title=title+" "+ss;
       s=SF.DropWord(s); if (s.length()==0){++i;s=args[i];s=SF.TrimString(s);}
       ss=SF.FirstWord(s);
       file[j]=ss;lastmod[j]=0; title=title+" "+ss;++j;if(j>=MAX_NOF_FILES){System.out.println("ERROR: maximum number of files"+j+" exceeded, recompile with larger MAX_NOF_FILES\n\n");System.exit(0);}
       s=SF.DropWord(s); if (s.length()==0&&i<args.length-1){++i;s=args[i];s=SF.TrimString(s);}
       }noffiles=j;

 display myPanel = new display();
         myPanel.initChart();

 Frame myFrame = new Frame(title);
       myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
       myFrame.add(myPanel.bRot);                                          //fügt dem JFrame den Button hinzu
       myFrame.pack();
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
 } //main
   
   public void mouseClicked(MouseEvent e) {
   //  System.out.println("mouse clicked");
    }
    public void mousePressed(MouseEvent e) {
  //   System.out.println("mouse pressed");
     xy[0]=e.getX();
     xy[1]=e.getY();
    }

    public void mouseReleased(MouseEvent e) {
//     System.out.println("mouse released");
     xy[2]=e.getX();
     xy[3]=e.getY();
        double xa= chart.getXAxis().getAxisStart();
        double ya= chart.getYAxis().getAxisStart();
        double xe= chart.getXAxis().getAxisEnd();
        double ye= chart.getYAxis().getAxisEnd();
        double h= chart.getHeight();
        double w= chart.getWidth();
    if(e.getButton()==1)
  { if(xy[1]==xy[3]&&xy[0]==xy[2])
   {//print position of curser in scientific notation
    //System.out.println(xy[0]+" "+xy[1]+" "+xy[2]+" "+xy[3]);
    //System.out.println(h+" "+w);
    //System.out.println(xa+" "+xe+" "+ya+" "+ye);
    double posx=xa+5./3.*(xe-xa)*(xy[0]/w-0.2);
    double posy=ya+5./3.*(ye-ya)*(1-xy[1]/h-0.2);
    System.out.println("x="+posx+" y="+posy);
   }
   else
   {// rescale picture
     // this sets the axis range
    double posx1=xa+5./3.*(xe-xa)*(xy[0]/w-0.2);
    double posy1=ya+5./3.*(ye-ya)*(1-xy[1]/h-0.2);
    double posx2=xa+5./3.*(xe-xa)*(xy[2]/w-0.2);
    double posy2=ya+5./3.*(ye-ya)*(1-xy[3]/h-0.2);
    int reload=0;
if (posx1<posx2){
         //chart.getXAxis().setAutoScale(false);
         chart.getXAxis().setAxisStart(posx1);
         chart.getXAxis().setAxisEnd(posx2);
                                  }
    else {//chart.getXAxis().setAutoScale(true);
         chart.getXAxis().setAxisStart(xxa);
         chart.getXAxis().setAxisEnd(xxe);
         reload=1;
         }
    if (posy1<posy2){
         //chart.getYAxis().setAutoScale(false);
         chart.getYAxis().setAxisStart(posy1);
         chart.getYAxis().setAxisEnd(posy2);
                                  }
    else {//chart.getYAxis().setAutoScale(true);
         chart.getYAxis().setAxisStart(yya);
         chart.getYAxis().setAxisEnd(yye);
         reload=1;
         }
     if(reload==1){reload_data();}
   }
  }
  else // getbutton =2 or 3 ... scroll  picture
  {double posx1=xa+5./3.*(xe-xa)*(xy[0]/w-0.2);
   double posy1=ya+5./3.*(ye-ya)*(1-xy[1]/h-0.2);
   double posx2=xa+5./3.*(xe-xa)*(xy[2]/w-0.2);
   double posy2=ya+5./3.*(ye-ya)*(1-xy[3]/h-0.2);
   double scrollx=posx2-posx1;
   double scrolly=posy2-posy1;
         chart.getXAxis().setAxisStart(xa-scrollx);
         chart.getXAxis().setAxisEnd(xe-scrollx);
         chart.getYAxis().setAxisStart(ya-scrolly);
         chart.getYAxis().setAxisEnd(ye-scrolly);
  } repaint();
 }
  public void mouseEntered(MouseEvent e) {//    System.out.println("mouse entered");
                                         }
  public void mouseExited(MouseEvent e) {   //   System.out.println("mouse exited ");
                                        }
  public void keyPressed(KeyEvent e) {}
  public void keyReleased(KeyEvent e) {}
  public void keyTyped(KeyEvent e) {
                                    if (e.getKeyChar()=='-'||e.getKeyChar()=='-'){chart.setLineVisible(!chart.getLineVisible());}
                                    if (e.getKeyChar()=='_'||e.getKeyChar()=='_'){chart.setXAxisVisible(!chart.isXAxisVisible());}
                                    if (e.getKeyChar()=='|'||e.getKeyChar()=='|'){chart.setYAxisVisible(!chart.isYAxisVisible());}
                                    if (e.getKeyChar()=='l'||e.getKeyChar()=='L'){chart.setLegendVisible(!chart.isLegendVisible());}
                                    if (e.getKeyChar()=='s'||e.getKeyChar()=='S'){bRot.setVisible(!bRot.isVisible());}
                                    if (e.getKeyChar()=='g'||e.getKeyChar()=='G'){chart.getXAxis().setGridVis(!chart.getXAxis().getGridVis());
                                                                                  chart.getYAxis().setGridVis(!chart.getYAxis().getGridVis());}
                                    if (e.getKeyChar()=='x'){chart.getXAxis().setLabelPrecision(chart.getXAxis().getLabelPrecision()+1);}
                                    if (e.getKeyChar()=='X'){chart.getXAxis().setLabelPrecision(chart.getXAxis().getLabelPrecision()-1);}
                                    if (e.getKeyChar()=='y'){chart.getYAxis().setLabelPrecision(chart.getYAxis().getLabelPrecision()+1);}
                                    if (e.getKeyChar()=='Y'){chart.getYAxis().setLabelPrecision(chart.getYAxis().getLabelPrecision()-1);}
                                    if (e.getKeyChar()=='t'){chart.getXAxis().setNumMinTicks(chart.getXAxis().getNumMinTicks()+1);
                                                             chart.getYAxis().setNumMinTicks(chart.getYAxis().getNumMinTicks()+1);}
                                    if (e.getKeyChar()=='T'){chart.getXAxis().setNumMinTicks(chart.getXAxis().getNumMinTicks()-1);
                                                             chart.getYAxis().setNumMinTicks(chart.getYAxis().getNumMinTicks()-1);}
                                    //if (e.getKeyChar()=='p'){ chart.getXAxis().setLogScaling(!chart.getXAxis().getLogScaling());}
                                    //if (e.getKeyChar()=='q'){ chart.getYAxis().setLogScaling(!chart.getYAxis().getLogScaling());}
                                    

                                   repaint();//System.out.println("Key pressed ");
                                   }


 public LineChart chart;
 static myStringfunc SF=new myStringfunc();
 Button bRot;
 Thread myThread = null;
 static double vals[]={0.,1.};
 static String[] file;
 static int noffiles;
 static int[] colx;
 static int[] coly;
 static long[] lastmod;
 static String [] legend; 
 static String xText = "";
 static String yText = "";
 static File fileIni;
 static FileInputStream ff;
 static int xy[]={0,0,0,0};
 static double xxa,xxe,yya,yye; // axes autolimits
  
 public display() {
 chart = new LineChart("display");
 bRot=new Button("save display.jpg"); //erstellt einen Button
 addMouseListener(this);
 addKeyListener(this);
 }

 
 public void start(){ myThread = new Thread (this); myThread.start();
                      myThread.setPriority(1); // be nice 1 is minimum priority
                    }

 public void stop(){myThread = null;}

 public void run(){ 
    while(myThread!=null){
    try{Thread.sleep(500);
       }catch(Exception ignored){}
       // here do something


  int filechanged=0;
           for (int i=0;i<noffiles;++i)
                  {fileIni = new File(file[i]);
                   if(fileIni.lastModified()!=lastmod[i]){lastmod[i]=fileIni.lastModified();filechanged=1;}
                  }
       if(filechanged==1)
      { 
  reload_data(); repaint();
}}}

protected void reload_data()
{
try{
    xxa=1e100;xxe=-1e100;yya=1e100;yye=-1e100;
 for (int i=0;i<noffiles;++i)
 {String s="";
  legend[i]="false";
      Dataset ds = chart.getDataset(file[i]+s.valueOf(i));
      ds.getData().removeAllElements();
 int nofelements=0;
 fileIni = new File(file[i]);
 ff = new FileInputStream(fileIni);
    //ffnen der Datei
    DataInputStream inStream = new DataInputStream(ff);
    String strLine;
    String sx;
    String sy;
    int clx = colx[i];
    int cly = coly[i];   
//    display app = new display();
//    app.setSize(640, 640);

//    Dataset ds = chart.getDataset("xy");
    ds.getGc().setMarkerStyle(Gc.MK_DIAMOND);
    ds.getGc().setMarkerSize(7);
    ds.getGc().setLineWidth(1);
//    ds.getGc().setFillColor(Color.blue);

     //Auslesen der Datei
    while (inStream.available() > 0)
    {
      strLine = inStream.readLine();
      if ((strLine.length() == 0)
        ||(SF.TrimString(strLine).substring(0, 1).equalsIgnoreCase("#")))
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
     //      System.out.println(strLine);
     // replace tabs by spaces
      strLine=strLine.replaceAll("[\t\n\u000B\u0009\f]"," ");
     sx=SF.NthWord(strLine,clx);
     sy=SF.NthWord(strLine,cly);

      Double p = new Double(0.0);
//      p.valueOf(strLine);
//    double[] myDatax = {p.parseDouble(sx)};
//    double[] myDatay = {p.parseDouble(sy)};
//   if(colxnotfound==0&&colynotfound==0){}else{System.out.println("d"+sx+"d"+sy+"d");}
//   if(sx.compareTo("nan")!=0&&sy.compareTo("nan")!=0)
   if(sx.length()!=0&&sy.length()!=0){
//System.out.println("d"+sx+"d"+sy+"d");
   try{
       sx=sx.replace('D','E');
       sy=sy.replace('D','E');
      Datum d = new Datum(p.parseDouble(sx),p.parseDouble(sy),null);
      ds.addDatum(d);++nofelements;
      if(xxa>p.parseDouble(sx)){xxa=p.parseDouble(sx);}
      if(yya>p.parseDouble(sy)){yya=p.parseDouble(sy);}
      if(xxe<p.parseDouble(sx)){xxe=p.parseDouble(sx);}
      if(yye<p.parseDouble(sy)){yye=p.parseDouble(sy);}
      }
      catch(NumberFormatException e){;}
    }}//while instream is available
    if(nofelements==0){Datum d = new Datum(0.0,0.0,null); // add (0,0) if no other points are found
                       ds.addDatum(d);}
    ff.close();
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
    }
// double[] myDatax = {1, 3, 2, 3,33};
// double[] myDatay = {123, 432, 223, 345,33};
// app.setVisible(true);

   chart.getXAxis().setTitleString(xText); 
   chart.getYAxis().setTitleString(yText); 
 }
 catch(EOFException e)
    {System.out.println("EOF: " + e.getLocalizedMessage());
      //EntSession.CWatch("Unplanned 'End Of File' in DSN-Konfigurationsdatei!");
    }
    catch (FileNotFoundException e)
    { System.out.println("File not found: " + e.getLocalizedMessage());
      //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }
    //Sonstiger Dateifehler
    catch (IOException e)
    { System.out.println("Dateifehler: " + e.getLocalizedMessage());
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }
       
 }



 protected void initChart(){ 
 // Give it a title
 //     chart.getBackground().setTitleFont(new Font("Serif", Font.PLAIN, 24));
 //     chart.getBackground().setTitleString("Comparing Apples and Oranges");
chart.getXAxis().setTitleString("hallo"); 
chart.getYAxis().setTitleString(yText);
chart.getXAxis().setMinTickVis(true);
chart.getXAxis().setNumMinTicks(5);
chart.getYAxis().setMinTickVis(true);
chart.getYAxis().setNumMinTicks(5);

 String s="";
 for (int i=0;i<noffiles;++i)
   {
   chart.addDataset(file[i]+s.valueOf(i),vals,vals);
   }  

    bRot.addActionListener(new ActionListener(){
    public void actionPerformed(ActionEvent ed){
    try{
         FileOutputStream fos=new FileOutputStream("display.jpg");
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


}// class display