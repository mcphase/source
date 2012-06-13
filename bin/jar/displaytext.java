/*
 * Copyright (c) 1995 - 2008 Martin Rotter All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *   - Neither the name of Sun Microsystems nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */ 

/* displaytext.java requires no other files. */

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.lang.ref.WeakReference;


public class displaytext extends JPanel implements FileListener {
    protected JTextArea textArea;
    protected JScrollPane scrollPane;
    private final static String newline = "\n";
    private static String filename;
    static myStringfunc SF=new myStringfunc();
    static int w,h;

    public displaytext(String filen) {
        super(new GridBagLayout());

       // Create the monitor
       FileMonitor monitor = new FileMonitor (1000);

      // Add some files to listen for
        monitor.addFile (new File (filen));

      // Add a listener
      monitor.addListener (this);

        textArea = new JTextArea(h, w);
        textArea.setEditable(false);
        textArea.setFont(new Font("Courier New",0,12));
        scrollPane = new JScrollPane(textArea);

      //Add Components to this panel.
        GridBagConstraints c = new GridBagConstraints();
        c.gridwidth = GridBagConstraints.REMAINDER;

        c.fill = GridBagConstraints.HORIZONTAL;
//        add(textField, c);

        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        add(scrollPane, c);

      try {FileInputStream fstream = new FileInputStream(filen);
           DataInputStream in = new DataInputStream(fstream);
           while (in.available() !=0)
           {
//            System.out.println (in.readLine());
            textArea.append(in.readLine()); 
            textArea.append(newline); 
           }
         in.close();fstream.close();
          }
         catch (Exception e) { System.err.println("File input error"); }


    }

    public displaytext(String filen,int rows,int columns) {
        super(new GridBagLayout());

       // Create the monitor
       FileMonitor monitor = new FileMonitor (1000);

      // Add some files to listen for
        monitor.addFile (new File (filen));

      // Add a listener
      monitor.addListener (this);

        textArea = new JTextArea(rows, columns);
        textArea.setEditable(false);
        textArea.setFont(new Font("Courier New",0,12));
        scrollPane = new JScrollPane(textArea);

      //Add Components to this panel.
        GridBagConstraints c = new GridBagConstraints();
        c.gridwidth = GridBagConstraints.REMAINDER;

        c.fill = GridBagConstraints.HORIZONTAL;
//        add(textField, c);

        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        add(scrollPane, c);

      try {FileInputStream fstream = new FileInputStream(filen);
           DataInputStream in = new DataInputStream(fstream);
           while (in.available() !=0)
           {
//            System.out.println (in.readLine());
            textArea.append(in.readLine()); 
            textArea.append(newline); 
           }
         in.close();fstream.close();
          }
         catch (Exception e) { System.err.println("File input error"); }

    }

    public void fileChanged (File file)
    { //System.out.println ("File changed: " + file);
     JViewport viewport= new JViewport();
     viewport=scrollPane.getViewport();
  
    try {FileInputStream fstream = new FileInputStream(file);
           DataInputStream in = new DataInputStream(fstream);
           textArea.setText("");
   
           while (in.available() !=0)
           {
//            System.out.println (in.readLine());
            textArea.append(in.readLine()); 
            textArea.append(newline); 
            }
         in.close();fstream.close();
         scrollPane.setViewport(viewport); 
         }
         catch (Exception e) { System.err.println("File input error"); System.exit(1);}
    }

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event dispatch thread.
     */
    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame(filename);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Add contents to the window.
        frame.getContentPane().add(new displaytext(filename));
        
	//frame.setSize(new Dimension(w, h));
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
//        frame.setLocation(screenSize.width-w*10, screenSize.height-h*10);
        frame.setLocation(screenSize.width-w*8, 0);
    
      //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
 


    public static void main(String args[]) {
      String s;String ss;
      if (args.length<1)
       {System.out.println("- too few arguments...\n");
        System.out.println("  program displaytext - show and watch text file by viewing a text box on screen\n\n");
        System.out.println("use as:  display [options] filename\n\n");
        System.out.println("	 filename ..... filename of textfile\n");
        System.out.println("	 option -w 100 ..... width of display\n");
        System.out.println("	 option -h 10  ..... height of display\n\n");
        System.exit(0);
       }
      int j=0;int k=0;
     w = 20; h =5; // default width and hight
     Double p = new Double(0.0);
      s=args[0];s=SF.TrimString(s); // command line arguments are treated here
       //look if options are present
       while(SF.TrimString(s).substring(0, 1).equalsIgnoreCase("-"))
          {// yes there are options
           if(SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-w")) // option "-w width"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
              ss=SF.FirstWord(s);
              w=p.valueOf(SF.DataCol(ss)).intValue();
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
           else if(SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-h")) // option "-h height"
            {s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
              ss=SF.FirstWord(s);
              h=p.valueOf(SF.DataCol(ss)).intValue();
             s=SF.DropWord(s); if (s.length()==0){++k;s=args[k];s=SF.TrimString(s);}
            }
          }
       
        filename=s;
        //Schedule a job for the event dispatch thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }
}



