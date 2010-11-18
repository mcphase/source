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

    public displaytext(String filen) {
        super(new GridBagLayout());

       // Create the monitor
       FileMonitor monitor = new FileMonitor (1000);

      // Add some files to listen for
        monitor.addFile (new File (filen));

      // Add a listener
      monitor.addListener (this);

        textArea = new JTextArea(5, 20);
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
         catch (Exception e) { System.err.println("File input error"); }
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

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
 


    public static void main(String args[]) {
      if (args.length<1)
       {System.out.println("- too few arguments...\n");
        System.out.println("  program displaytext - show and watch text file by viewing a text box on screen\n\n");
        System.out.println("use as:  display filename\n\n");
        System.out.println("	 filename ..... filename of textfile\n\n");
        System.exit(0);
       }


        filename=args[0];
        //Schedule a job for the event dispatch thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }
}



