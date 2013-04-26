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
import java.util.regex.*;
import javax.swing.text.BadLocationException;
import javax.swing.text.Utilities; 


public class displaytext extends JPanel implements FileListener {
    protected JTextArea textArea, numberArea;
    protected JScrollPane scrollPane;
	protected JButton saveButton;
	
    private boolean areaLineBreak = false;
    private boolean numberLineBreak = true; 
	
    private final static String newline = "\n";
    private static String filename;
    static myStringfunc SF=new myStringfunc();
    static int w,h;
	
	static boolean relative = false;
	
	protected int xval, yval, /* Test */xmax, ymax;

    public displaytext(String filen) {
        super(new GridBagLayout());

       // Create the monitor
       FileMonitor monitor = new FileMonitor (1000);

      // Add some files to listen for
        monitor.addFile (new File (filen));

      // Add a listener
      monitor.addListener (this);

        textArea = new JTextArea(h, w);
        textArea.setEditable(true);
        textArea.setFont(new Font("Courier New",0,12));
        scrollPane = new JScrollPane(textArea);
		
        numberArea = new JTextArea();
		numberArea.setMargin(new Insets(0, 4, 0, 4));
        numberArea.setBackground(new Color(240, 240, 255));
        numberArea.setForeground(new Color(180, 180, 180));
        numberArea.setEditable(false);
		scrollPane.setRowHeaderView(numberArea);
		
		// Load scrollBar values
		xval = scrollPane.getHorizontalScrollBar().getValue();
		yval = scrollPane.getVerticalScrollBar().getValue();
		/* Test */
		xmax = scrollPane.getHorizontalScrollBar().getMaximum();
		ymax = scrollPane.getVerticalScrollBar().getMaximum();
		
		scrollPane.getVerticalScrollBar().addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseReleased(java.awt.event.MouseEvent e) {
		// Load scrollBar values
		xval = scrollPane.getHorizontalScrollBar().getValue();
		yval = scrollPane.getVerticalScrollBar().getValue();
		/* Test */
		xmax = scrollPane.getHorizontalScrollBar().getMaximum();
		ymax = scrollPane.getVerticalScrollBar().getMaximum();
			}
		});
		
		saveButton = new JButton("save");
		saveButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				try {
					byte[] toSave = textArea.getText().getBytes();
					FileOutputStream fos = new FileOutputStream(filename);
					fos.write(toSave);
					fos.close();
				} catch (Exception ex) {
					System.out.println("FILE SAVING ERROR! " + ex);
				}
			}
		});

      //Add Components to this panel.
        GridBagConstraints c = new GridBagConstraints();
        c.gridwidth = GridBagConstraints.REMAINDER;

        c.fill = GridBagConstraints.HORIZONTAL;
//        add(textField, c);

        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        add(scrollPane, c);
		
		//Add Components to this panel.
        GridBagConstraints g = new GridBagConstraints();
        g.gridwidth = GridBagConstraints.REMAINDER;

        g.fill = GridBagConstraints.HORIZONTAL;
//        add(textField, c);

        g.fill = GridBagConstraints.BOTH;
		add(saveButton, g);

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
		 
		 gibNummern();


    }
	
	 public void gibNummern() {
        StringBuilder sb = new StringBuilder();
        int zeilenZahl = textArea.getLineCount();
        int textLength = textArea.getText().length();

        if (!areaLineBreak & zeilenZahl != 0) {
            for (int i = 0; i < zeilenZahl; i++) {
                sb.append(String.valueOf(i + 1) + "\n");
            }
        }

        try {
            if (numberLineBreak && areaLineBreak) {
                Pattern p = Pattern.compile(System
                        .getProperty("line.separator"));
                sb.append("1\n");
                int n = 1;
                for (int i = 0; i < textLength; i++) {
                    int lineEnd = Utilities.getRowEnd(textArea, i);
                    int lineStart = Utilities.getRowEnd(textArea, i + 1);
                    if (lineEnd < lineStart) {
                        String s = textArea.getText().substring(lineEnd,
                                lineStart);
                        Matcher m = p.matcher(s);
                        boolean result = m.find();
                        if (!result) {
                            sb.append("\n");
                        } else {
                            n++;
                            sb.append(String.valueOf(n) + "\n");
                        }
                    }
                }

            } else if (!numberLineBreak && areaLineBreak) {
                int n = 0;

                for (int i = 0; i < textLength; i++) {
                    int lineEnd = Utilities.getRowEnd(textArea, i);
                    int lineStart = Utilities.getRowEnd(textArea, i + 1);
                    if (lineEnd < lineStart) {
                        n++;
                        sb.append(String.valueOf(n) + "\n");
                    }
                }

                sb.append(String.valueOf(++n) + "\n");
            }
        } catch (BadLocationException e1) {
            e1.printStackTrace();
        }
        numberArea.setText(sb.toString());
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
		
		// Load scrollBar values
		xval = scrollPane.getHorizontalScrollBar().getValue();
		yval = scrollPane.getVerticalScrollBar().getValue();
		/* Test */
		xmax = scrollPane.getHorizontalScrollBar().getMaximum();
		ymax = scrollPane.getVerticalScrollBar().getMaximum();
		
		scrollPane.getVerticalScrollBar().addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseReleased(java.awt.event.MouseEvent e) {
		// Load scrollBar values
		xval = scrollPane.getHorizontalScrollBar().getValue();
		yval = scrollPane.getVerticalScrollBar().getValue();
		/* Test */
		xmax = scrollPane.getHorizontalScrollBar().getMaximum();
		ymax = scrollPane.getVerticalScrollBar().getMaximum();
			}
		});

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
	
	static int Max(int i1, int i2) {
		return (i1 > i2) ? i1 : i2;
	}
	
	static int Min(int i1, int i2) {
		return (i1 < i2) ? i1 : i2;
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
		 /* %#* Nikolai */
		 JScrollBar hbar = scrollPane.getHorizontalScrollBar();
		 JScrollBar vbar = scrollPane.getVerticalScrollBar();
		 if (relative) {
			double xrel, yrel;
			xrel = (double)xval / xmax;
			yrel = (double)yval / ymax;
			hbar.setValues((int)(xrel * hbar.getMaximum()), 1, 0, hbar.getMaximum());
			vbar.setValues((int)(yrel * vbar.getMaximum()), 1, 0, vbar.getMaximum());
		 } else {
			hbar.setValues(Min(xval, hbar.getMaximum()), 1, 0, hbar.getMaximum());
			vbar.setValues(Min(yval, vbar.getMaximum()), 1, 0, vbar.getMaximum());
		 }
		 xmax = hbar.getMaximum();
		 ymax = vbar.getMaximum();
		 scrollPane.setHorizontalScrollBar(hbar);
		 scrollPane.setVerticalScrollBar(vbar);
		 /* ----------- */
         }
         catch (Exception e) { System.err.println("File input error"); System.exit(1);}
		 gibNummern();
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
        System.out.println("  program displaytext - show and watch text file by viewing a text box on screen\n");
		System.out.println("                        NOTE: file will be automatically reloaded when it is changed\n\n");
        System.out.println("use as:  displaytext [options] filename\n\n");
        System.out.println("	 filename ..... filename of textfile\n");
        System.out.println("	 option -w 100 ..... width of display\n");
        System.out.println("	 option -h 10  ..... height of display\n\n");
        System.out.println("	 option -r ..... reset scrollbar relative to old position when reloading\n\n");
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
		   else if (SF.TrimString(s).substring(0, 2).equalsIgnoreCase("-r")) {
				relative = true;
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



