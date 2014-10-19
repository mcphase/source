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

/* searchspacestatus.java requires no other files. */

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;


public class searchspacestatus extends JPanel implements ActionListener {
    protected JButton Button;
    protected displaytext text;
    protected JScrollPane scrollPane;
    private final static String newline = "\n";
    private static String filename="results/searchspace.status";

    public searchspacestatus() {
        super(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.gridwidth = GridBagConstraints.REMAINDER;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;

        Button = new JButton ("Stop searchspace");
        Button.addActionListener(this);
        text = new displaytext(filename,5,50);
        scrollPane = new JScrollPane(text);
        add(Button);
        add(scrollPane, c);
    }

    public void actionPerformed(ActionEvent e)
    {try {FileOutputStream fstream = new FileOutputStream(filename);
          PrintStream out= new PrintStream(fstream );   
            out.println ("exiting searchspace");
            out.close();
          }
         catch (Exception f) { System.err.println("Error writing to file"); } 
    }

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event dispatch thread.
     */
    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("searchspace status");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Add contents to the window.
        frame.getContentPane().add(new searchspacestatus());

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
 


    public static void main(String args[]) {
//       Runtime.getRuntime().addShutdownHook(new Thread() {
//                                                           public void run() {
//                                                                              File f = new File(filename);
//                                                                              if (!f.delete())throw new IllegalArgumentException("Delete: deletion failed");
//                                                                             }
//                                                          }
//                                             ); 
        //Schedule a job for the event dispatch thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }
}



