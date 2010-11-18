/* 
 * Viewer.java
 */

import java.awt.*;
import java.awt.event.*;
import java.awt.Toolkit;
import java.util.Timer;
import java.util.TimerTask;

public class Viewer extends Frame {
	private Image image;

	public Viewer(String fileName,int seconds) {
                if(seconds>0) {new ReminderBeep(seconds);}
 		Toolkit toolkit = Toolkit.getDefaultToolkit();
		image = toolkit.getImage(fileName);
		MediaTracker mediaTracker = new MediaTracker(this);
		mediaTracker.addImage(image, 0);
		try
		{
			mediaTracker.waitForID(0);
		}
		catch (InterruptedException ie)
		{
			System.err.println(ie);
			System.exit(1);
		}
		addWindowListener(new WindowAdapter() {
      		public void windowClosing(WindowEvent e) {
        		System.exit(0);
      		}
		});
		setSize(image.getWidth(null), image.getHeight(null));
		setTitle(fileName);
		setVisible(true);
	}

	public void paint(Graphics graphics) {
		graphics.drawImage(image, 0, 0, null);
	}

	public static void main(String[] args) {
         int seconds=0; int i=0;
         // check command string
         if (args[0].equals("-h")){System.out.println("Viewer - program to view images - use as: \njava Viewer [options] imagefile\n\n");
                                   System.out.println("Options: -h               print this help\n");
                                   System.out.println("Options: -timeout 500     end view after 500 seconds\n");
                                   System.exit(0);                     }
         if (args[0].equals("-timeout")){seconds=Integer.valueOf(args[1]).intValue();i=2;}
 		new Viewer(args[i],seconds);
  	}

  public class ReminderBeep {
    Toolkit toolkit;
    Timer timer;

  public ReminderBeep(int seconds) {
    toolkit = Toolkit.getDefaultToolkit();
    timer = new Timer();
    timer.schedule(new RemindTask(), seconds * 1000);
  }

  class RemindTask extends TimerTask {
    public void run() {
//      System.out.println("Time's up!");
//      toolkit.beep();
      //timer.cancel(); //Not necessary because we call System.exit
      System.exit(0); //Stops the AWT thread (and everything else)
    }
  }


 }
 }