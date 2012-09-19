//-----------------------------------------------------------------
/**
 * Provides tooltips for AWT Components.
 * @author R. Mark Volkmann, Object Computing, Inc.
 */
import java.awt.*;
import java.awt.event.*;
import java.util.Hashtable;



public class ToolTipManager implements MouseListener, Runnable {

   /**
    * The background color of the tooltips.
    */
   private static final Color BACKGROUND_COLOR = new Color(200, 200, 255);

   /**
    * The number of milliseconds that the mouse must remain
    * over a Component before its tooltip is displayed.
    */
   private static final int WAIT_TIME = 1000; // milliseconds

   /**
    * Only one ToolTipManager object per application is needed
    * so a singleton is created.
    */
   private static ToolTipManager singleton = new ToolTipManager();

   /**
    * The Component over which the mouse is currently positioned.
    */
   private Component currentComponent;

   /**
    * Mapping from Components to their tooltip text.
    */
   private Hashtable componentToTipTable = new Hashtable();

   /**
    * The Component used to display the tooltip.
    */
   private Label label = new Label();

   /**
    * A thread used to determine whether the cursor has remained
    * over a Component long enough to display its tooltip.
    */
   private Thread timerThread = new Thread(this);

   /**
    * The Window in which the tooltip is displayed.
    */
   private Window window;

   /**
    * Creates a ToolTipManager.
    * This is only invoked from this class since only a singleton is needed.
    */
   private ToolTipManager() {
       label.setBackground(BACKGROUND_COLOR);
       timerThread.start();
   }

   /**
    * Creates the Window in which tooltips will be displayed.
    */
   private void createWindow() {
       // Find the Frame parent of currentComponent.
       Component top = currentComponent;
       while (true) {
           Container parent = top.getParent();
           if (parent == null) break;
           top = parent;
       }
       

       // Use that Frame for the parent of the tooltip Window.
       window = new Window((Frame) top);
      

        window.add(label, BorderLayout.CENTER);
      
   }

    /**
     * Prevents the tooltip from being displayed if the timer hasn't expired
     * and hides the tooltip window if it has. 
     */
   private void hideTip() {
       // Prevent the tooltip from being displayed
       // if the timer hasn't expired.
       timerThread.interrupt();
       

       // Hide the tooltip Window if it is already displayed.
       if (window != null) {
           window.setVisible(false);
       }
   }


    /**
     * Adds a tooltip to a given Component.
     * @param component the Component
     * @param tip the tooltip String
     */
   public static void setToolTipText(Component component, String tip) {
       singleton.componentToTipTable.put(component, tip);        
       component.addMouseListener(singleton);
   }


    /**
     * Determines which Component the mouse cursor is over and
     * starts a thread to wait the appropriate amount of time
     * before displaying its tooltip.
     */
   public void mouseEntered(MouseEvent e) {
       currentComponent = (Component) e.getSource();
      

       // Notify the timerThread that the cursor is
       // now over a Component that has a tooltip.
       synchronized (this) {
           notify();
       }
   }


     /**
      * Interrupts the thread that is waiting to display the tooltip
      * and hides the tooltip window.
      */
   public void mouseExited(MouseEvent e) {
       if (e.getSource() == currentComponent) {
           hideTip();
       }
   }


     /**
      * Interrupts the thread that is waiting to display the tooltip
      * and hides the tooltip window.


       */
   public void mousePressed(MouseEvent e) {
       if (e.getSource() == currentComponent) {
           hideTip();
       }
   }


    // We don't need to do anything in these methods.
   public void mouseReleased(MouseEvent e) {}
   public void mouseClicked(MouseEvent e) {}


    /**
     * Removes a tooltip from a given Component.
     * @param component the Component
     */
   public static void removeToolTipText(Component component) {
       singleton.componentToTipTable.remove(component);        
       component.removeMouseListener(singleton);
   }


    /**
     * Implementation of the timer thread.
     * This waits the appropriate amount of time and then displays the tooltip.
     */
   public synchronized void run() {
       while (true) {
           try {
               // Wait until notified that the cursor is over a Component
               // that has a tooltip.
               synchronized (this) {
                   wait();
               }
               

               Thread.sleep(WAIT_TIME);
                   

               // Get the tooltip text to be displayed.
               String tip =
                   (String) componentToTipTable.get(currentComponent);
               label.setText(tip);
               

               if (window == null) {
                   createWindow();
               }
               window.pack();
    

                // Place the tooltip directly below its Component.
               Rectangle bounds = currentComponent.getBounds();
               Point location = currentComponent.getLocationOnScreen();
               window.setLocation(location.x, location.y + bounds.height); 
                

                window.setVisible(true);
           } catch (InterruptedException e) {
               // Don't break out of loop!  Just wait for the mouse
               // to move over another component that has a tooltip.
           }
       }
   }
}

