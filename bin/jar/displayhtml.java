import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;

import javax.swing.text.*;
import javax.swing.text.html.*;
import javax.swing.event.*;

public class displayhtml extends JFrame implements HyperlinkListener {
    protected JScrollPane scrollPane;
	protected JEditorPane htmlPane;
	public int w, h;
	protected ReloadThread reloader;
	
	static String[] params;
	
	public displayhtml() {
		htmlPane = new JEditorPane("text/html", "");
		htmlPane.addHyperlinkListener(this);
		htmlPane.setEditable(false);
		
		scrollPane = new JScrollPane(htmlPane);
		scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		scrollPane.setMinimumSize(new Dimension(10, 10));
		
		reloader = new ReloadThread(this);
		
		this.addWindowListener(new HTMLCloser());
		this.setLayout(new BorderLayout());
		
		this.add(scrollPane, BorderLayout.CENTER);
		
		w = 0;
		h = 0;
	}
	
	public static void main(String[] args) {
		if (args.length<1) {
			System.out.println("- too few arguments...\n");
			System.out.println("  program displayhtml - show html file by viewing a editorpane on screen\n\n");
			System.out.println("use as:  displayhtml <options> filename\n\n");
			System.out.println("         filename ..... filename of html-file (for online-files, please");
			System.out.println("                        don't forget to set \"http://\" before starting program)\n\n");
			System.out.println("options: -w <value> ..... sets the width of displaywindow");
			System.out.println("         -h <value> ..... sets the height of displaywindow");
			System.exit(0);
		}
		
		params = args;
		
		// Set save input...
		if (!(params[params.length - 1].startsWith("http://") || params[params.length - 1].startsWith("file:\\"))) {
			try {
				params[params.length - 1] = new File(params[params.length - 1]).getCanonicalPath();
			} catch (IOException e) {
			}
			params[params.length - 1] = "file:\\" + params[params.length - 1];
		}
		
		displayhtml dh = new displayhtml();
		dh.reloadSite(args);
		
		dh.pack();
		dh.setVisible(true);
		
		dh.reloader.run();
	}
	
	public void reloadSite(String[] args) {
		int k = 0;
		
		// Interpreting argumengts
		while(args[k].substring(0, 1).equals("-")) {
			if (args[k].equalsIgnoreCase("-w")) {
				if (w <= 0) w = Integer.parseInt(args[k+1]);
				++k;
			}
			if (args[k].equalsIgnoreCase("-h")) {
				if (h <= 0) h = Integer.parseInt(args[k+1]);
				++k;
			}
			++k;
		}
		
		// Setting size?
		if (w <= 0) w = this.getWidth();
		if (h <= 0) h = this.getHeight();
		
		// Loading document...
		try {
			Document doc = htmlPane.getDocument();
			doc.putProperty(Document.StreamDescriptionProperty, null);
			htmlPane.setPage(args[args.length - 1]);
		} catch (IOException e) {
			System.out.println("ERROR LOADING HTML-DOCUMENT! " + e);
			System.exit(2);
		} catch (Exception e) {
			System.out.println("ERROR! " + e);
			System.exit(1);
		}
		
		this.setSize(w, h);
	}

	public void hyperlinkUpdate(HyperlinkEvent e) {
        if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
            JEditorPane pane = (JEditorPane) e.getSource();
            if (e instanceof HTMLFrameHyperlinkEvent) {
                HTMLFrameHyperlinkEvent  evt = (HTMLFrameHyperlinkEvent)e;
                HTMLDocument doc = (HTMLDocument)pane.getDocument();
                doc.processHTMLFrameHyperlinkEvent(evt);
            } else {
                try {
                    pane.setPage(e.getURL());
                } catch (Throwable t) {
                    t.printStackTrace();
                }
            }
        }
    }
	
	public class HTMLCloser extends WindowAdapter {
		public void windowClosing(WindowEvent event) {
			JFrame fr = (JFrame)event.getSource();
			fr.setVisible(false);
			fr.dispose();
			System.exit(0);
		}
	}
	
	public class FilesWatcher {
		private final String FilesRegex = "(\\s|\"|\'){1}[\\w\\\\/.]*?(\\s|\"|\'){1}";
		
		protected String FileName;
		protected LinkedList<String> FileNames;
		protected LinkedList<Long> FileMods;
		
		public FilesWatcher(String FileName) {
			this.FileName = FileName;
			FileNames = new LinkedList<String>();
			FileMods = new LinkedList<Long>();
			
			reloadData();
		}
		
		protected String TrimEnd(String In, boolean TrimQuotes) {	
			String Trim = "^\\s+";
			String TrimWithQuote = "^[\\s\"\']+";
			return In.replaceAll((TrimQuotes) ? TrimWithQuote : Trim, "");
		}
		
		protected String TrimStart(String In, boolean TrimQuotes) {	
			String Trim = "\\s+$";
			String TrimWithQuote = "[\\s\'\"]+$";
			return In.replaceAll((TrimQuotes) ? TrimWithQuote : Trim, "");
		}
		
		public void reloadData() {
			FileNames.clear();
			FileMods.clear();
			try {
				Scanner c = new Scanner(new File(FileName));
				File src = new File(FileName);
				String file = "";
				File f;
				Pattern p = Pattern.compile(FilesRegex);
				while(c.hasNext()) {
					file = c.next();
					String tmp = "";
					Matcher m = p.matcher(file);
					while (m.find()) {
						tmp = TrimStart(TrimEnd(m.group(), true), true);
						f = new File(src.getParent(), tmp);
						if (f.isFile()) {
							FileNames.add(tmp);
							FileMods.add(f.lastModified());
						}
					}
				}
				c.close();
			} catch (FileNotFoundException e) {
				System.out.println("FILE ERROR! " + e);
			}
		}
		
		public boolean isModified() {
			Iterator<String> itNames = FileNames.iterator();
			Iterator<Long> itMods = FileMods.iterator();
			
			File f;
			String name = "";
			long mod = 0;
			
			File src = new File(FileName);
			
			while (itNames.hasNext() && itMods.hasNext()) {
				name = itNames.next();
				mod = itMods.next();
				
				f = new File(src.getParent(), name);
				if (f.lastModified() != mod) {
					return true;
				}
			}
			return false;
		}
	}
	
	public class ReloadThread extends Thread {
		protected displayhtml dh;
		protected long oldMod;
		protected FilesWatcher fw;
		
		public ReloadThread(displayhtml html) {
			dh = html;
			oldMod = -1000;
			fw = null;
		}
	
		public void run() {
            setPriority(MIN_PRIORITY);
			
			String fileName = params[params.length - 1];
			if (fileName.startsWith("file:\\")) fileName = fileName.substring(6, fileName.length());
			fw = new FilesWatcher(fileName);
			File f;
			
			while (true) {
				try {
					sleep(500);
			
					dh.w = dh.getWidth();
					dh.h = dh.getHeight();
				
					f = new File(fileName);
					if (oldMod != f.lastModified() || fw.isModified()) {
						dh.reloadSite(params);
						oldMod = f.lastModified();
						fw.reloadData();
						System.out.println("SITE SUCCESSFULLY RELOADED");
					}
				} catch (InterruptedException e) {
					System.out.println("ERROR! Reloading interrupted! " + e);
				}
			}
		}
	}
}