package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JSeparator;

@SuppressWarnings("serial")
public class FileMenu extends JMenu {

	private final static String name = "File";
	private JMenuItem exitItem;
	private JMenuItem newItem;
	private JMenuItem importItem;
	private JMenuItem exportItem;

	public FileMenu() {

		super(name);
		
		this.setName(name);
		this.setMnemonic(KeyEvent.VK_F);
		
		newItem = new JMenuItem("New");
		this.add(newItem);
		newItem.addActionListener(new ListenMenuItemNew());
		
		this.add(new JSeparator());
		
		importItem = new JMenuItem("Import");
		this.add(importItem);
		importItem.addActionListener(new ListenMenuItemImport());
		exportItem = new JMenuItem("Export");
		this.add(exportItem);
		exportItem.addActionListener(new ListenMenuItemExport());
		
		this.add(new JSeparator());
		
		exitItem = new JMenuItem("Exit");
		this.add(exitItem);
		exitItem.addActionListener(new ListenMenuItemExit());
		
	}// END: constructor

	private class ListenMenuItemNew implements ActionListener {
		public void actionPerformed(ActionEvent ev) {
			System.out.println("TODO");
		}
	}//END: ListenMenuItemNew
	
	private class ListenMenuItemImport implements ActionListener {
		public void actionPerformed(ActionEvent ev) {
			System.out.println("TODO");
		}
	}//END: ListenMenuItemImport
	
	private class ListenMenuItemExport implements ActionListener {
		public void actionPerformed(ActionEvent ev) {
			System.out.println("TODO");
		}
	}//END: ListenMenuItemExport
	
	private class ListenMenuItemExit implements ActionListener {
		public void actionPerformed(ActionEvent ev) {
			System.exit(0);
		}
	}//END: ListenMenuItemExit
	
}// END: class
