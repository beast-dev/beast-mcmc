package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JSeparator;

@SuppressWarnings("serial")
public class HelpMenu extends JMenu {

	private final static String name = "Help";
	private JMenuItem helpItem;

	public HelpMenu() {

		super(name);
		
		this.setName(name);
		this.setMnemonic(KeyEvent.VK_H);
		helpItem = new JMenuItem("Help...");
		helpItem.addActionListener(new ListenMenuHelp());
		
		this.add(new JSeparator());
		this.add(helpItem);

	}// END: constructor

	private class ListenMenuHelp implements ActionListener {
		public void actionPerformed(ActionEvent ev) {
			System.out.println("TODO");
		}
	}//END: ListenMenuHelp
	
}// END: class
