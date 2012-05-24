package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

@SuppressWarnings("serial")
public class EditMenu extends JMenu {

	private final static String name = "Edit";
	private JMenuItem preferencesItem;

	public EditMenu() {

		super(name);
		
		this.setName(name);
		this.setMnemonic(KeyEvent.VK_E);
		preferencesItem = new JMenuItem("Preferences");
		this.add(preferencesItem);
		preferencesItem.addActionListener(new ListenMenuItemPreferences());
		
	}// END: constructor

	private class ListenMenuItemPreferences implements ActionListener {
		public void actionPerformed(ActionEvent ev) {
			System.out.println("TODO");
		}
	}//END: ListenMenuItemPreferences
	
}// END: class
