package dr.app.bss;

import jam.framework.AbstractFrame;
import jam.framework.MenuFactory;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

public class DefaultHelpMenuFactory implements MenuFactory {

	public DefaultHelpMenuFactory() {
	}// END: Constructor

	@Override
	public void populateMenu(JMenu menu, AbstractFrame frame) {

		menu.setMnemonic('H');
		
		JMenuItem item;
		
		// Setup About
		item = new JMenuItem();
		item.setText("About " + BeagleSequenceSimulatorApp.LONG_NAME);
		item.addActionListener(new ListenAboutMenuItem());
		menu.add(item);
		
		
		
	}//END: populateMenu
	
	@Override
	public String getMenuName() {
		return "Help";
	}

	@Override
	public int getPreferredAlignment() {
		return LEFT;
	}

	private class ListenAboutMenuItem implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			 AboutDialog ad = new AboutDialog();
             ad.setVisible(true);

		}// END: actionPerformed
	}// END: ListenAboutMenuItem
	
}//END: class
