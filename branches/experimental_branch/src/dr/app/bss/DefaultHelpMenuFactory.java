package dr.app.bss;

import jam.framework.AbstractFrame;
import jam.framework.MenuFactory;

import java.awt.Desktop;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class DefaultHelpMenuFactory implements MenuFactory {

	public DefaultHelpMenuFactory() {
	}// END: Constructor

	@Override
	public void populateMenu(JMenu menu, AbstractFrame frame) {

		menu.setMnemonic('H');

		JMenuItem item;

		// Setup About
		item = new JMenuItem();
		item.setText("About " + BeagleSequenceSimulatorApp.LONG_NAME + "...");
		item.addActionListener(new ListenAboutMenuItem());
		menu.add(item);

		// Setup Help
		item = new JMenuItem();
		item.setText("Help...");
		item.addActionListener(new ListenHelpMenuItem());
		menu.add(item);

	}// END: populateMenu

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

	private class ListenHelpMenuItem implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			try {
				Desktop.getDesktop()
						.browse(new URI(
								"http://rega.kuleuven.be/cev/ecv/software/buss"));
			} catch (IOException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this link in your system's standard browser.");
			} catch (URISyntaxException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this link in your system's standard browser.");
			}

		}// END: actionPerformed
	}// END: ListenAboutMenuItem

}// END: class
