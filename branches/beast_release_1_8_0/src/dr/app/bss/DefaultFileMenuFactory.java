package dr.app.bss;

import jam.framework.AbstractFrame;
import jam.framework.Application;
import jam.framework.MenuBarFactory;
import jam.framework.MenuFactory;

import java.awt.event.KeyEvent;

import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class DefaultFileMenuFactory implements MenuFactory {

	public DefaultFileMenuFactory() {
	}// END: Constructor

	@Override
	public void populateMenu(JMenu menu, AbstractFrame frame) {

		JMenuItem item;
		Action action;

		Application application = Application.getApplication();
		menu.setMnemonic('F');

		if (frame instanceof FileMenuHandler) {

			// Setup Open
			action = ((FileMenuHandler) frame).getLoadSettingsAction();
			if (action != null) {
				item = new JMenuItem(action);
				item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_L,
						MenuBarFactory.MENU_MASK));
				menu.add(item);
			}
			
			// Setup Save As
			action = ((FileMenuHandler) frame).getSaveSettingsAction();
			if (action != null) {
				item = new JMenuItem(action);
				item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,
						MenuBarFactory.MENU_MASK));
				menu.add(item);
			}

		}// END: instanceof check

		menu.addSeparator();

		// Setup Exit
		item = new JMenuItem(application.getExitAction());
		menu.add(item);

	}// END: populateMenu

	@Override
	public String getMenuName() {
		return "File";
	}// END: getMenuName

	@Override
	public int getPreferredAlignment() {
		return LEFT;
	}// END: getPreferredAlignment

}// END: class
