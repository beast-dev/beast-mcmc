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

			// Setup Generate XML
			action = ((FileMenuHandler) frame).getGenerateXMLAction();
			if (action != null) {
				item = new JMenuItem(action);
				item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_G,
						MenuBarFactory.MENU_MASK));
				menu.add(item);
			}

			// Setup Import
			action = frame.getImportAction();
			if (action != null) {
				item = new JMenuItem(action);
				item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I,
						MenuBarFactory.MENU_MASK));
				menu.add(item);
			}

		}// END: instanceof check

		// TODO: setup Recent files
		menu.addSeparator();

		if (application.getRecentFileMenu() != null) {
			JMenu subMenu = application.getRecentFileMenu();
			menu.add(subMenu);

			menu.addSeparator();
		}
		
//		menu.addSeparator();

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
