package dr.app.tracer.application;

import org.virion.jam.framework.MenuFactory;
import org.virion.jam.framework.AbstractFrame;
import org.virion.jam.framework.Application;
import org.virion.jam.framework.MenuBarFactory;

import javax.swing.*;
import java.awt.event.KeyEvent;

import dr.app.beauti.BeautiFrame;

/**
 * @author rambaut
 *         Date: Dec 26, 2004
 *         Time: 11:01:06 AM
 */
public class TracerDefaultFileMenuFactory implements MenuFactory {


    public TracerDefaultFileMenuFactory() {
    }

    public String getMenuName() {
        return "File";
    }

    public void populateMenu(JMenu menu, AbstractFrame frame) {

        JMenuItem item;

        Application application = Application.getApplication();
        menu.setMnemonic('F');

        item = new JMenuItem(application.getNewAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, MenuBarFactory.MENU_MASK));
        menu.add(item);

	    if (frame instanceof TracerFileMenuHandler) {
	        item = new JMenuItem(frame.getImportAction());
	        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I, MenuBarFactory.MENU_MASK));
	        menu.add(item);

	        menu.addSeparator();

	        item = new JMenuItem(((TracerFileMenuHandler)frame).getExportDataAction());
	        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK));
	        menu.add(item);

	        item = new JMenuItem(((TracerFileMenuHandler)frame).getExportPDFAction());
	        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK + KeyEvent.ALT_MASK));
	        menu.add(item);
	    } else {
	        // If the frame is not a TracerFileMenuHandler then create a dummy set of disabled menu options.

	        item = new JMenuItem("Import Trace File...");
	        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I, MenuBarFactory.MENU_MASK));
	        item.setEnabled(false);
	        menu.add(item);

	        menu.addSeparator();

	        item = new JMenuItem("Export Data...");
	        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK));
	        item.setEnabled(false);
	        menu.add(item);

	        item = new JMenuItem("Export PDF...");
	        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK + KeyEvent.ALT_MASK));
	        item.setEnabled(false);
	        menu.add(item);
	    }

        menu.addSeparator();

        item = new JMenuItem(frame.getPrintAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_P, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(application.getPageSetupAction());
        menu.add(item);

        menu.addSeparator();

        item = new JMenuItem(application.getExitAction());
        menu.add(item);
    }

    public int getPreferredAlignment() {
        return LEFT;
    }
}