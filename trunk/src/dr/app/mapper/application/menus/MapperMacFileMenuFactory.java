package dr.app.mapper.application.menus;

import jam.framework.MenuFactory;
import jam.framework.AbstractFrame;
import jam.framework.Application;
import jam.framework.MenuBarFactory;

import javax.swing.*;
import java.awt.event.KeyEvent;
import java.awt.event.ActionEvent;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MapperMacFileMenuFactory implements MenuFactory {

    public MapperMacFileMenuFactory() {
    }

    public String getMenuName() {
        return "File";
    }

    public void populateMenu(JMenu menu, AbstractFrame frame) {

        Application application = Application.getApplication();
        JMenuItem item;

        item = new JMenuItem(application.getNewAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(application.getOpenAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, MenuBarFactory.MENU_MASK));
        menu.add(item);

        if (application.getRecentFileMenu() != null) {
            JMenu subMenu = application.getRecentFileMenu();
            menu.add(subMenu);
        }

        if (frame instanceof MapperFileMenuHandler) {
            menu.addSeparator();

            item = new JMenuItem(((MapperFileMenuHandler)frame).getImportMeasurementsAction());
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I, MenuBarFactory.MENU_MASK));
            menu.add(item);

            item = new JMenuItem(((MapperFileMenuHandler)frame).getImportLocationsAction());
            menu.add(item);

            item = new JMenuItem(((MapperFileMenuHandler)frame).getImportTreesAction());
            menu.add(item);

            menu.addSeparator();

            item = new JMenuItem(((MapperFileMenuHandler)frame).getExportDataAction());
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK));
            menu.add(item);

            item = new JMenuItem(((MapperFileMenuHandler)frame).getExportPDFAction());
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK + KeyEvent.ALT_MASK));
            menu.add(item);
        } else {
            // If the frame is not a BeautiFrame then create a dummy set of disabled menu options.
            // At present the only situation where this may happen is in Mac OS X when no windows
            // are open and the menubar is created by the hidden frame.

            menu.addSeparator();

            item = new JMenuItem("Import Strains...");
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I, MenuBarFactory.MENU_MASK));
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem("Import Measurements...");
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem("Import Locations...");
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem("Import Trees...");
            item.setEnabled(false);
            menu.add(item);

            menu.addSeparator();

            item = new JMenuItem("Export Data...");
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK));
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem("Export PDF...");
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK));
            item.setEnabled(false);
            menu.add(item);
        }

        menu.addSeparator();

        item = new JMenuItem(frame.getCloseWindowAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_W, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(frame.getSaveAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(frame.getSaveAsAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, MenuBarFactory.MENU_MASK + ActionEvent.SHIFT_MASK));
        menu.add(item);

        menu.addSeparator();

        item = new JMenuItem(frame.getPrintAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_P, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(application.getPageSetupAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_P, MenuBarFactory.MENU_MASK + ActionEvent.SHIFT_MASK));
        menu.add(item);

    }

    public int getPreferredAlignment() {
        return LEFT;
    }
}