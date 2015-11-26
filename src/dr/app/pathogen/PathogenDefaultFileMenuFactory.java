/*
 * PathogenDefaultFileMenuFactory.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.pathogen;

import jam.framework.AbstractFrame;
import jam.framework.Application;
import jam.framework.MenuBarFactory;
import jam.framework.MenuFactory;

import javax.swing.*;
import java.awt.event.KeyEvent;

/**
 * @author rambaut
 *         Date: Dec 26, 2004
 *         Time: 11:01:06 AM
 */
public class PathogenDefaultFileMenuFactory implements MenuFactory {


    public PathogenDefaultFileMenuFactory() {
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

        item = new JMenuItem(application.getOpenAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(frame.getSaveAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(frame.getSaveAsAction());
        menu.add(item);

        menu.addSeparator();

        // On Windows and Linux platforms, each window has its own menu so items which are not needed
        // are simply missing. In contrast, on Mac, the menu is for the application so items should
        // be enabled/disabled as frames come to the front.
        if (frame instanceof PathogenFrame) {
//            Action action = frame.getImportAction();
//            if (action != null) {
//                item = new JMenuItem(action);
//                item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_I, MenuBarFactory.MENU_MASK));
//                menu.add(item);
//
//                menu.addSeparator();
//            }

            item = new JMenuItem(((PathogenFrame)frame).getExportTreeAction());
            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK));
            menu.add(item);

//            item = new JMenuItem(((TemporalSamplerFrame)frame).getExportGraphicAction());
//            item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, MenuBarFactory.MENU_MASK + KeyEvent.ALT_MASK));
//            menu.add(item);

            item = new JMenuItem(((PathogenFrame)frame).getExportDataAction());
            menu.add(item);

        } else {
            // do nothing
        }

        menu.addSeparator();

        item = new JMenuItem(frame.getPrintAction());
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_P, MenuBarFactory.MENU_MASK));
        menu.add(item);

        item = new JMenuItem(application.getPageSetupAction());
        menu.add(item);

        menu.addSeparator();

        if (application.getRecentFileMenu() != null) {
            JMenu subMenu = application.getRecentFileMenu();
            menu.add(subMenu);

            menu.addSeparator();
        }

        item = new JMenuItem(application.getExitAction());
        menu.add(item);
    }

    public int getPreferredAlignment() {
        return LEFT;
    }
}