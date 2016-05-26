/*
 * DefaultFileMenuFactory.java
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
