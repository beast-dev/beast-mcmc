/*
 * DefaultHelpMenuFactory.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.framework;

import javax.swing.*;
import java.awt.event.KeyEvent;

/**
 * @author rambaut
 *         Date: Dec 26, 2004
 *         Time: 11:02:20 AM
 */
public class DefaultHelpMenuFactory implements MenuFactory {
	public String getMenuName() {
		return "Help";
	}

	public void populateMenu(JMenu menu, AbstractFrame frame) {

        menu.setMnemonic('H');

		JMenuItem item;

		Application application = Application.getApplication();

		item = new JMenuItem(application.getAboutAction());
		menu.add(item);

		if (frame.getHelpAction() != null) {
			item = new JMenuItem(frame.getHelpAction());
			item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_SLASH, MenuBarFactory.MENU_MASK));
			menu.add(item);

			menu.addSeparator();
		}

		if (application.getHelpAction() != null) {
			item = new JMenuItem(application.getHelpAction());
			menu.add(item);

			menu.addSeparator();
		}

		if (application.getWebsiteAction() != null) {
			item = new JMenuItem("Website");
		}
	}

	public int getPreferredAlignment() {
		return RIGHT;
	}
}
