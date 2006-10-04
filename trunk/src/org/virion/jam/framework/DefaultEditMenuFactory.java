/*
 * DefaultEditMenuFactory.java
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
 *         Time: 11:01:40 AM
 */
public class DefaultEditMenuFactory implements MenuFactory {
	public String getMenuName() {
		return "Edit";
	}

	public void populateMenu(JMenu menu, AbstractFrame frame) {

		JMenuItem item;

        menu.setMnemonic('E');

		item = new JMenuItem(frame.getCutAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X, MenuBarFactory.MENU_MASK));
		menu.add(item);

		item = new JMenuItem(frame.getCopyAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, MenuBarFactory.MENU_MASK));
		menu.add(item);

		item = new JMenuItem(frame.getPasteAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, MenuBarFactory.MENU_MASK));
		menu.add(item);

		item = new JMenuItem(frame.getDeleteAction());
		menu.add(item);

		item = new JMenuItem(frame.getSelectAllAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, MenuBarFactory.MENU_MASK));
		menu.add(item);

		menu.addSeparator();

		item = new JMenuItem(frame.getFindAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, MenuBarFactory.MENU_MASK));
		menu.add(item);

	}

	public int getPreferredAlignment() {
		return LEFT;
	}
}
