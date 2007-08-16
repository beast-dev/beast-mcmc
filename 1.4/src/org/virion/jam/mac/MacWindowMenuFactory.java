/*
 * MacWindowMenuFactory.java
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

package org.virion.jam.mac;

import org.virion.jam.framework.MenuFactory;
import org.virion.jam.framework.AbstractFrame;
import org.virion.jam.framework.Application;
import org.virion.jam.framework.MenuBarFactory;

import javax.swing.*;
import java.awt.event.KeyEvent;

/**
 * @author rambaut
 *         Date: Dec 26, 2004
 *         Time: 11:03:39 AM
 */
public class MacWindowMenuFactory implements MenuFactory {
	public String getMenuName() {
		return "Window";
	}

	public void populateMenu(JMenu menu, AbstractFrame frame) {

		Application application = Application.getApplication();

		JMenuItem item;

		item = new JMenuItem(frame.getZoomWindowAction());
		menu.add(item);

		item = new JMenuItem(frame.getMinimizeWindowAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, MenuBarFactory.MENU_MASK));
		menu.add(item);

	}

	public int getPreferredAlignment() {
		return RIGHT;
	}
}
