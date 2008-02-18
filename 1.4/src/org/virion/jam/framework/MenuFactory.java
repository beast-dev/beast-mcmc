/*
 * MenuFactory.java
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

public interface MenuFactory {

	public final static int LEFT = 0;
	public final static int CENTER = 1;
	public final static int RIGHT = 2;

	/**
	 * Give the name of this menu. If multiple MenuFactories are
	 * registered with the same name, then these will be appended
	 * into a single actual menu.
	 */
	String getMenuName();

	/**
	 * This method should populate the menu with menu items. Reference
	 * can be made to the frame in order to get Actions.
	 * @param menu
	 * @param frame
	 */
    void populateMenu(JMenu menu, AbstractFrame frame);

	/**
	 * Returns the preferred alignment of the menu in the menu bar. This
	 * should be one of MenuFactory.LEFT, MenuFactory.CENTER or MenuFactory.RIGHT.
	 * @return the alignment
	 */
	int getPreferredAlignment();
}