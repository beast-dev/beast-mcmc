/*
 * ToolbarOptions.java
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

package org.virion.jam.toolbar;

/**
 * @author rambaut
 *         Date: Oct 18, 2005
 *         Time: 10:23:01 PM
 */
public final class ToolbarOptions {

	public static final int ICON_AND_TEXT = 0;
	public static final int ICON_ONLY = 1;
	public static final int TEXT_ONLY = 2;

	public ToolbarOptions(int display, boolean smallSize) {
		this.display = display;
		this.smallSize = smallSize;
	}

	public int getDisplay() {
		return display;
	}

	public boolean getSmallSize() {
		return smallSize;
	}

	private int display = ICON_AND_TEXT;
	private boolean smallSize = false;
}
