/*
 * ToolbarAction.java
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

import javax.swing.*;

/**
 * @author rambaut
 *         Date: Oct 18, 2005
 *         Time: 10:10:52 PM
 */
public abstract class ToolbarAction extends AbstractAction {

	protected ToolbarAction(String label, String toolTipText, Icon icon) {
		this(label, toolTipText, icon, null, null);
	}

	protected ToolbarAction(String label, String toolTipText, Icon icon, Icon disabledIcon, Icon pressedIcon) {
		super(label, icon);

		this.label = label;
		this.toolTipText = toolTipText;
		this.icon = icon;
		this.disabledIcon = disabledIcon;
		this.pressedIcon = pressedIcon;
	}

	public String getLabel() {
		return label;
	}

	public Icon getIcon() {
		return icon;
	}

	public Icon getDisabledIcon() {
		return disabledIcon;
	}

	public Icon getPressedIcon() {
		return pressedIcon;
	}

	public String getToolTipText() {
		return toolTipText;
	}

	private String label;
	private String toolTipText;
	private Icon icon;
	private Icon disabledIcon;
	private Icon pressedIcon;
}
