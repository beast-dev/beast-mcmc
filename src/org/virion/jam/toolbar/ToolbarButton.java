/*
 * ToolbarButton.java
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
 *         Time: 10:09:21 PM
 */
public class ToolbarButton extends JButton implements ToolbarItem {

	public ToolbarButton(ToolbarAction action) {
		super(action);

		setHorizontalTextPosition(SwingConstants.CENTER);
		setVerticalTextPosition(SwingConstants.BOTTOM);
		putClientProperty("JButton.buttonType", "toolbar");
		setBorderPainted(false);

	    setToolTipText(action.getToolTipText());

		setDisabledIcon(action.getDisabledIcon());
		setPressedIcon(action.getPressedIcon());
	}

	public void setToolbarOptions(ToolbarOptions options) {
		switch (options.getDisplay()) {
			case ToolbarOptions.ICON_AND_TEXT:
				setText(action.getLabel());
				setIcon(action.getIcon());
				setDisabledIcon(action.getDisabledIcon());
				setPressedIcon(action.getPressedIcon());
				break;
			case ToolbarOptions.ICON_ONLY:
				setText(null);
				setIcon(action.getIcon());
				setDisabledIcon(action.getDisabledIcon());
				setPressedIcon(action.getPressedIcon());
				break;
			case ToolbarOptions.TEXT_ONLY:
				setText(action.getLabel());
				setIcon(null);
				setDisabledIcon(null);
				setPressedIcon(null);
				break;
		}
	}

	public void setAction(Action action) {
		super.setAction(action);
		if (action instanceof ToolbarAction) {
			this.action = (ToolbarAction)action;
		}
	}

	private ToolbarAction action;
}
