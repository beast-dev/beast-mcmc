/*
 * StatusBar.java
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

package org.virion.jam.panels;

import javax.swing.*;
import java.awt.*;

/**
 * @author rambaut
 *         Date: Oct 12, 2004
 *         Time: 12:18:09 AM
 */
public class StatusBar extends StatusPanel {
	public StatusBar(String initialText) {
		super(initialText);

		setBorder(BorderFactory.createCompoundBorder(
		    BorderFactory.createMatteBorder(0, 0, 1, 0, Color.gray),
		    BorderFactory.createEmptyBorder(2, 12, 2, 12)));
  //      panel.setBackground(new Color(0.0F, 0.0F, 0.0F, 0.05F));

	}

	public void paintComponent(Graphics g) {
	    super.paintComponent(g);
	    g.setColor(new Color(0.0F, 0.0F, 0.0F, 0.05F));
	    g.fillRect(0, 0, getWidth(), getHeight());
	}

}
