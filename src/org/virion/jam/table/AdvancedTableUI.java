/*
 * AdvancedTableUI.java
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

package org.virion.jam.table;

import javax.swing.plaf.basic.BasicTableUI;
import javax.swing.event.MouseInputListener;
import javax.swing.*;
import java.awt.event.MouseEvent;
import java.awt.*;

/**
 * @author rambaut
 *         Date: Oct 20, 2004
 *         Time: 10:16:52 PM
 */
public class AdvancedTableUI extends BasicTableUI {

    public void installUI(JComponent c) {
        super.installUI(c);

        if (org.virion.jam.mac.Utils.isMacOSX()) {
            c.setFont(new Font("Lucida Grande", Font.PLAIN, 9));
        }
    }

	protected MouseInputListener createMouseInputListener() {
		return new AdvancedTableUI.AdvancedMouseInputHandler();
	}

	class AdvancedMouseInputHandler extends MouseInputHandler {
		public void mousePressed(MouseEvent e) {
			Point origin = e.getPoint();
			int row = table.rowAtPoint(origin);
			int column = table.columnAtPoint(origin);
			if (row != -1 && column != -1) {
				if (table.isCellSelected(row, column)) {
					e.consume();
				}
			}

			super.mousePressed(e);
		}
	}
}
