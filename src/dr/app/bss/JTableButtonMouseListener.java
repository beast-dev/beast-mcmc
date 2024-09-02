/*
 * JTableButtonMouseListener.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.bss;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JButton;
import javax.swing.JTable;

/**
 * @author Filip Bielejec
 */
public class JTableButtonMouseListener extends MouseAdapter {

	private final JTable table;

	public JTableButtonMouseListener(JTable table) {
		this.table = table;
	}// END: Constructor

	public void mouseClicked(MouseEvent e) {

		int column = table.getColumnModel().getColumnIndexAtX(e.getX());
		int row = e.getY() / table.getRowHeight();

		if (row < table.getRowCount() && row >= 0
				&& column < table.getColumnCount() && column >= 0) {

			Object value = table.getValueAt(row, column);
			if (value instanceof JButton) {

				((JButton) value).doClick();

			}// END: JButton check

		}// END: placement check
	}// END: mouseClicked

}// END: class