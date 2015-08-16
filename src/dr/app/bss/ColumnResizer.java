/*
 * ColumnResizer.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.bss;

import java.awt.Component;

import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class ColumnResizer {

	public static void adjustColumnPreferredWidths(JTable table) {

		TableColumnModel columnModel = table.getColumnModel();

		for (int col = 0; col < table.getColumnCount(); col++) {

			int maxwidth = 0;
			for (int row = 0; row < table.getRowCount(); row++) {

				TableCellRenderer rend = table.getCellRenderer(row, col);
				Object value = table.getValueAt(row, col);
				Component comp = rend.getTableCellRendererComponent(table,
						value, false, false, row, col);
				maxwidth = Math.max(comp.getPreferredSize().width, maxwidth);

			}// END: row loop

			TableColumn column = columnModel.getColumn(col);
			TableCellRenderer headerRenderer = column.getHeaderRenderer();
			if (headerRenderer == null) {
				headerRenderer = table.getTableHeader().getDefaultRenderer();
			}

			Object headerValue = column.getHeaderValue();
			Component headerComp = headerRenderer
					.getTableCellRendererComponent(table, headerValue, false,
							false, 0, col);
			maxwidth = Math.max(maxwidth, headerComp.getPreferredSize().width);
			column.setPreferredWidth(maxwidth);
		}// END: col loop
	}//END: adjustColumnPreferredWidths

}// END: class