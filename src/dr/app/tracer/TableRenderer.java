/*
 * TableRenderer.java
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

package dr.app.tracer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.DefaultTableCellRenderer;
import java.awt.*;


public class TableRenderer extends DefaultTableCellRenderer {
	/**
	 *
	 */
	private static final long serialVersionUID = -3815955688412224510L;
	Color bg1 = Color.white;
	Color bg2 = new Color(0xEE, 0xEE, 0xFF);
	int alignment;

	public TableRenderer(int alignment, int inset) {
		super();
		this.alignment = alignment;
		if (alignment == SwingConstants.LEFT) {
			setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, inset, 0, 0)));
		} else if (alignment == SwingConstants.RIGHT) {
			setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 0, 0, inset)));
		}
	}

	public Component getTableCellRendererComponent(
		JTable table, Object value, boolean isSelected, 
		boolean hasFocus, int row, int column) {

		setOpaque(true);
		setText(value.toString());


		// if cell is selected, set background color to default cell selection background color
		if (isSelected) {
			setBackground( table.getSelectionBackground() ) ;
			setForeground( table.getSelectionForeground() );
		} else {
			// otherwise, set cell background color to our custom color
			if (row % 2 == 0) {
				setBackground(bg1);
			} else {
				setBackground(bg2);
			}

			// set cell's foreground to default cell foreground color
			setForeground(table.getForeground());
		}

		// draw border on cell if it has focus
		//if (hasFocus) {
		//	setBorder( UIManager.getBorder("Table.focusCellHighlightBorder") );
		//}

		// position cell text at center
		setHorizontalAlignment(alignment);
		return this;
	}
}
