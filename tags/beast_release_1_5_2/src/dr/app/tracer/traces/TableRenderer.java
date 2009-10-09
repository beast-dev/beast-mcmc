/*
 * TableRenderer.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.tracer.traces;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.DefaultTableCellRenderer;
import java.awt.*;


class TableRenderer extends DefaultTableCellRenderer {
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
