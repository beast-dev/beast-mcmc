/*
 * JTableComboBoxCellRenderer.java
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

import javax.swing.DefaultListCellRenderer;
import javax.swing.JComboBox;
import javax.swing.JList;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;


@SuppressWarnings("serial")
public class JTableComboBoxCellRenderer extends JComboBox implements
		TableCellRenderer {

	private int columnIndex;
	
	private DefaultListCellRenderer comboBoxRenderer = new DefaultListCellRenderer() {

		@Override
		public Component getListCellRendererComponent(JList list, Object value,
				int index, boolean isSelected, boolean cellHasFocus) {

			super.getListCellRendererComponent(list, value, index, isSelected,
					cellHasFocus);

			if (value != null) {

				switch (columnIndex) {

				case PartitionTableModel.DATA_INDEX:

					TreesTableRecord record = (TreesTableRecord) value;
					this.setText(record.getName());
					break;

				default:
					break;

				}// END: switch				
				
			}//END: null check

			return this;
		}

	};

	public JTableComboBoxCellRenderer(int columnIndex) {

		super();
		setOpaque(true);
		this.setRenderer(comboBoxRenderer);
		this.columnIndex = columnIndex;
		
	}// END: Constructor

	public Component getTableCellRendererComponent(JTable table, Object value,
			boolean isSelected, boolean hasFocus, int row, int column) {

		if (isSelected) {

			this.setForeground(table.getSelectionForeground());
			this.setBackground(table.getSelectionBackground());

		} else {

			this.setForeground(table.getForeground());
			this.setBackground(table.getBackground());

		}

		// Select the current value
		setSelectedItem(value);

		if (value != null) {
			removeAllItems();
			addItem(value);
		}

		return this;
	}
}// END: JTableComboBoxCellRenderer class