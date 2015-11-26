/*
 * JTableComboBoxCellEditor.java
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

import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JTable;


@SuppressWarnings("serial")
public class JTableComboBoxCellEditor extends DefaultCellEditor {

	private PartitionDataList dataList = null;
	
	public JTableComboBoxCellEditor(PartitionDataList dataList) {
		super(new JComboBox());
		this.dataList = dataList;
	}// END: Constructor
	
	@SuppressWarnings("unchecked")
	public Component getTableCellEditorComponent(JTable table,
			Object value, boolean isSelected, int row, int column) {

		((JComboBox)editorComponent).removeAllItems();

		if (column == PartitionTableModel.DATA_INDEX) {

			for (TreesTableRecord record : dataList.recordsList) {
				((JComboBox)editorComponent).addItem(record);
			}// END: fill loop
			
		} else if (column == PartitionTableModel.DATA_TYPE_INDEX) {

			for (String dataType : PartitionData.dataTypes) {
				((JComboBox)editorComponent).addItem(dataType);
			}// END: fill loop

		} else {

			// do nothing

		}// END: column check

		((JComboBox)editorComponent).setSelectedItem(value);
		delegate.setValue(value);

		return editorComponent;
	}// END: getTableCellEditorComponent

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: setDataList
	
}// END: JTableComboBoxCellEditor class