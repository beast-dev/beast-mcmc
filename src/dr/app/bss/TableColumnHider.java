/*
 * TableColumnHider.java
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

import java.util.HashMap;
import java.util.Map;

import javax.swing.JTable;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class TableColumnHider {

	private JTable table;
	private TableColumnModel tcm;
	private Map<String, IndexedColumn> hidden = new HashMap<String, IndexedColumn>();

	public TableColumnHider(JTable table) {
		this.table = table;
		this.tcm = this.table.getColumnModel();
	}// END: Constructor

	public void hide(String columnName) {
		int index = tcm.getColumnIndex(columnName);
		TableColumn column = tcm.getColumn(index);
		IndexedColumn ic = new IndexedColumn(index, column);
		if (hidden.put(columnName, ic) != null) {
			throw new IllegalArgumentException("Duplicate column name.");
		}
		tcm.removeColumn(column);
	}// END: hide

	public void show(String columnName) {
		IndexedColumn ic = hidden.remove(columnName);
		if (ic != null) {
			tcm.addColumn(ic.column);
			int lastColumn = tcm.getColumnCount() - 1;
			if (ic.index < lastColumn) {
				tcm.moveColumn(lastColumn, ic.index);
			}
		}
	}// END: show

	private static class IndexedColumn {

		private Integer index;
		private TableColumn column;

		public IndexedColumn(Integer index, TableColumn column) {
			this.index = index;
			this.column = column;
		}
	}// END: IndexedColumn

}// END: class
