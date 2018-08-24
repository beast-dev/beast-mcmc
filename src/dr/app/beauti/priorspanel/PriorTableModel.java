/*
 * PriorTableModel.java
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

package dr.app.beauti.priorspanel;

import dr.app.beauti.options.Parameter;

import javax.swing.table.AbstractTableModel;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
class PriorTableModel extends AbstractTableModel {

    private static final long serialVersionUID = -8864178122484971872L;

    String[] columnNames = {"Parameter", "Prior", "Bound", "Description"};
    private PriorsPanel priorsPanel;

    public PriorTableModel(PriorsPanel priorsPanel) {
        this.priorsPanel = priorsPanel;
    }

    public int getColumnCount() {
        return columnNames.length;
    }

    public int getRowCount() {
        return priorsPanel.parameters.size();
    }

    public Object getValueAt(int row, int col) {
        Parameter param = priorsPanel.parameters.get(row);
        switch (col) {
            case 0:
                return param.getName();
            case 1:
                return (param.isLinked ? "Linked as [" + param.linkedName + "]" : param.priorType.getPriorString(param));
            case 2:
                return param.priorType.getPriorBoundString(param);
            case 3:
                return param.getDescription();
        }
        return null;
    }

    public String getColumnName(int column) {
        return columnNames[column];
    }

    public Class getColumnClass(int c) {
        return getValueAt(0, c).getClass();
    }

    public boolean isCellEditable(int row, int col) {
        Parameter param = priorsPanel.parameters.get(row);
        return col == 1 && !param.isPriorFixed;
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer();

        buffer.append(getColumnName(0));
        for (int j = 1; j < getColumnCount(); j++) {
            buffer.append("\t");
            buffer.append(getColumnName(j));
        }
        buffer.append("\n");

        for (int i = 0; i < getRowCount(); i++) {
            buffer.append(getValueAt(i, 0));
            for (int j = 1; j < getColumnCount(); j++) {
                buffer.append("\t");
                buffer.append(getValueAt(i, j));
            }
            buffer.append("\n");
        }

        return buffer.toString();
    }
}
