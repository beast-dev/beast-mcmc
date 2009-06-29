/*
 * ClockModelPanel.java
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

package dr.app.beauti.modelsPanel;

import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.options.*;
import dr.app.beauti.*;
import dr.evolution.datatype.DataType;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.*;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: ClockModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class ClockModelPanel extends BeautiPanel implements Exportable {

    JTable dataTable = null;
    DataTableModel dataTableModel = null;
    
    BeautiFrame frame = null;

    BeautiOptions options = null;

    public ClockModelPanel(BeautiFrame parent) {

        this.frame = parent;

        dataTableModel = new DataTableModel();
        dataTable = new JTable(dataTableModel);

        dataTable.getTableHeader().setReorderingAllowed(false);
        dataTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = dataTable.getColumnModel().getColumn(3);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        col.setCellRenderer(comboBoxRenderer);
        
        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);
        
        JScrollPane scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));        
        add(scrollPane, BorderLayout.SOUTH);
    }

    
    
    private void fireDataChanged() {
        frame.setDirty();
    }

    private void modelsChanged() {        
        TableColumn col = dataTable.getColumnModel().getColumn(1);
        col.setCellEditor(new DefaultCellEditor(new JComboBox(ClockType.values())));
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;
        

        modelsChanged();

        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
        
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Clock Model Name", "Molecular Clock Model"};      

        public DataTableModel() {
        }

        public int getColumnCount() {           
            return columnNames.length;            
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.getPartitionClockModels().size();
        }

        public Object getValueAt(int row, int col) {
        	PartitionClockModel model = options.getPartitionClockModels().get(row);
            switch (col) {
                case 0:
                    return model.getName();
                case 1: 
                	return model.getClockType();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
        	PartitionClockModel model = options.getPartitionClockModels().get(row);
            switch (col) {
            	case 0:
            		String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                    	model.setName(name);
                    }
            		break;
            	case 1:
            		model.setClockType((ClockType) aValue);
                    break;
            }
            fireDataChanged();
        }

        public boolean isCellEditable(int row, int col) {
            return true;
        }


        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            if (getRowCount() == 0) {
                return Object.class;
            }
            return getValueAt(0, c).getClass();
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

}