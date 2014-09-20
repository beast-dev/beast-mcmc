/*
 * OperatorsPanel.java
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

package dr.app.oldbeauti;

import dr.app.gui.table.RealNumberCellEditor;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.ArrayList;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: OperatorsPanel.java,v 1.12 2005/07/11 14:07:25 rambaut Exp $
 */
public class OperatorsPanel extends JPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = -3456667023451785854L;
    JScrollPane scrollPane = new JScrollPane();
    JTable operatorTable = null;
    OperatorTableModel operatorTableModel = null;

    JCheckBox autoOptimizeCheck = null;

    public ArrayList operators = new ArrayList();

    BeautiFrame frame = null;

    public OperatorsPanel(BeautiFrame parent) {

        this.frame = parent;

        operatorTableModel = new OperatorTableModel();
        operatorTable = new JTable(operatorTableModel);

        operatorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        operatorTable.getTableHeader().setReorderingAllowed(false);
        operatorTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

//		operatorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_LAST_COLUMN);

        operatorTable.getColumnModel().getColumn(0).setPreferredWidth(50);

        operatorTable.getColumnModel().getColumn(1).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(1).setPreferredWidth(180);

        operatorTable.getColumnModel().getColumn(2).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(2).setPreferredWidth(140);

        operatorTable.getColumnModel().getColumn(3).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(3).setCellEditor(
                new RealNumberCellEditor(0, Double.POSITIVE_INFINITY));
        operatorTable.getColumnModel().getColumn(3).setPreferredWidth(50);

        operatorTable.getColumnModel().getColumn(4).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(4).setCellEditor(
                new RealNumberCellEditor(0, Double.MAX_VALUE));
        operatorTable.getColumnModel().getColumn(4).setPreferredWidth(50);

        operatorTable.getColumnModel().getColumn(5).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(5).setPreferredWidth(400);

        scrollPane = new JScrollPane(operatorTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

        scrollPane.setOpaque(false);

        autoOptimizeCheck = new JCheckBox("Auto Optimize - This option will attempt to tune the operators to maximum efficiency. Turn off to tune the operators manually.");
        autoOptimizeCheck.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        toolBar1.add(autoOptimizeCheck);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        add(toolBar1, "North");
        add(scrollPane, "Center");
    }

    public final void operatorsChanged() {
        frame.operatorsChanged();
    }

    public void setOptions(BeautiOptions options) {
        autoOptimizeCheck.setSelected(options.autoOptimize);

        operators = options.selectOperators();

        operatorTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {

        options.autoOptimize = autoOptimizeCheck.isSelected();

    }

    public JComponent getExportableComponent() {
        return operatorTable;
    }

    class OperatorTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -575580804476182225L;
        String[] columnNames = {"In use", "Operates on", "Type", "Tuning", "Weight", "Description"};

        public OperatorTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return operators.size();
        }

        public Object getValueAt(int row, int col) {
            BeastGenerator.Operator op = (BeastGenerator.Operator) operators.get(row);
            switch (col) {
                case 0:
                    return op.inUse;
                case 1:
                    return op.name;
                case 2:
                    return op.type;
                case 3:
                    if (op.isTunable()) {
                        return op.tuning;
                    } else {
                        return "n/a";
                    }
                case 4:
                    return op.weight;
                case 5:
                    return op.getDescription();
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            BeastGenerator.Operator op = (BeastGenerator.Operator) operators.get(row);
            switch (col) {
                case 0:
                    op.inUse = (Boolean) aValue;
                    break;
                case 3:
                    op.tuning = (Double) aValue;
                    op.tuningEdited = true;
                    break;
                case 4:
                    op.weight = (Double) aValue;
                    break;
            }
            operatorsChanged();
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            BeastGenerator.Operator op = (BeastGenerator.Operator) operators.get(row);

            switch (col){
                case 0:// Check box
                    editable = true;
                    break;
                case 3:
                    editable = op.inUse && op.isTunable();
                    break;
                case 4:
                    editable = op.inUse;
                    break;
                default:
                    editable = false;
            }

            return editable;
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

    class OperatorTableCellRenderer extends TableRenderer{

        public OperatorTableCellRenderer( int alignment, Insets insets){
            super(alignment, insets);
        }

        public Component getTableCellRendererComponent(JTable aTable,
                                                       Object value,
                                                       boolean aIsSelected,
                                                       boolean aHasFocus,
                                                       int aRow, int aColumn) {

            if (value == null) return this;

            Component renderer = super.getTableCellRendererComponent(aTable,
                    value,
                    aIsSelected,
                    aHasFocus,
                    aRow, aColumn);

            BeastGenerator.Operator op = (BeastGenerator.Operator) operators.get(aRow);
            if (! op.inUse && aColumn > 0)
                renderer.setForeground(Color.gray);
            else
                renderer.setForeground(Color.black);
            return this;
        }

    }
}
