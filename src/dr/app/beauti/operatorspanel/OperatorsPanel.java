/*
 * OperatorsPanel.java
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

package dr.app.beauti.operatorspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Operator;
import dr.app.beauti.types.OperatorSetType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.RealNumberCellEditor;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: OperatorsPanel.java,v 1.12 2005/07/11 14:07:25 rambaut Exp $
 */
public class OperatorsPanel extends BeautiPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = -3456667023451785854L;
    JScrollPane scrollPane = new JScrollPane();
    JTable operatorTable = null;
    OperatorTableModel operatorTableModel = null;

    JCheckBox autoOptimizeCheck = null;

    JComboBox operatorSetCombo = new JComboBox(new OperatorSetType[] {
            OperatorSetType.DEFAULT,
            OperatorSetType.FIXED_TREE_TOPOLOGY,
            OperatorSetType.NEW_TREE_MIX,
            OperatorSetType.ADAPTIVE_MULTIVARIATE,
            OperatorSetType.CUSTOM
    });

    public List<Operator> operators = new ArrayList<Operator>();

    private BeautiOptions options;
    BeautiFrame frame = null;

    public OperatorsPanel(BeautiFrame parent) {

        this.frame = parent;

        operatorTableModel = new OperatorTableModel();
        operatorTable = new JTable(operatorTableModel);

        operatorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        operatorTable.getTableHeader().setReorderingAllowed(false);
//        operatorTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        operatorTable.getColumnModel().getColumn(0).setMinWidth(40);
        operatorTable.getColumnModel().getColumn(0).setMaxWidth(40);

        operatorTable.getColumnModel().getColumn(1).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(1).setMinWidth(200);

        operatorTable.getColumnModel().getColumn(2).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(2).setMinWidth(140);

        operatorTable.getColumnModel().getColumn(3).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(3).setCellEditor(
                new RealNumberCellEditor(0, Double.POSITIVE_INFINITY));
        operatorTable.getColumnModel().getColumn(3).setMinWidth(40);

        operatorTable.getColumnModel().getColumn(4).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(4).setCellEditor(
                new RealNumberCellEditor(0, Double.MAX_VALUE));
        operatorTable.getColumnModel().getColumn(4).setMinWidth(40);

        operatorTable.getColumnModel().getColumn(5).setCellRenderer(
                new OperatorTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        operatorTable.getColumnModel().getColumn(5).setMinWidth(380);

        scrollPane = new JScrollPane(operatorTable,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

        scrollPane.setOpaque(false);

        autoOptimizeCheck = new JCheckBox("Auto Optimize");
        autoOptimizeCheck.setToolTipText("<html>This option will attempt to tune the operators<br>" +
                "to maximum efficiency. Turn off to tune the<br>" +
                "operators manually.</html>");
        autoOptimizeCheck.setOpaque(false);
        autoOptimizeCheck.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                options.autoOptimize = autoOptimizeCheck.isSelected();
                if (!autoOptimizeCheck.isSelected()) options.operatorAnalysis = true;
                frame.setDirty();
            }
        });

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setBorder(BorderFactory.createEmptyBorder());
        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        toolBar1.add(autoOptimizeCheck);
        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));
        final JLabel label = new JLabel("Operator mix: ");
        toolBar1.add(label);
        PanelUtils.setupComponent(operatorSetCombo);
        toolBar1.add(operatorSetCombo);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);


        operatorSetCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        options.operatorSetType = (OperatorSetType)operatorSetCombo.getSelectedItem();
                        operators = options.selectOperators();
                        operatorTableModel.fireTableDataChanged();
                        operatorsChanged();
                    }
                }
        );
    }

    public final void operatorsChanged() {
        frame.setDirty();
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        operatorSetCombo.setSelectedItem(options.operatorSetType);

        autoOptimizeCheck.setSelected(options.autoOptimize);
        operators = options.selectOperators();
        operatorTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
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
            Operator op = operators.get(row);
            switch (col) {
                case 0:
                    return op.isUsed();
                case 1:
                    return op.getName();
                case 2:
                    return op.getOperatorType();
                case 3:
                    if (op.isTunable()) {
                        return op.getTuning();
                    } else {
                        return "n/a";
                    }
                case 4:
                    return op.getWeight();
                case 5:
                    return op.getDescription();
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            Operator op = operators.get(row);
            switch (col) {
                case 0:
                    op.setUsed((Boolean) aValue);
                    operatorSetCombo.setSelectedItem(OperatorSetType.CUSTOM);
                    break;
                case 3:
                    op.setTuning((Double) aValue);
                    break;
                case 4:
                    op.setWeight((Double) aValue);
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

            Operator op = operators.get(row);

            switch (col) {
                case 0:// Check box
                    // if the paramter is fixed then 'in use' can't be turned on
                    editable = !op.isParameterFixed();
                    break;
                case 3:
                    editable = op.isUsed() && op.isTunable();
                    break;
                case 4:
                    editable = op.isUsed();
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

    class OperatorTableCellRenderer extends TableRenderer {

        public OperatorTableCellRenderer(int alignment, Insets insets) {
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

            Operator op = operators.get(aRow);
            if (!op.isUsed() && aColumn > 0) {
                renderer.setForeground(Color.gray);
            } else {
                renderer.setForeground(Color.black);
            }
            if (op.isParameterFixed()) {
               setToolTipText(
                       "This parameter is set to a fixed value. To turn \r" +
                       "this move on, select a prior in the Priors tab");
            }

            return this;
        }

    }

}