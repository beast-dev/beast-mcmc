/*
 * DataPanel.java
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

package dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.DataPartition;
import dr.app.beauti.options.PartitionModel;
import dr.evolution.datatype.DataType;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class DataPanel extends JPanel implements Exportable {

    JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    UnlinkModelsAction unlinkModelsAction = new UnlinkModelsAction();
    LinkModelsAction linkModelAction = new LinkModelsAction();

    SelectModelDialog selectModelDialog = null;

    BeautiFrame frame = null;

    BeautiOptions options = null;

    public DataPanel(BeautiFrame parent, Action importDataAction, Action removeDataAction) {

        this.frame = parent;

        dataTableModel = new DataTableModel();
        dataTable = new JTable(dataTableModel);

        dataTable.getTableHeader().setReorderingAllowed(false);
        dataTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });


        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        JButton button = new JButton(unlinkModelsAction);
        unlinkModelsAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(linkModelAction);
        linkModelAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(importDataAction);
        actionPanel1.setRemoveAction(removeDataAction);

        removeDataAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
        add(controlPanel1, BorderLayout.SOUTH);
    }


    private void fireDataChanged() {
        frame.dataChanged();
    }

    private void modelsChanged() {
        Object[] modelArray = options.getPartitionModels().toArray();
        TableColumn col = dataTable.getColumnModel().getColumn(4);

        col.setCellEditor(new DefaultCellEditor(new JComboBox(modelArray)));
        //col.setCellRenderer(new ComboBoxRenderer(modelArray));
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

    public void selectionChanged() {
        int[] selRows = dataTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
        frame.dataSelectionChanged(hasSelection);
        unlinkModelsAction.setEnabled(hasSelection);
        linkModelAction.setEnabled(selRows != null && selRows.length > 1);
    }

    public void removeSelection() {
        int[] selRows = dataTable.getSelectedRows();
        Set<DataPartition> partitionsToRemove = new HashSet<DataPartition>();
        for (int row : selRows) {
            partitionsToRemove.add(options.dataPartitions.get(row));
        }

        options.dataPartitions.removeAll(partitionsToRemove);
        dataTableModel.fireTableDataChanged();

        fireDataChanged();
    }

    public void unlinkModels() {
        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            DataPartition partition = options.dataPartitions.get(row);
            if (!partition.getPartitionModel().getName().equals(partition.getName())) {
                PartitionModel model = new PartitionModel(partition);
                options.addPartitionModel(model);
                partition.setPartitionModel(model);
            }
        }

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void linkModels() {
        int[] selRows = dataTable.getSelectedRows();
        DataType dateType = null;
        for (int row : selRows) {
            DataPartition partition = options.dataPartitions.get(row);
            if (dateType == null) {
                dateType = partition.getPartitionModel().dataType;
            } else {
                if (partition.getPartitionModel().dataType != dateType) {
                    JOptionPane.showMessageDialog(this, "Can only link the models for data partitions \n" +
                            "of the same data type (e.g., nucleotides)",
                            "Unable to link models",
                            JOptionPane.ERROR_MESSAGE);

                }
            }
        }

        java.util.List<PartitionModel> models = options.getPartitionModels(dateType);
        Object[] modelArray = models.toArray();

        if (selectModelDialog == null) {
            selectModelDialog = new SelectModelDialog(frame);
        }

        int result = selectModelDialog.showDialog(modelArray);
        if (result != JOptionPane.CANCEL_OPTION) {
            PartitionModel model = selectModelDialog.getModel();
            if (selectModelDialog.getMakeCopy()) {
                model = new PartitionModel(selectModelDialog.getName(), model);
                options.addPartitionModel(model);
            }

            for (int row : selRows) {
                DataPartition partition = options.dataPartitions.get(row);
                partition.setPartitionModel(model);
            }
        }

        fireDataChanged();
        repaint();
    }

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Name", "FileName", "Sites", "Sequence Type", "Partition Model"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.dataPartitions.size();
        }

        public Object getValueAt(int row, int col) {
            DataPartition partition = options.dataPartitions.get(row);
            switch (col) {
                case 0:
                    return partition.getName();
                case 1:
                    return partition.getFileName();
                case 2:
                    return "" + (partition.getToSite() - partition.getFromSite() + 1);
                case 3:
                    return partition.getAlignment().getDataType().getDescription();
                case 4:
                    return partition.getPartitionModel().getName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
            DataPartition partition = options.dataPartitions.get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        partition.setName(name);
                    }
                    break;
                case 4:
                    partition.setPartitionModel((PartitionModel)aValue);
                    break;
            }
            fireDataChanged();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            switch (col) {
                case 0:// name
                    editable = true;
                    break;
                case 4:// model selection menu
                    editable = true;
                    break;
                default:
                    editable = false;
            }

            return editable;
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

    public class ComboBoxRenderer extends JComboBox implements TableCellRenderer {
        public ComboBoxRenderer(Object[] items) {
            super(items);
        }

        public Component getTableCellRendererComponent(JTable table, Object value,
                                                       boolean isSelected, boolean hasFocus,
                                                       int row, int column) {

            if (isSelected) {
                this.setForeground(table.getSelectionForeground());
                this.setBackground(table.getSelectionBackground());
            } else {
                this.setForeground(table.getForeground());
                this.setBackground(table.getBackground());
            }

            setSelectedItem(value);
            return this;
        }
    }

    public class UnlinkModelsAction extends AbstractAction {
        public UnlinkModelsAction() {
            super("Unlink Models");
            setToolTipText("Use this tool to use a different model for each selected data partition");
        }

        public void actionPerformed(ActionEvent ae) {
            unlinkModels();
        }
    }


    public class LinkModelsAction extends AbstractAction {
        public LinkModelsAction() {
            super("Link Models");
            setToolTipText("Use this tool to set all the selected partitions to the same model");
        }

        public void actionPerformed(ActionEvent ae) {
            linkModels();
        }
    }
}