/*
 * AncestralStatesPanel.java
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

package dr.app.beauti.ancestralstatespanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.gui.table.TableEditorStopper;
import dr.evolution.datatype.DataType;
import jam.framework.Exportable;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class AncestralStatesPanel extends BeautiPanel implements Exportable {

    private static final long serialVersionUID = 2778103564318492601L;

    private static final int MINIMUM_TABLE_WIDTH = 140;

    JTable partitionTable = null;
    PartitionTableModel partitionTableModel = null;
    BeautiOptions options = null;

    JPanel optionsPanelParent;
    AbstractPartitionData currentPartition = null;
    Map<AbstractPartitionData, AncestralStatesOptionsPanel> optionsPanels = new HashMap<AbstractPartitionData, AncestralStatesOptionsPanel>();
    TitledBorder optionsBorder;

    BeautiFrame frame = null;
    boolean settingOptions = false;

    public AncestralStatesPanel(BeautiFrame parent) {

        super();

        this.frame = parent;

        partitionTableModel = new PartitionTableModel();
        partitionTable = new JTable(partitionTableModel);

        partitionTable.getTableHeader().setReorderingAllowed(false);
        partitionTable.getTableHeader().setResizingAllowed(false);
//        modelTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        final TableColumnModel model = partitionTable.getColumnModel();
        final TableColumn tableColumn0 = model.getColumn(0);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(partitionTable);

        partitionTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        partitionTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        JScrollPane scrollPane = new JScrollPane(partitionTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setOpaque(false);

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.setOpaque(false);
        panel.add(scrollPane, BorderLayout.CENTER);
        panel.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, 0));

        optionsPanelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        optionsPanelParent.setOpaque(false);
        optionsBorder = new TitledBorder(null, "Ancestral state options:", TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.ABOVE_TOP);
        optionsPanelParent.setBorder(optionsBorder);

        setCurrentPartition(null);

        JScrollPane scrollPane2 = new JScrollPane(optionsPanelParent, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane2.setOpaque(false);
        scrollPane2.setBorder(null);
        scrollPane2.getViewport().setOpaque(false);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel, scrollPane2);
        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.CENTER);
    }

    private void resetPanel() {
        if (!options.hasData()) {
            currentPartition = null;
            optionsPanels.clear();
            optionsPanelParent.removeAll();
        }
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        resetPanel();

        settingOptions = true;

        int selRow = partitionTable.getSelectedRow();
        partitionTableModel.fireTableDataChanged();
        if (options.getDataPartitions().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            partitionTable.getSelectionModel().setSelectionInterval(selRow, selRow);

            setCurrentPartition(options.getDataPartitions().get(selRow));
        }

        AncestralStatesOptionsPanel panel = optionsPanels.get(currentPartition);
        if (panel != null) {
            panel.setupPanel();
        }

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
    }

    public void fireModelChanged() {
        frame.setDirty();
    }

    private void selectionChanged() {
        if (settingOptions) return;

        int selRow = partitionTable.getSelectedRow();

        if (selRow >= options.getDataPartitions().size()) {
            selRow = 0;
            partitionTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        if (selRow >= 0) {
            setCurrentPartition(options.getDataPartitions().get(selRow));
//            frame.modelSelectionChanged(!isUsed(selRow));
        }
    }

    /**
     * Sets the current partition that this panel is displaying
     *
     * @param partition the new partition to display
     */
    private void setCurrentPartition(AbstractPartitionData partition) {
        if (partition != null) {
            if (currentPartition != null) optionsPanelParent.removeAll();

            AncestralStatesOptionsPanel panel = optionsPanels.get(partition);
            if (panel == null) {
                panel = new AncestralStatesOptionsPanel(this, options, partition);
                optionsPanels.put(partition, panel);
            }

            currentPartition = partition;

            panel.setupPanel();

            optionsPanelParent.add(panel);

            updateBorder();
        }
    }

    private void updateBorder() {

        String title;

        switch (currentPartition.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                title = "Nucleotide";
                break;
            case DataType.AMINO_ACIDS:
                title = "Amino Acid";
                break;
            case DataType.TWO_STATES:
                title = "Binary";
                break;
            case DataType.GENERAL:
                title = "Discrete Traits";
                break;
            case DataType.CONTINUOUS:
                title = "Continuous Traits";
                break;
            case DataType.MICRO_SAT:
                title = "Microsatellite";
                break;
            default:
                throw new IllegalArgumentException("Unsupported data type");

        }
        optionsBorder.setTitle(title + " Partition - " + currentPartition.getName());
        repaint();
    }

    public JComponent getExportableComponent() {
        return this;
    }

    class PartitionTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Partition"};

        public PartitionTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.getDataPartitions().size();
        }

        public Object getValueAt(int row, int col) {
            AbstractPartitionData partition = options.getDataPartitions().get(row);
            switch (col) {
                case 0:
                    return partition.getName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public boolean isCellEditable(int row, int col) {
            return false;
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