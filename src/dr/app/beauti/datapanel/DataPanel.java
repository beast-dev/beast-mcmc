/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.SnAPhyl.datapanel;

import dr.app.SnAPhyl.BeautiFrame;

import dr.app.beauti.BeautiPanel;
import dr.app.beauti.alignmentviewer.*;
import dr.app.beauti.datapanel.BeautiAlignmentBuffer;
import dr.app.beauti.options.*;
import dr.app.beauti.util.PanelUtils;
import dr.evolution.datatype.DataType;
import dr.evolution.alignment.Alignment;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;

import java.awt.*;
import java.awt.event.*;
import java.util.*;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class DataPanel extends BeautiPanel implements Exportable {

    JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    ShowAction showAction = new ShowAction();

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

        dataTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    showAlignment();
                }
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


            toolBar1.addSeparator();

            JButton button = new JButton(showAction);
            showAction.setEnabled(false);
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

    private void showAlignment() {

        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            JFrame frame = new JFrame();
            frame.setSize(800, 600);

            PartitionData partition = options.dataPartitions.get(row);
            Alignment alignment = partition.getAlignment();
            AlignmentViewer viewer = new AlignmentViewer();
            if (alignment.getDataType().getType() == DataType.NUCLEOTIDES) {
                viewer.setCellDecorator(new StateCellDecorator(new NucleotideDecorator(), false));
            } else if (alignment.getDataType().getType() == DataType.AMINO_ACIDS) {
                viewer.setCellDecorator(new StateCellDecorator(new AminoAcidDecorator(), false));
            } else {
                // no colouring
            }
            viewer.setAlignmentBuffer(new BeautiAlignmentBuffer(alignment));

            JPanel panel = new JPanel(new BorderLayout());
            panel.setOpaque(false);
            panel.add(viewer, BorderLayout.CENTER);

            JPanel infoPanel = new JPanel(new BorderLayout());
            infoPanel.setOpaque(false);
            panel.add(infoPanel, BorderLayout.SOUTH);

            frame.setContentPane(panel);
            frame.setVisible(true);
        }

    }

    private void fireDataChanged() {
        options.updateLinksBetweenPDPCMPSMPTMPTPP();
        options.updatePartitionClockTreeLinks();
        
        options.clockModelOptions.fixRateOfFirstClockPartition(); //TODO correct?

        frame.setDirty();
    }


    public void selectionChanged() {
        int[] selRows = dataTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
        frame.dataSelectionChanged(hasSelection);

        showAction.setEnabled(hasSelection);
//        unlinkAllAction.setEnabled(hasSelection);
//        linkAllAction.setEnabled(selRows != null && selRows.length > 1);
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;



        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {

    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void removeSelection() {
        int[] selRows = dataTable.getSelectedRows();
        Set<PartitionData> partitionsToRemove = new HashSet<PartitionData>();
        for (int row : selRows) {
            partitionsToRemove.add(options.dataPartitions.get(row));
        }

        // TODO: would probably be a good idea to check if the user wants to remove the last partition
        options.dataPartitions.removeAll(partitionsToRemove);

        

        if (options.dataPartitions.size() == 0) {
            // all data partitions removed so reset the taxa
            options.reset();
            frame.statusLabel.setText("");
            frame.setAllOptions();
            frame.getExportAction().setEnabled(false);
        }

        dataTableModel.fireTableDataChanged();

        fireDataChanged();
    }

    public void selectAll() {
        dataTable.selectAll();
    }
    

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Name", "FileName", "Taxa", "Sites", "Sequence Type"};

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
            PartitionData partition = options.dataPartitions.get(row);
            switch (col) {
                case 0:
                    return partition.getName();
                case 1:
                    return partition.getFileName();
                case 2:
                    return "" + partition.getTaxaCount();
                case 3:
                    return "" + partition.getSiteCount(); // sequence length
                case 4:
                    return partition.getAlignment().getDataType().getDescription();

                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
            PartitionData partition = options.dataPartitions.get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        partition.setName(name);
                    }
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

    public class ShowAction extends AbstractAction {
        public ShowAction() {
            super("Show");
            setToolTipText("Display the selected alignments");
        }

        public void actionPerformed(ActionEvent ae) {
            showAlignment();
        }
    }

}