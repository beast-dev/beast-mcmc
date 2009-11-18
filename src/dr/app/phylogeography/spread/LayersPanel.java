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

package dr.app.phylogeography.spread;

import dr.app.phylogeography.builder.*;
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
 * @version $Id$
 */
public class LayersPanel extends JPanel implements Exportable {
    public final static BuilderFactory[] builderFactories = {
            DiscreteTreeBuilder.FACTORY
    };

    private JScrollPane scrollPane = new JScrollPane();
    private JTable layerTable = null;
    private LayerTableModel layerTableModel = null;

    private SpreadFrame frame = null;

    private LayerBuilderDialog layerBuilderDialog = null;

    private final SpreadDocument document;

    public LayersPanel(final SpreadFrame parent, final SpreadDocument document) {

        this.frame = parent;
        this.document = document;

        layerTableModel = new LayerTableModel();
        layerTable = new JTable(layerTableModel);

        layerTable.getTableHeader().setReorderingAllowed(false);
        layerTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

//        TableColumn col = layerTable.getColumnModel().getColumn(5);
//        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
//        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
//        col.setCellRenderer(comboBoxRenderer);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(layerTable);

        layerTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        layerTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    editSelection();
                }
            }
        });

        scrollPane = new JScrollPane(layerTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));

//        JButton button = new JButton(unlinkModelsAction);
//        unlinkModelsAction.setEnabled(false);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);

        Action addLayerAction = new AddLayerAction();
        Action removeLayerAction = new RemoveLayerAction();
        Action editLayerAction = new EditLayerAction();

        ActionPanel actionPanel1 = new ActionPanel(true);
        actionPanel1.setAddAction(addLayerAction);
        actionPanel1.setRemoveAction(removeLayerAction);
        actionPanel1.setActionAction(editLayerAction);

        removeLayerAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
        add(controlPanel1, BorderLayout.SOUTH);

        document.addListener(new SpreadDocument.Listener() {
            public void dataChanged() {
            }

            public void settingsChanged() {
                layerTableModel.fireTableDataChanged();
            }
        });
    }

    public void selectionChanged() {
        int[] selRows = layerTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
        frame.dataSelectionChanged(hasSelection);
    }

    public JComponent getExportableComponent() {
        return layerTable;
    }

    public void removeSelection() {
        int[] selRows = layerTable.getSelectedRows();
        Set<Builder> buildersToRemove = new HashSet<Builder>();
        for (int row : selRows) {
            buildersToRemove.add(document.getLayerBuilders().get(row));
        }
//
//        // TODO: would probably be a good idea to check if the user wants to remove the last partition
        document.getLayerBuilders().removeAll(buildersToRemove);
        document.fireSettingsChanged();
    }

    public void editSelection() {
        int selRow = layerTable.getSelectedRow();
        if (selRow >= 0) {
            Builder builder = document.getLayerBuilders().get(selRow);
            editSettings(builder);
        }
    }

    private void editSettings(Builder builder) {
        if (layerBuilderDialog == null) {
            layerBuilderDialog = new LayerBuilderDialog(frame);
        }

        int result = layerBuilderDialog.showDialog(builder);
        layerBuilderDialog.getBuilder(); // force update of builder settings
        if (result != JOptionPane.CANCEL_OPTION) {
            document.fireSettingsChanged();
        }
    }

    public void addLayer() {
        if (layerBuilderDialog == null) {
            layerBuilderDialog = new LayerBuilderDialog(frame);
        }

        int result = layerBuilderDialog.showDialog(builderFactories);
        if (result != JOptionPane.CANCEL_OPTION) {
            Builder builder = layerBuilderDialog.getBuilder();
            document.addLayerBuilder(builder);
        }
    }

    public void selectAll() {
        layerTable.selectAll();
    }

    class LayerTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Name", "Layer Type"};

        public LayerTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return document.getLayerBuilders().size();
        }

        public Object getValueAt(int row, int col) {
            Builder builder = document.getLayerBuilders().get(row);
            switch (col) {
                case 0:
                    return builder.getName();
                case 1:
                    return builder.getBuilderName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
            Builder builder = document.getLayerBuilders().get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        builder.setName(name);
                    }
                    break;
            }
            document.fireSettingsChanged();
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

    public class AddLayerAction extends AbstractAction {
        public AddLayerAction() {
            super("Add");
            setToolTipText("Use this button to create a new layer");
        }

        public void actionPerformed(ActionEvent ae) {
            addLayer();
        }
    }

    public class RemoveLayerAction extends AbstractAction {
        public RemoveLayerAction() {
            super("Remove");
            setToolTipText("Use this button to remove a selected layer from the table");
        }

        public void actionPerformed(ActionEvent ae) {
            removeSelection();
        }
    }

    public class EditLayerAction extends AbstractAction {
        public EditLayerAction() {
            super("Edit");
            setToolTipText("Use this button to edit a selected layer in the table");
        }

        public void actionPerformed(ActionEvent ae) {
            editSelection();
        }
    }


}