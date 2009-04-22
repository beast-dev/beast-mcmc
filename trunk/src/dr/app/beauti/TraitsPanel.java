/*
 * SamplesPanel.java
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
import dr.evolution.util.*;
import dr.gui.table.TableSorter;
import org.virion.jam.framework.Exportable;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;
import org.virion.jam.panels.ActionPanel;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class TraitsPanel extends BeautiPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = 5283922195494563924L;

    JTable traitsTable = null;
    TraitsTableModel traitsTableModel = null;

    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    ClearLocationsAction clearLocationsAction = new ClearLocationsAction();
    GuessLocationsAction guessLocationsAction = new GuessLocationsAction();

    BeautiFrame frame = null;

    BeautiOptions options = null;

    double[] heights = null;

    GuessDatesDialog guessDatesDialog = null;

    public TraitsPanel(BeautiFrame parent) {

        this.frame = parent;

        traitsTableModel = new TraitsTableModel();
        TableSorter sorter = new TableSorter(traitsTableModel);
        traitsTable = new JTable(sorter);
        sorter.setTableHeader(traitsTable.getTableHeader());

        traitsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        traitsTable.getColumnModel().getColumn(0).setPreferredWidth(80);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(traitsTable);


        dataTableModel = new DataTableModel();
        sorter = new TableSorter(dataTableModel);
        dataTable = new JTable(sorter);
        sorter.setTableHeader(dataTable.getTableHeader());

        dataTable.getTableHeader().setReorderingAllowed(false);
        dataTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        dataTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        dataTable.getColumnModel().getColumn(0).setPreferredWidth(80);

        TableColumn col = dataTable.getColumnModel().getColumn(1);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        col.setCellRenderer(comboBoxRenderer);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        JScrollPane scrollPane1 = new JScrollPane(traitsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane1.setOpaque(false);

        JScrollPane scrollPane2 = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane2.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);

        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        JButton button = new JButton(clearLocationsAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);
        button = new JButton(guessLocationsAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);
//        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addTraitAction);
        actionPanel1.setRemoveAction(removeTraitAction);

        removeTraitAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane1, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(toolBar1, BorderLayout.NORTH);
        panel2.add(scrollPane2, BorderLayout.CENTER);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                panel1, panel2);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.CENTER);
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        setupTable();

        dataTableModel.fireTableDataChanged();
    }

    private void setupTable() {
        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
//        options.datesUnits = unitsCombo.getSelectedIndex();
//        options.datesDirection = directionCombo.getSelectedIndex();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void selectionChanged() {
        // nothing to do
    }

    public void clearLocations() {
//        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
//            java.util.Date origin = new java.util.Date(0);
//
//            double d = 0.0;
//
//            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
//            options.taxonList.getTaxon(i).setAttribute("date", date);
//        }
//
//        // adjust the dates to the current timescale...
//        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    public void guessLocations() {

//        if (guessDatesDialog == null) {
//            guessDatesDialog = new GuessDatesDialog(frame);
//        }
//
//        int result = guessDatesDialog.showDialog();
//
//        if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
//            return;
//        }
//
//        DateGuesser guesser = options.dateGuesser;
//
//        guesser.guessDates = true;
//        guessDatesDialog.setupGuesser(guesser);
//
//        String warningMessage = null;
//
//        guesser.guessDates(options.taxonList);
//
//        if (warningMessage != null) {
//            JOptionPane.showMessageDialog(this, "Warning: some dates may not be set correctly - \n" + warningMessage,
//                    "Error guessing dates",
//                    JOptionPane.WARNING_MESSAGE);
//        }
//
//        // adjust the dates to the current timescale...
//        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    private void addLocation() {

    }

    private void removeLocation() {

    }

    public class ClearLocationsAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = -7281309694753868635L;

        public ClearLocationsAction() {
            super("Clear Locations");
            setToolTipText("Use this tool to remove sampling locations from each taxon");
        }

        public void actionPerformed(ActionEvent ae) {
            clearLocations();
        }
    }

    public class GuessLocationsAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = 8514706149822252033L;

        public GuessLocationsAction() {
            super("Guess Locations");
            setToolTipText("Use this tool to guess the sampling locations from the taxon labels");
        }

        public void actionPerformed(ActionEvent ae) {
            guessLocations();
        }
    }

    AbstractAction addTraitAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            addLocation();
        }
    };


    AbstractAction removeTraitAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            removeLocation();
        }
    };


    class TraitsTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Trait", "Type"};

        public TraitsTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            if (options.traits == null) return 0;

            return options.traits.size();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.traits.get(row);
                case 1:
                    return options.traitTypes.get(options.traits.get(row));
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            if (col == 0) {
                options.traits.set(row, aValue.toString());
            } else if (col == 1) {
                options.traitTypes.put(options.traits.get(row), aValue.getClass());
            }
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return true;
            return false;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
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

    class DataTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Taxon", "Value"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            if (options.taxonList == null) return 0;

            return options.taxonList.getTaxonCount();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.taxonList.getTaxonId(row);
                case 1:
                    Location location = options.taxonList.getTaxon(row).getLocation();
                    if (location != null) {
                        return location.getId();
                    } else {
                        return "-";
                    }
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            if (col == 0) {
                options.taxonList.getTaxon(row).setId(aValue.toString());
            } else if (col == 1) {
                Location location = options.taxonList.getTaxon(row).getLocation();
                if (location != null) {
                    options.taxonList.getTaxon(row).setLocation(location);
                }
            }
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return false;
            if (col == 1) return true;
            return false;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
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