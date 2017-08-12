/*
 * TraitsPanel.java
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

package dr.app.beauti.traitspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.datapanel.DataPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class TraitsPanel extends BeautiPanel implements Exportable {

    private static final long serialVersionUID = 5283922195494563924L;

    private static final int MINIMUM_TABLE_WIDTH = 140;

    private static final String ADD_TRAITS_TOOLTIP = "<html>Define a new trait for the current taxa</html>";
    private static final String IMPORT_TRAITS_TOOLTIP = "<html>Import one or more traits for these taxa from a tab-delimited<br>" +
            "file. Taxa should be in the first column and the trait names<br>" +
            "in the first row</html>";
    private static final String GUESS_TRAIT_VALUES_TOOLTIP = "<html>This attempts to extract values for this trait that are<br>" +
            "encoded in the names of the selected taxa.</html>";
    private static final String SET_TRAIT_VALUES_TOOLTIP = "<html>This sets values for this trait for all<br>" +
            "the selected taxa.</html>";
    private static final String CLEAR_TRAIT_VALUES_TOOLTIP = "<html>This clears all the values currently assigned to taxa for<br>" +
            "this trait.</html>";
    private static final String CREATE_TRAIT_PARTITIONS_TOOLTIP = "<html>Create a data partition for the selected traits.</html>";

    public final JTable traitsTable;
    private final TraitsTableModel traitsTableModel;

    private final JTable dataTable;
    private final DataTableModel dataTableModel;

    private final BeautiFrame frame;
    private final DataPanel dataPanel;

    private BeautiOptions options = null;

    private TraitData currentTrait = null; // current trait

    private CreateTraitDialog createTraitDialog = null;
    private GuessTraitDialog guessTraitDialog = null;
    private TraitValueDialog traitValueDialog = null;

    AddTraitAction addTraitAction = new AddTraitAction();
    Action importTraitsAction;
    CreateTraitPartitionAction createTraitPartitionAction = new CreateTraitPartitionAction();
    GuessTraitsAction guessTraitsAction = new GuessTraitsAction();
    SetValueAction setValueAction = new SetValueAction();

    public TraitsPanel(BeautiFrame parent, DataPanel dataPanel, Action importTraitsAction) {

        this.frame = parent;
        this.dataPanel = dataPanel;

        traitsTableModel = new TraitsTableModel();
//        TableSorter sorter = new TableSorter(traitsTableModel);
//        predictorsTable = new JTable(sorter);
//        sorter.setTableHeader(predictorsTable.getTableHeader());
        traitsTable = new JTable(traitsTableModel);

        traitsTable.getTableHeader().setReorderingAllowed(false);
        traitsTable.getTableHeader().setResizingAllowed(false);
//        predictorsTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = traitsTable.getColumnModel().getColumn(1);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer(TraitData.TraitType.values());
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        col.setCellRenderer(comboBoxRenderer);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(traitsTable);

        traitsTable.getSelectionModel().setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        traitsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                traitSelectionChanged();
//                dataTableModel.fireTableDataChanged();
            }
        });

        JScrollPane scrollPane1 = new JScrollPane(traitsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane1.setOpaque(false);

        dataTableModel = new DataTableModel();
        TableSorter sorter = new TableSorter(dataTableModel);
        dataTable = new JTable(sorter);
        sorter.setTableHeader(dataTable.getTableHeader());

        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        dataTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        dataTable.getColumnModel().getColumn(0).setPreferredWidth(80);

        col = dataTable.getColumnModel().getColumn(1);
        comboBoxRenderer = new ComboBoxRenderer();
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        col.setCellRenderer(comboBoxRenderer);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

//        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
//            public void valueChanged(ListSelectionEvent evt) {
//                traitSelectionChanged();
//            }
//        });

        JScrollPane scrollPane2 = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane2.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setBorder(BorderFactory.createEmptyBorder());

        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));

        JButton button;

        button = new JButton(addTraitAction);
        PanelUtils.setupComponent(button);
        button.setToolTipText(ADD_TRAITS_TOOLTIP);
        toolBar1.add(button);

        this.importTraitsAction = importTraitsAction;
        button = new JButton(importTraitsAction);
        PanelUtils.setupComponent(button);
        button.setToolTipText(IMPORT_TRAITS_TOOLTIP);
        toolBar1.add(button);

        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));

        button = new JButton(guessTraitsAction);
        PanelUtils.setupComponent(button);
        button.setToolTipText(GUESS_TRAIT_VALUES_TOOLTIP);
        toolBar1.add(button);

        button = new JButton(setValueAction);
        PanelUtils.setupComponent(button);
        button.setToolTipText(SET_TRAIT_VALUES_TOOLTIP);
        toolBar1.add(button);

        // Don't see the need for a clear values button
//        button = new JButton(new ClearTraitAction());
//        PanelUtils.setupComponent(button);
//        button.setToolTipText(CLEAR_TRAIT_VALUES_TOOLTIP);
//        toolBar1.add(button);
        button = new JButton(createTraitPartitionAction);
        PanelUtils.setupComponent(button);
        button.setToolTipText(CREATE_TRAIT_PARTITIONS_TOOLTIP);
        toolBar1.add(button);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addTraitAction);
        actionPanel1.setRemoveAction(removeTraitAction);
        actionPanel1.setAddToolTipText(ADD_TRAITS_TOOLTIP);

        removeTraitAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane1, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);
        panel1.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, 0));

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(toolBar1, BorderLayout.NORTH);
        panel2.add(scrollPane2, BorderLayout.CENTER);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                panel1, panel2);
        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.CENTER);
        add(toolBar1, BorderLayout.NORTH);

    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        updateButtons();

//        int selRow = predictorsTable.getSelectedRow();
//        traitsTableModel.fireTableDataChanged();
//
//        if (selRow < 0) {
//            selRow = 0;
//        }
//        predictorsTable.getSelectionModel().setSelectionInterval(selRow, selRow);


//        if (selectedTrait == null) {
//            predictorsTable.getSelectionModel().setSelectionInterval(0, 0);
//        }

        traitsTableModel.fireTableDataChanged();
        dataTableModel.fireTableDataChanged();

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
//        int selRow = predictorsTable.getSelectedRow();
//        if (selRow >= 0 && options.traitsOptions.selecetedTraits.size() > 0) {
//            selectedTrait = options.traitsOptions.selecetedTraits.get(selRow);
//        }
//        options.datesUnits = unitsCombo.getSelectedIndex();
//        options.datesDirection = directionCombo.getSelectedIndex();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void fireTraitsChanged() {
        if (currentTrait != null) {
//        if (currentTrait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString())) {
//            frame.setupStarBEAST();
//        } else
            if (currentTrait != null && currentTrait.getTraitType() == TraitData.TraitType.DISCRETE) {
                frame.updateDiscreteTraitAnalysis();
            }

//            if (selRow > 0) {
//                predictorsTable.getSelectionModel().setSelectionInterval(selRow-1, selRow-1);
//            } else if (selRow == 0 && options.traitsOptions.traits.size() > 0) { // options.traitsOptions.traits.size() after remove
//                predictorsTable.getSelectionModel().setSelectionInterval(0, 0);
//            }

            traitsTableModel.fireTableDataChanged();
//            options.updatePartitionAllLinks();
            frame.setDirty();
        }
    }

    private void traitSelectionChanged() {
        int selRow = traitsTable.getSelectedRow();
        if (selRow >= 0) {
            currentTrait = options.traits.get(selRow);
//            predictorsTable.getSelectionModel().setSelectionInterval(selRow, selRow);
            removeTraitAction.setEnabled(true);
//        } else {
//            currentTrait = null;
//            removeTraitAction.setEnabled(false);
        }

        if (options.traits.size() <= 0) {
            currentTrait = null;
            removeTraitAction.setEnabled(false);
        }
        dataTableModel.fireTableDataChanged();
//        traitsTableModel.fireTableDataChanged();
    }

    private void updateButtons() { //TODO: better to merge updateButtons() fireTraitsChanged() traitSelectionChanged() into one
        boolean hasData = options.hasData();

        addTraitAction.setEnabled(hasData);
        importTraitsAction.setEnabled(hasData);
        createTraitPartitionAction.setEnabled(hasData && options.traits.size() > 0);
        guessTraitsAction.setEnabled(hasData && options.traits.size() > 0);
        setValueAction.setEnabled(hasData && options.traits.size() > 0);
    }

    public void clearTraitValues(String traitName) {
        options.clearTraitValues(traitName);

        dataTableModel.fireTableDataChanged();
    }

    public void guessTrait() {
        if (options.taxonList == null) { // validation of check empty taxonList
            return;
        }

        if (currentTrait == null) {
            if (!addTrait()) {
                return; // if addTrait() cancel then false
            }
        }
        int result;
        do {
            TraitGuesser currentTraitGuesser = new TraitGuesser(currentTrait);
            if (guessTraitDialog == null) {
                guessTraitDialog = new GuessTraitDialog(frame);
            }
            guessTraitDialog.setDescription("Extract values for trait '" + currentTrait + "' from taxa labels");
            result = guessTraitDialog.showDialog();

            if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                return;
            }

            guessTraitDialog.setupGuesserFromDialog(currentTraitGuesser);

            try {
                int[] selRows = dataTable.getSelectedRows();
                if (selRows.length > 0) {
                    Taxa selectedTaxa = new Taxa();

                    for (int row : selRows) {
                        Taxon taxon = (Taxon) dataTable.getValueAt(row, 0);
                        selectedTaxa.addTaxon(taxon);
                    }
                    currentTraitGuesser.guessTrait(selectedTaxa);
                } else {
                    currentTraitGuesser.guessTrait(options.taxonList);
                }

            } catch (IllegalArgumentException iae) {
                JOptionPane.showMessageDialog(this, iae.getMessage(), "Unable to guess trait value", JOptionPane.ERROR_MESSAGE);
                result = -1;
            }

            dataTableModel.fireTableDataChanged();
        } while (result < 0);
    }

    public void setTraitValue() {
        if (options.taxonList == null) { // validation of check empty taxonList
            return;
        }

        int result;
        do {
            if (traitValueDialog == null) {
                traitValueDialog = new TraitValueDialog(frame);
            }

            int[] selRows = dataTable.getSelectedRows();

            if (selRows.length > 0) {
                traitValueDialog.setDescription("Set values for trait '" + currentTrait + "' for selected taxa");
            } else {
                traitValueDialog.setDescription("Set values for trait '" + currentTrait + "' for all taxa");
            }

            result = traitValueDialog.showDialog();

            if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                return;
            }

//            currentTrait.guessTrait = true; // ?? no use?
            String value = traitValueDialog.getTraitValue();

            try {
                if (selRows.length > 0) {
                    for (int row : selRows) {
                        Taxon taxon = (Taxon) dataTable.getValueAt(row, 0);
                        taxon.setAttribute(currentTrait.getName(), value);
                    }
                } else {
                    for (Taxon taxon : options.taxonList) {
                        taxon.setAttribute(currentTrait.getName(), value);
                    }
                }

            } catch (IllegalArgumentException iae) {
                JOptionPane.showMessageDialog(this, iae.getMessage(), "Unable to guess trait value", JOptionPane.ERROR_MESSAGE);
                result = -1;
            }

            dataTableModel.fireTableDataChanged();
        } while (result < 0);
    }

    public boolean addTrait() {
        return addTrait("Untitled");
    }

    public boolean addTrait(String traitName) {
        return addTrait(null, traitName, false);
    }

    public boolean addTrait(String message, String traitName, boolean isSpeciesTrait) {
        if (createTraitDialog == null) {
            createTraitDialog = new CreateTraitDialog(frame);
        }

        createTraitDialog.setSpeciesTrait(isSpeciesTrait);
        createTraitDialog.setTraitName(traitName);
        createTraitDialog.setMessage(message);

        int result = createTraitDialog.showDialog();
        if (result == JOptionPane.OK_OPTION) {
            frame.tabbedPane.setSelectedComponent(this);

            String name = createTraitDialog.getName();
            TraitData.TraitType type = createTraitDialog.getType();
            TraitData newTrait = new TraitData(options, name, "", type);
            currentTrait = newTrait;

            // The createTraitDialog will have already checked if the
            // user is overwriting an existing trait
            addTrait(newTrait);

            if (createTraitDialog.createTraitPartition()) {
                options.createPartitionForTraits(name, newTrait);
            }

            fireTraitsChanged();
            updateButtons();

        } else if (result == CreateTraitDialog.OK_IMPORT) {
            boolean done = frame.doImportTraits();
            if (done) {
                if (isSpeciesTrait) {
                    // check that we did indeed import a 'species' trait
                    if (!options.traitExists(TraitData.TRAIT_SPECIES)) {
                        JOptionPane.showMessageDialog(this,
                                "The imported trait file didn't contain a trait\n" +
                                        "called '" + TraitData.TRAIT_SPECIES + "', required for *BEAST.\n" +
                                        "Please edit it or select a different file.",
                                "Reserved trait name",
                                JOptionPane.WARNING_MESSAGE);

                        return false;
                    }
                }
                updateButtons();
            }
            return done;
        } else if (result == JOptionPane.CANCEL_OPTION) {
            return false;
        }

        return true;
    }

    public void createTraitPartition() {
        int[] selRows = traitsTable.getSelectedRows();
        java.util.List<TraitData> traits = new ArrayList<TraitData>();
        int discreteCount = 0;
        int continuousCount = 0;
        for (int row : selRows) {
            TraitData trait = options.traits.get(row);
            traits.add(trait);

            if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                discreteCount ++;
            }
            if (trait.getTraitType() == TraitData.TraitType.CONTINUOUS) {
                continuousCount ++;
            }
        }

        boolean success = false;
        if (discreteCount > 0) {
            if (continuousCount > 0)  {
                JOptionPane.showMessageDialog(TraitsPanel.this, "Don't mix discrete and continuous traits when creating partition(s).", "Mixed Trait Types", JOptionPane.ERROR_MESSAGE);
                return;
            }

            // with discrete traits, create a separate partition for each
            for (TraitData trait : traits) {
                java.util.List<TraitData> singleTrait = new ArrayList<TraitData>();
                singleTrait.add(trait);
                if (dataPanel.createFromTraits(singleTrait)) {
                    success = true;
                }
            }
        } else {
            // with
            success = dataPanel.createFromTraits(traits);
        }
        if (success) {
            frame.switchToPanel(BeautiFrame.DATA_PARTITIONS);
        }
    }

    public void addTrait(TraitData newTrait) {
        int selRow = options.addTrait(newTrait);
        traitsTable.getSelectionModel().setSelectionInterval(selRow, selRow);
    }

    private void removeTrait() {
        int selRow = traitsTable.getSelectedRow();
        removeTrait(traitsTable.getValueAt(selRow, 0).toString());
    }

    public void removeTrait(String traitName) {
        if (options.useStarBEAST && traitName.equalsIgnoreCase(TraitData.TRAIT_SPECIES)) {
            JOptionPane.showMessageDialog(this, "The trait named '" + traitName + "' is being used by *BEAST.\nTurn *BEAST off before deleting this trait.", "Trait in use", JOptionPane.ERROR_MESSAGE);
            return;
        }
        TraitData traitData = options.getTrait(traitName);
        if (options.getTraitPartitions(traitData).size() > 0) {
            JOptionPane.showMessageDialog(this, "The trait named '" + traitName + "' is being used in a partition.\nRemove the partition before deleting this trait.", "Trait in use", JOptionPane.ERROR_MESSAGE);
            return;
        }
        options.removeTrait(traitName);

        updateButtons();
        fireTraitsChanged();
        traitSelectionChanged();
    }

    public class ClearTraitAction extends AbstractAction {
        private static final long serialVersionUID = -7281309694753868635L;

        public ClearTraitAction() {
            super("Clear trait values");
        }

        public void actionPerformed(ActionEvent ae) {
            if (currentTrait != null) clearTraitValues(currentTrait.getName()); // Clear trait values
        }
    }

    public class GuessTraitsAction extends AbstractAction {
        private static final long serialVersionUID = 8514706149822252033L;

        public GuessTraitsAction() {
            super("Guess trait values");
        }

        public void actionPerformed(ActionEvent ae) {
            guessTrait();
        }
    }

    public class AddTraitAction extends AbstractAction {

        public AddTraitAction() {
            super("Add trait");
        }

        public void actionPerformed(ActionEvent ae) {
            addTrait();
        }
    }

    AbstractAction removeTraitAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            removeTrait();
        }
    };

    public class SetValueAction extends AbstractAction {

        public SetValueAction() {
            super("Set trait values");
            setToolTipText("Use this button to set the trait values of selected taxa");
        }

        public void actionPerformed(ActionEvent ae) {
            setTraitValue();
        }
    }

    public class CreateTraitPartitionAction extends AbstractAction {
        public CreateTraitPartitionAction() {
            super("Create partition from trait ...");
        }

        public void actionPerformed(ActionEvent ae) {
            createTraitPartition();
        }
    }

    class TraitsTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Trait", "Type"};

        public TraitsTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.traits.size();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.traits.get(row).getName();
                case 1:
                    return options.traits.get(row).getTraitType();
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            switch (col) {
                case 0:
                    String oldName = options.traits.get(row).getName();
                    options.traits.get(row).setName(aValue.toString());
                    Object value;
                    for (Taxon t : options.taxonList) {
                        value = t.getAttribute(oldName);
                        t.setAttribute(aValue.toString(), value);
                        // cannot remvoe attribute in Attributable inteface
                    }
                    fireTraitsChanged();
                    break;
                case 1:
                    options.traits.get(row).setTraitType((TraitData.TraitType) aValue);
                    break;
            }
        }

        public boolean isCellEditable(int row, int col) {
            return !options.useStarBEAST || !options.traits.get(row).getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES.toString());
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
            if (currentTrait == null) return 0;

            return options.taxonList.getTaxonCount();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.taxonList.getTaxon(row);
                case 1:
                    Object value = null;
                    if (currentTrait != null) {
                        value = options.taxonList.getTaxon(row).getAttribute(currentTrait.getName());
                    }
                    if (value != null) {
                        return value;
                    } else {
                        return "";
                    }
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            if (col == 1) {
//                Location location = options.taxonList.getTaxon(row).getLocation();
//                if (location != null) {
//                    options.taxonList.getTaxon(row).setLocation(location);
//                }
                if (currentTrait != null) {
                    options.taxonList.getTaxon(row).setAttribute(currentTrait.getName(), aValue);
                }

            }
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 1) {
//                Object t = options.taxonList.getTaxon(row).getAttribute(currentTrait.getName());
//                return (t != null);
                return true;
            } else {
                return false;
            }
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

