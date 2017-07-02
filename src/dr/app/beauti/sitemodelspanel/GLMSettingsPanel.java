/*
 * PriorSettingsPanel.java
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

package dr.app.beauti.sitemodelspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.util.BEAUTiImporter;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.app.util.OSType;
import jam.panels.ActionPanel;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * @author Andrew Rambaut
 */
public class GLMSettingsPanel extends JPanel {

    private static final int MINIMUM_TABLE_WIDTH = 140;

    private static final String IMPORT_PREDICTORS_TOOLTIP = "<html>Import one or more predictor matrices for these taxa from a tab-delimited<br>" +
            "file. Taxa should be in the first column and the trait names<br>" +
            "in the first row</html>";

    private final BeautiFrame frame;
    private JDialog dialog;
    private BeautiOptions options;

    public final JTable predictorsTable;
    private final PredictorsTableModel predictorsTableModel;

    public GLMSettingsPanel(BeautiFrame frame) {
        this.frame = frame;

        predictorsTableModel = new PredictorsTableModel();
        TableSorter sorter = new TableSorter(predictorsTableModel);
        predictorsTable = new JTable(sorter);
        sorter.setTableHeader(predictorsTable.getTableHeader());

        predictorsTable.getTableHeader().setReorderingAllowed(false);
        predictorsTable.getTableHeader().setResizingAllowed(false);
//        predictorsTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = predictorsTable.getColumnModel().getColumn(1);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer(TraitData.TraitType.values());
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        col.setCellRenderer(comboBoxRenderer);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(predictorsTable);

        predictorsTable.getSelectionModel().setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        predictorsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
//                traitSelectionChanged();
//                dataTableModel.fireTableDataChanged();
            }
        });

        JScrollPane scrollPane1 = new JScrollPane(predictorsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane1.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);

        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));

        JButton button;

        button = new JButton(importPredictorsAction);
        PanelUtils.setupComponent(button);
        button.setToolTipText(IMPORT_PREDICTORS_TOOLTIP);
        toolBar1.add(button);

        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(importPredictorsAction);
        actionPanel1.setRemoveAction(removePredictorAction);
        actionPanel1.setAddToolTipText(IMPORT_PREDICTORS_TOOLTIP);

        importPredictorsAction.setEnabled(true);
        removePredictorAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane1, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);
//        panel1.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, 0));

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(panel1, BorderLayout.CENTER);
        add(toolBar1, BorderLayout.NORTH);
    }

    /**
     * Set the trait to be controlled
     *                                                                                                       q
     * @param trait
     */
    public void setTrait(final String trait) {
//        this.parameter = parameter;
//
//        priorCombo = new JComboBox();
//        for (PriorType priorType : PriorType.getPriorTypes(parameter)) {
//            priorCombo.addItem(priorType);
//        }
//
//        if (parameter.priorType != null) {
//            priorCombo.setSelectedItem(parameter.priorType);
//        }
//
//        setupComponents(); // setArguments here
//
//
//        priorCombo.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent e) {
//                setupComponents();
//                dialog.pack();
//                dialog.repaint();
//            }
//        });
//
//        for (PriorOptionsPanel optionsPanel : optionsPanels.values()) {
//            optionsPanel.removeAllListeners();
//            optionsPanel.addListener(new PriorOptionsPanel.Listener() {
//                public void optionsPanelChanged() {
//                    setupChart();
//                    dialog.pack();
//                    dialog.repaint();
//                }
//            });
//        }
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;
    }

    public void getOptions(BeautiOptions options) {
    }

    private void setupComponents() {
        removeAll();

        OptionsPanel optionsPanel = new OptionsPanel(12, (OSType.isMac() ? 6 : 24));

        JPanel panel1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel1.add(optionsPanel);


        repaint();
    }

    public void setDialog(final JDialog dialog) {
        this.dialog = dialog;
    }

    public void firePredictorsChanged() {
//        if (currentTrait != null) {
////        if (currentTrait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString())) {
////            frame.setupStarBEAST();
////        } else
//            if (currentTrait != null && currentTrait.getTraitType() == TraitData.TraitType.DISCRETE) {
//                frame.updateDiscreteTraitAnalysis();
//            }
//
////            if (selRow > 0) {
////                predictorsTable.getSelectionModel().setSelectionInterval(selRow-1, selRow-1);
////            } else if (selRow == 0 && options.traitsOptions.traits.size() > 0) { // options.traitsOptions.traits.size() after remove
////                predictorsTable.getSelectionModel().setSelectionInterval(0, 0);
////            }
//
//            traitsTableModel.fireTableDataChanged();
//            options.updatePartitionAllLinks();
//            frame.setDirty();
//        }
        predictorsTableModel.fireTableDataChanged();
    }

    private void updateButtons() {
//        boolean hasData = options.hasData();
//
//        addTraitAction.setEnabled(hasData);
//        importTraitsAction.setEnabled(hasData);
    }

//    public boolean addTrait() {
//        return addTrait("Untitled");
//    }
//
//    public boolean addTrait(String traitName) {
//        return addTrait(null, traitName, false);
//    }
//
//    public boolean addTrait(String message, String traitName, boolean isSpeciesTrait) {
//        if (createTraitDialog == null) {
//            createTraitDialog = new CreateTraitDialog(frame);
//        }
//
//        createTraitDialog.setSpeciesTrait(isSpeciesTrait);
//        createTraitDialog.setTraitName(traitName);
//        createTraitDialog.setMessage(message);
//
//        int result = createTraitDialog.showDialog();
//        if (result == JOptionPane.OK_OPTION) {
//            frame.tabbedPane.setSelectedComponent(this);
//
//            String name = createTraitDialog.getName();
//            TraitData.TraitType type = createTraitDialog.getType();
//            TraitData newTrait = new TraitData(options, name, "", type);
//            currentTrait = newTrait;
//
//            // The createTraitDialog will have already checked if the
//            // user is overwriting an existing trait
//            addTrait(newTrait);
//
//            if (createTraitDialog.createTraitPartition()) {
//                options.createPartitionForTraits(name, newTrait);
//            }
//
//            fireTraitsChanged();
//            updateButtons();
//
//        } else if (result == CreateTraitDialog.OK_IMPORT) {
//            boolean done = frame.doImportTraits();
//            if (done) {
//                if (isSpeciesTrait) {
//                    // check that we did indeed import a 'species' trait
//                    if (!options.traitExists(TraitData.TRAIT_SPECIES)) {
//                        JOptionPane.showMessageDialog(this,
//                                "The imported trait file didn't contain a trait\n" +
//                                        "called '" + TraitData.TRAIT_SPECIES + "', required for *BEAST.\n" +
//                                        "Please edit it or select a different file.",
//                                "Reserved trait name",
//                                JOptionPane.WARNING_MESSAGE);
//
//                        return false;
//                    }
//                }
//                updateButtons();
//            }
//            return done;
//        } else if (result == JOptionPane.CANCEL_OPTION) {
//            return false;
//        }
//
//        return true;
//    }
//
//    public void createTraitPartition() {
//        int[] selRows = predictorsTable.getSelectedRows();
//        java.util.List<TraitData> traits = new ArrayList<TraitData>();
//        int discreteCount = 0;
//        int continuousCount = 0;
//        for (int row : selRows) {
//            TraitData trait = options.traits.get(row);
//            traits.add(trait);
//
//            if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
//                discreteCount ++;
//            }
//            if (trait.getTraitType() == TraitData.TraitType.CONTINUOUS) {
//                continuousCount ++;
//            }
//        }
//
//        boolean success = false;
//        if (discreteCount > 0) {
//            if (continuousCount > 0)  {
//                JOptionPane.showMessageDialog(TraitsPanel.this, "Don't mix discrete and continuous traits when creating partition(s).", "Mixed Trait Types", JOptionPane.ERROR_MESSAGE);
//                return;
//            }
//
//            // with discrete traits, create a separate partition for each
//            for (TraitData trait : traits) {
//                java.util.List<TraitData> singleTrait = new ArrayList<TraitData>();
//                singleTrait.add(trait);
//                if (dataPanel.createFromTraits(singleTrait)) {
//                    success = true;
//                }
//            }
//        } else {
//            // with
//            success = dataPanel.createFromTraits(traits);
//        }
//        if (success) {
//            frame.switchToPanel(BeautiFrame.DATA_PARTITIONS);
//        }
//    }
//
//    public void addTrait(TraitData newTrait) {
//        int selRow = options.addTrait(newTrait);
//        predictorsTable.getSelectionModel().setSelectionInterval(selRow, selRow);
//    }

    private boolean importPredictors() {
        File[] files = frame.selectImportFiles("Import Predictors File...", false, new FileNameExtensionFilter[] {
                new FileNameExtensionFilter("Tab-delimited text files", "txt", "tab", "dat") });

        if (files != null && files.length != 0) {
            try {
                BEAUTiImporter beautiImporter = new BEAUTiImporter(frame, options);
                beautiImporter.importPredictors(files[0]);
            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
                return false;
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe.getMessage(),
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
                return false;
            } catch (Exception ex) {
                ex.printStackTrace(System.err);
                JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                        "Error reading file",
                        JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace();
                return false;
            }
        } else {
            return false;
        }


        firePredictorsChanged();
//        setAllOptions();

        return true;
    }

    private void removePredictor() {
        int selRow = predictorsTable.getSelectedRow();
        removePredictor(predictorsTable.getValueAt(selRow, 0).toString());
    }

    public void removePredictor(String traitName) {
//        TraitData traitData = options.getTrait(traitName);
//        if (options.getTraitPartitions(traitData).size() > 0) {
//            JOptionPane.showMessageDialog(this, "The trait named '" + traitName + "' is being used in a partition.\nRemove the partition before deleting this trait.", "Trait in use", JOptionPane.ERROR_MESSAGE);
//            return;
//        }
//        options.removeTrait(traitName);
//
//        updateButtons();
//        fireTraitsChanged();
//        traitSelectionChanged();
    }

    AbstractAction importPredictorsAction = new AbstractAction("Import Predictors...") {
        public void actionPerformed(ActionEvent ae) {
            importPredictors();
        }
    };

    AbstractAction removePredictorAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            removePredictor();
        }
    };

    class PredictorsTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Trait", "Type"};

        public PredictorsTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
//            if (options == null) return 0;
//            return options.traits.size();
            return 0;
        }

        public Object getValueAt(int row, int col) {
//            switch (col) {
//                case 0:
//                    return options.traits.get(row).getName();
//                case 1:
//                    return options.traits.get(row).getTraitType();
//            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
//            switch (col) {
//                case 0:
//                    String oldName = options.traits.get(row).getName();
//                    options.traits.get(row).setName(aValue.toString());
//                    Object value;
//                    for (Taxon t : options.taxonList) {
//                        value = t.getAttribute(oldName);
//                        t.setAttribute(aValue.toString(), value);
//                        // cannot remvoe attribute in Attributable inteface
//                    }
//                    fireTraitsChanged();
//                    break;
//                case 1:
//                    options.traits.get(row).setTraitType((TraitData.TraitType) aValue);
//                    break;
//            }
        }

        public boolean isCellEditable(int row, int col) {
//            return !options.useStarBEAST || !options.traits.get(row).getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES.toString());
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