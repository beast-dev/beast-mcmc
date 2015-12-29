/*
 * SpeciesSetPanel.java
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

package dr.app.beauti.taxonsetspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitData;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * It is specified to *BEAST and used to replace Taxon Sets panel,
 * because *BEAST calibration is only allowed in species tree
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class SpeciesSetPanel extends TaxonSetPanel {

    private final String[] columnToolTips = {"The set of species defined for the calibration",
            "Enforce the selected species set to be monophyletic on the specified tree"};
    protected final String TAXA = "Species";
    protected final String TAXON = "Species set";

    protected SpeciesSetsTableModel speciesSetsTableModel = new SpeciesSetsTableModel();

    public SpeciesSetPanel(BeautiFrame parent) {
        super.frame = parent;

        setText(true);

        initTaxonSetsTable(speciesSetsTableModel, columnToolTips);

        initTableColumn();

        initPanel(addSpeciesSetAction, removeSpeciesSetAction);
    }

    protected void initTableColumn() {
        final TableColumnModel tableColumnModel = taxonSetsTable.getColumnModel();
        TableColumn tableColumn = tableColumnModel.getColumn(0);
        tableColumn.setCellRenderer(new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        tableColumn.setMinWidth(20);

        tableColumn = tableColumnModel.getColumn(1);
        tableColumn.setPreferredWidth(10);
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        resetPanel();

        if (options.speciesSets == null) {
            addSpeciesSetAction.setEnabled(false);
            removeSpeciesSetAction.setEnabled(false);
        } else if (options.starBEASTOptions.getEmptySpeciesIndex() < 0) {
            addSpeciesSetAction.setEnabled(true);
        }
    }
    
    protected void taxonSetsTableSelectionChanged() {
        treeModelsChanged();

        int[] rows = taxonSetsTable.getSelectedRows();
        if (rows.length == 0) {
            removeSpeciesSetAction.setEnabled(false);
        } else if (rows.length == 1) {
            currentTaxonSet = options.speciesSets.get(rows[0]);
            setCurrentTaxonSet(currentTaxonSet);
            removeSpeciesSetAction.setEnabled(true);
        } else {
            setCurrentTaxonSet(null);
            removeSpeciesSetAction.setEnabled(true);
        }
    }
    
    protected void taxonSetChanged() {
        currentTaxonSet.removeAllTaxa();
        for (Taxon anIncludedTaxa : includedTaxa) {
            currentTaxonSet.addTaxon(anIncludedTaxa);
        }

        setupTaxonSetsComboBoxes();

        if (options.speciesSetsMono.get(currentTaxonSet) != null &&
                options.speciesSetsMono.get(currentTaxonSet) &&
                !checkCompatibility(currentTaxonSet)) {
            options.speciesSetsMono.put(currentTaxonSet, Boolean.FALSE);
        }

        frame.setDirty();
    }

    protected void resetPanel() {
        if (!options.hasData() || options.speciesSets == null || options.speciesSets.size() < 1) {
            setCurrentTaxonSet(null);
        }
    }

    protected void setCurrentTaxonSet(Taxa taxonSet) {

        currentTaxonSet = taxonSet;

        includedTaxa.clear();
        excludedTaxa.clear();

        if (currentTaxonSet != null) {
            for (int i = 0; i < taxonSet.getTaxonCount(); i++) {
                includedTaxa.add(taxonSet.getTaxon(i));
            }
            Collections.sort(includedTaxa);

            Set<String> allSpecies = TraitData.getStatesListOfTrait(options.taxonList, TraitData.TRAIT_SPECIES);

            for (String sp : allSpecies) {
                excludedTaxa.add(new Taxon(sp));
            }
            excludedTaxa.removeAll(includedTaxa);
            Collections.sort(excludedTaxa);
        }

        setTaxonSetTitle();

        setupTaxonSetsComboBoxes();

        includedTaxaTableModel.fireTableDataChanged();
        excludedTaxaTableModel.fireTableDataChanged();
    }

    protected void setupTaxonSetsComboBox(JComboBox comboBox, List<Taxon> availableTaxa) {
        comboBox.removeAllItems();

        comboBox.addItem(TAXON.toLowerCase() + "...");
        for (Taxa taxa : options.speciesSets) {
            if (taxa != currentTaxonSet) {
                if (isCompatible(taxa, availableTaxa)) {
                    comboBox.addItem(taxa);
                }
            }
        }
    }
    
    protected boolean checkCompatibility(Taxa taxa) {
        for (Taxa taxa2 : options.speciesSets) {
            if (taxa2 != taxa && options.speciesSetsMono.get(taxa2)) {
                if (taxa.containsAny(taxa2) && !taxa.containsAll(taxa2) && !taxa2.containsAll(taxa)) {
                    JOptionPane.showMessageDialog(frame,
                            "You cannot enforce monophyly on this " + TAXON.toLowerCase() + " \n" +
                                    "because it is not compatible with another " + TAXON.toLowerCase() + ",\n" +
                                    taxa2.getId() + ", for which monophyly is\n" + "enforced.",
                            "Warning",
                            JOptionPane.WARNING_MESSAGE);
                    return false;
                }
            }
        }
        return true;
    }
    
    protected void treeModelsChanged() { }

    Action addSpeciesSetAction = new AbstractAction("+") {

        public void actionPerformed(ActionEvent ae) {
            taxonSetCount++;

            String newSpeciesSetName = "untitled" + taxonSetCount;
            Taxa newSpeciesSet = new Taxa(newSpeciesSetName); // cannot use currentTaxonSet

            options.speciesSets.add(newSpeciesSet);
            Collections.sort(options.speciesSets);

            options.speciesSetsMono.put(newSpeciesSet, Boolean.FALSE);

            setCurrentTaxonSet(newSpeciesSet);

            taxonSetChanged();

            speciesSetsTableModel.fireTableDataChanged();

            int sel = options.getSpeciesIndex(newSpeciesSetName);
            if (sel < 0) {
                taxonSetsTable.setRowSelectionInterval(0, 0);
            } else {
                taxonSetsTable.setRowSelectionInterval(sel, sel);
            }
        }
    };

    Action removeSpeciesSetAction = new AbstractAction("-") {

        public void actionPerformed(ActionEvent ae) {
            int row = taxonSetsTable.getSelectedRow();
            if (row != -1) {
                Taxa taxa = options.speciesSets.remove(row);
                options.speciesSetsMono.remove(taxa);
            }
            taxonSetChanged();

            speciesSetsTableModel.fireTableDataChanged();

            if (row >= options.speciesSets.size()) {
                row = options.speciesSets.size() - 1;
            }
            if (row >= 0) {
                taxonSetsTable.setRowSelectionInterval(row, row);
            } else {
                setCurrentTaxonSet(null);
            }
        }
    };


    /**
     * The table on the left side of panel
     */
    protected class SpeciesSetsTableModel extends TaxonSetPanel.TaxonSetsTableModel {

        String[] columnNames = {"Species Sets", "Monophyletic?"};

        public SpeciesSetsTableModel() {
            super();
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.speciesSets.size();
        }

        public Object getValueAt(int rowIndex, int columnIndex) {
            Taxa taxonSet = options.speciesSets.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    return taxonSet.getId();
                case 1:
                    return options.speciesSetsMono.get(taxonSet);

                default:
                    throw new IllegalArgumentException("unknown column, " + columnIndex);
            }
        }

        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            Taxa taxonSet = options.speciesSets.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    taxonSet.setId(aValue.toString());
                    options.renameTMRCAStatistic(taxonSet);
                    setTaxonSetTitle();
                    break;

                case 1:
                    if ((Boolean) aValue) {
                        Taxa taxa = options.speciesSets.get(rowIndex);
                        if (checkCompatibility(taxa)) {
                            options.speciesSetsMono.put(taxonSet, (Boolean) aValue);
                        }
                    } else {
                        options.speciesSetsMono.put(taxonSet, (Boolean) aValue);
                    }
                    break;

                default:
                    throw new IllegalArgumentException("unknown column, " + columnIndex);
            }
        }
    }
}
