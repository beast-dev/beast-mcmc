/*
 * SpeciesSetPanel.java
 *
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

package dr.app.beauti.taxonsetspanel;

import dr.app.beauti.BeautiApp;
import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitData;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * It is specified to *BEAST and used to replace Taxon Sets panel,
 * because *BEAST calibration is only allowed in species tree
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class SpeciesSetPanel extends BeautiPanel implements Exportable {

    private final String SP_SET_DEFAULT = "species set...";
    private final String[] columnToolTips = {"The set of species defined for the calibration",
            "Enforce the selected species set to be monophyletic on the specified tree"};

    BeautiFrame frame = null;
    BeautiOptions options = null;

    private JTable speciesSetsTable = null;
    private final TableColumnModel tableColumnModel;
    private SpeciesSetsTableModel speciesSetsTableModel = null;

    private JPanel speciesSetEditingPanel = null;

    private List<String> currentSpeciesSet = null; // 1st element is name, 2nd is monophyletic

    // these 2 lists below only contain species
    private final List<String> includedSpecies = new ArrayList<String>();
    private final List<String> excludedSpecies = new ArrayList<String>();

    private JTable excludedSpeciesTable = null;
    private SpeciesTableModel excludedSpeciesTableModel = null;
    private JComboBox excludedSpeciesSetsComboBox = null;
    private boolean excludedSelectionChanging = false;

    private JTable includedSpeciesTable = null;
    private SpeciesTableModel includedSpeciesTableModel = null;
    private JComboBox includedSpeciesSetsComboBox = null;
    private boolean includedSelectionChanging = false;

    private static int speciesSetCount = 0;

    public SpeciesSetPanel(BeautiFrame parent) {

        this.frame = parent;

        Icon includeIcon = null, excludeIcon = null;
        try {
            includeIcon = new ImageIcon(IconUtils.getImage(BeautiApp.class, "images/include.png"));
            excludeIcon = new ImageIcon(IconUtils.getImage(BeautiApp.class, "images/exclude.png"));
        } catch (Exception e) {
            // do nothing
        }

        // species Sets
        speciesSetsTableModel = new SpeciesSetsTableModel();
        speciesSetsTable = new JTable(speciesSetsTableModel) {
            //Implement table header tool tips.
            protected JTableHeader createDefaultTableHeader() {
                return new JTableHeader(columnModel) {
                    public String getToolTipText(MouseEvent e) {
                        Point p = e.getPoint();
                        int index = columnModel.getColumnIndexAtX(p.x);
                        int realIndex = columnModel.getColumn(index).getModelIndex();
                        return columnToolTips[realIndex];
                    }
                };
            }
        };

        speciesSetsTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        tableColumnModel = speciesSetsTable.getColumnModel();
        TableColumn tableColumn = tableColumnModel.getColumn(0);
        tableColumn.setCellRenderer(new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        tableColumn.setMinWidth(20);

        tableColumn = tableColumnModel.getColumn(1);
        tableColumn.setPreferredWidth(10);

        speciesSetsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                speciesSetsTableSelectionChanged();
            }
        });

        JScrollPane scrollPane1 = new JScrollPane(speciesSetsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addSpeciesSetAction);
        actionPanel1.setRemoveAction(removeSpeciesSetAction);

        addSpeciesSetAction.setEnabled(false);
        removeSpeciesSetAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.add(actionPanel1);

        // Excluded Species List
        excludedSpeciesTableModel = new SpeciesTableModel(false);
        excludedSpeciesTable = new JTable(excludedSpeciesTableModel);

        excludedSpeciesTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        excludedSpeciesTable.getColumnModel().getColumn(0).setMinWidth(20);

        excludedSpeciesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                excludedSpeciesTableSelectionChanged();
            }
        });

        JScrollPane scrollPane2 = new JScrollPane(excludedSpeciesTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        Box panel1 = new Box(BoxLayout.X_AXIS);
        panel1.add(new JLabel("Select: "));
        panel1.setOpaque(false);
        excludedSpeciesSetsComboBox = new JComboBox(new String[]{SP_SET_DEFAULT});
        excludedSpeciesSetsComboBox.setOpaque(false);
        panel1.add(excludedSpeciesSetsComboBox);

        JPanel buttonPanel = createAddRemoveButtonPanel(includeSpeciesAction, includeIcon,
                "Include selected species in the species set", excludeSpeciesAction, excludeIcon,
                "Exclude selected species from the species set", BoxLayout.Y_AXIS);

        // Included Species List
        includedSpeciesTableModel = new SpeciesTableModel(true);
        includedSpeciesTable = new JTable(includedSpeciesTableModel);

        includedSpeciesTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        includedSpeciesTable.getColumnModel().getColumn(0).setMinWidth(20);

        includedSpeciesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                includedSpeciesTableSelectionChanged();
            }
        });
        includedSpeciesTable.doLayout();

        JScrollPane scrollPane3 = new JScrollPane(includedSpeciesTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        Box panel2 = new Box(BoxLayout.X_AXIS);
        panel2.add(new JLabel("Select: "));
        panel2.setOpaque(false);
        includedSpeciesSetsComboBox = new JComboBox(new String[]{SP_SET_DEFAULT});
        includedSpeciesSetsComboBox.setOpaque(false);
        panel2.add(includedSpeciesSetsComboBox);

        speciesSetEditingPanel = new JPanel();
        speciesSetEditingPanel.setBorder(BorderFactory.createTitledBorder(""));
        speciesSetEditingPanel.setOpaque(false);
        speciesSetEditingPanel.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0.5;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(12, 12, 4, 0);
        speciesSetEditingPanel.add(scrollPane2, c);

        c.gridx = 0;
        c.gridy = 1;
        c.weightx = 0.5;
        c.weighty = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0, 12, 12, 0);
        speciesSetEditingPanel.add(panel1, c);

        c.gridx = 1;
        c.gridy = 0;
        c.weightx = 0;
        c.weighty = 1;
        c.gridheight = 2;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(12, 2, 12, 4);
        speciesSetEditingPanel.add(buttonPanel, c);

        c.gridx = 2;
        c.gridy = 0;
        c.weightx = 0.5;
        c.weighty = 1;
        c.gridheight = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(12, 0, 4, 12);
        speciesSetEditingPanel.add(scrollPane3, c);

        c.gridx = 2;
        c.gridy = 1;
        c.weightx = 0.5;
        c.weighty = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0, 0, 12, 12);
        speciesSetEditingPanel.add(panel2, c);

        JPanel panel3 = new JPanel();
        panel3.setOpaque(false);
        panel3.setLayout(new GridBagLayout());
        c = new GridBagConstraints();

        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0.5;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0, 0, 2, 12);
        panel3.add(scrollPane1, c);

        c.gridx = 0;
        c.gridy = 1;
        c.weightx = 0;
        c.weighty = 0;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.WEST;
        c.insets = new Insets(2, 0, 0, 12);
        panel3.add(actionPanel1, c);

        c.gridx = 1;
        c.gridy = 0;
        c.weightx = 0.5;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0, 0, 0, 0);
        panel3.add(speciesSetEditingPanel, c);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(panel3, BorderLayout.CENTER);

        includedSpeciesTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    includeSelectedSpecies();
                }
            }
        });
        excludedSpeciesTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    excludeSelectedSpecies();
                }
            }
        });

        includedSpeciesTable.addFocusListener(new FocusAdapter() {
            public void focusGained(FocusEvent focusEvent) {
                excludedSpeciesTable.clearSelection();
            }
        });
        excludedSpeciesTable.addFocusListener(new FocusAdapter() {
            public void focusGained(FocusEvent focusEvent) {
                includedSpeciesTable.clearSelection();
            }
        });

        includedSpeciesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (!includedSelectionChanging && includedSpeciesSetsComboBox.getSelectedIndex() != 0) {
                    includedSpeciesSetsComboBox.setSelectedIndex(0);
                }
            }
        });
        includedSpeciesSetsComboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                includedSelectionChanging = true;
                includedSpeciesTable.clearSelection();
                if (includedSpeciesSetsComboBox.getSelectedIndex() > 0) {
                    Taxa taxa = (Taxa) includedSpeciesSetsComboBox.getSelectedItem();
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        int index = includedSpecies.indexOf(taxon);
                        includedSpeciesTable.getSelectionModel().addSelectionInterval(index, index);

                    }
                }
                includedSelectionChanging = false;
            }
        });

        excludedSpeciesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (!excludedSelectionChanging && excludedSpeciesSetsComboBox.getSelectedIndex() != 0) {
                    excludedSpeciesSetsComboBox.setSelectedIndex(0);
                }
            }
        });
        excludedSpeciesSetsComboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                excludedSelectionChanging = true;
                excludedSpeciesTable.clearSelection();
                if (excludedSpeciesSetsComboBox.getSelectedIndex() > 0) {
                    Taxa taxa = (Taxa) excludedSpeciesSetsComboBox.getSelectedItem();
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        int index = excludedSpecies.indexOf(taxon);
                        excludedSpeciesTable.getSelectionModel().addSelectionInterval(index, index);

                    }
                }
                excludedSelectionChanging = false;
            }
        });

        speciesSetsTable.doLayout();
        includedSpeciesTable.doLayout();
        excludedSpeciesTable.doLayout();
    }

    private void speciesSetChanged() {
        for (int i = 2; i < currentSpeciesSet.size(); i++) {
            currentSpeciesSet.remove(i);
        }

        for (String anIncludedSp : includedSpecies) {
            currentSpeciesSet.add(anIncludedSp);
        }

        setupSpeciesSetsComboBoxes();

        if (currentSpeciesSet != null && Boolean.parseBoolean(currentSpeciesSet.get(1)) &&
                !checkCompatibility(currentSpeciesSet)) {
            currentSpeciesSet.set(1, Boolean.FALSE.toString()); // monophyletic false
        }

        frame.setDirty();
    }

    private void resetPanel() {
        if (!options.hasData() || options.speciesSets == null || options.speciesSets.size() < 1) {
            setCurrentSpeciesSet(null);
        }
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        resetPanel();

        if (options.speciesSets == null) {
            addSpeciesSetAction.setEnabled(false);
            removeSpeciesSetAction.setEnabled(false);
        } else {
            addSpeciesSetAction.setEnabled(true);
        }

        speciesSetsTableSelectionChanged();
        speciesSetsTableModel.fireTableDataChanged();

        validateSpeciesSets();
    }

    private void validateSpeciesSets() {
//        if (speciesSetsTable.getRowCount() > 0) {
//            for (Taxa taxonSet : options.speciesSets) {
//                if (taxonSet.getTaxonCount() < 1) {
//                    JOptionPane.showMessageDialog(this, "Species set " + taxonSet.getId() + " is empty, "
//                            + "\nplease go back to Species Sets panel to select included species.",
//                                "Empty species set error", JOptionPane.ERROR_MESSAGE);
//                }
//            }
//        }
    }

    public void getOptions(BeautiOptions options) {    }

    public JComponent getExportableComponent() {
        return speciesSetsTable;
    }

    private void speciesSetsTableSelectionChanged() {

        int[] rows = speciesSetsTable.getSelectedRows();
        if (rows.length == 0) {
            removeSpeciesSetAction.setEnabled(false);
        } else if (rows.length == 1) {
            currentSpeciesSet = options.speciesSets.get(rows[0]);
            setCurrentSpeciesSet(currentSpeciesSet);
            removeSpeciesSetAction.setEnabled(true);
        } else {
            setCurrentSpeciesSet(null);
            removeSpeciesSetAction.setEnabled(true);
        }
    }

    Action addSpeciesSetAction = new AbstractAction("+") {

        public void actionPerformed(ActionEvent ae) {
            speciesSetCount++;
            // initialize currentSpeciesSet with 1st element is name, 2nd is monophyletic
            currentSpeciesSet = new ArrayList<String>();
            currentSpeciesSet.add("untitled" + speciesSetCount); // name
            currentSpeciesSet.add(Boolean.FALSE.toString()); // monophyletic false

            options.speciesSets.add(currentSpeciesSet);
//            Collections.sort(options.speciesSets);

            speciesSetsTableModel.fireTableDataChanged();

            int sel = options.speciesSets.indexOf(currentSpeciesSet);
            speciesSetsTable.setRowSelectionInterval(sel, sel);

            speciesSetChanged();
        }
    };

    Action removeSpeciesSetAction = new AbstractAction("-") {

        public void actionPerformed(ActionEvent ae) {
            int row = speciesSetsTable.getSelectedRow();
            if (row != -1) {
                Taxa taxa = options.taxonSets.remove(row);
                options.taxonSetsMono.remove(taxa);
                options.taxonSetsIncludeStem.remove(taxa);
            }
            speciesSetChanged();

            speciesSetsTableModel.fireTableDataChanged();

            if (row >= options.taxonSets.size()) {
                row = options.taxonSets.size() - 1;
            }
            if (row >= 0) {
                speciesSetsTable.setRowSelectionInterval(row, row);
            } else {
                setCurrentSpeciesSet(null);
            }
        }
    };

    private void setCurrentSpeciesSet(List<String> speciesSet) {

        this.currentSpeciesSet = speciesSet;

        includedSpecies.clear();
        excludedSpecies.clear();

        if (currentSpeciesSet != null) {
            includedSpecies.addAll(speciesSet.subList(2, speciesSet.size()));
            Collections.sort(includedSpecies);

            Set<String> allSpecies = TraitData.getStatesListOfTrait(options.taxonList, TraitData.TRAIT_SPECIES);

            excludedSpecies.addAll(allSpecies);
            excludedSpecies.removeAll(includedSpecies);
            Collections.sort(excludedSpecies);
        }

        setSpeciesSetTitle();

        setupSpeciesSetsComboBoxes();

        includedSpeciesTableModel.fireTableDataChanged();
        excludedSpeciesTableModel.fireTableDataChanged();
    }

    private void setSpeciesSetTitle() {

        if (currentSpeciesSet == null) {
            speciesSetEditingPanel.setBorder(BorderFactory.createTitledBorder(""));
            speciesSetEditingPanel.setEnabled(false);
        } else {
            speciesSetEditingPanel.setEnabled(true);
            speciesSetEditingPanel.setBorder(BorderFactory.createTitledBorder("Species Set: " + currentSpeciesSet.get(0)));
        }
    }


    private void setupSpeciesSetsComboBoxes() {
        setupSpeciesSetsComboBox(excludedSpeciesSetsComboBox, excludedSpecies);
        excludedSpeciesSetsComboBox.setSelectedIndex(0);
        setupSpeciesSetsComboBox(includedSpeciesSetsComboBox, includedSpecies);
        includedSpeciesSetsComboBox.setSelectedIndex(0);
    }

    private void setupSpeciesSetsComboBox(JComboBox comboBox, List<String> availableTaxa) {
        comboBox.removeAllItems();

        comboBox.addItem(SP_SET_DEFAULT);
        for (List<String> spSet : options.speciesSets) {
            if (spSet.containsAll(currentSpeciesSet) &&
                    currentSpeciesSet.containsAll(spSet.subList(2, spSet.size()))) { // this means spSet = currentSpeciesSet
                if (isCompatible(spSet, availableTaxa)) {
                    comboBox.addItem(spSet);
                }
            }
        }
    }

    /**
     * Returns true if speciesSet are all found in availableSpecies
     */
    private boolean isCompatible(List<String> speciesSet, List<String> availableSpecies) {
        // 1st element is name, 2nd is monophyletic
        for (int i = 2; i < speciesSet.size(); i++) {
            String sp = speciesSet.get(i);
            if (!availableSpecies.contains(sp)) {
                return false;
            }
        }
        return true;
    }

    private boolean checkCompatibility(List<String> speciesSet) {
//        for (List<String> speciesSet2 : options.speciesSets) {
//            if (speciesSet2.get(0).equalsIgnoreCase(speciesSet.get(0)) && Boolean.parseBoolean(speciesSet2.get(1))) {
//                if (containsAny(speciesSet, speciesSet2) && !containsAll(speciesSet, speciesSet2)
//                        && !containsAll(speciesSet2, speciesSet)) {
//                    JOptionPane.showMessageDialog(frame,
//                            "You cannot enforce monophyly on this species set \n" +
//                                    "because it is not compatible with another species \n" +
//                                    "set, " + speciesSet2.get(0) + ", for which monophyly is\n" +
//                                    "enforced.",
//                            "Warning",
//                            JOptionPane.WARNING_MESSAGE);
//                    return false;
//                }
//            }
//        }
        return true;
    }

    /**
     * The table on the left side of panel
     */
    class SpeciesSetsTableModel extends AbstractTableModel {

        String[] columnNames = {"Species Sets", "Monophyletic?"};
        public SpeciesSetsTableModel() {
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
            List<String> speciesSet = options.speciesSets.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    return speciesSet.get(0); // name
                case 1:
                    return speciesSet.get(1); // monophyletic
                default:
                    throw new IllegalArgumentException("unknown column, " + columnIndex);
            }
        }

        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            List<String> speciesSet = options.speciesSets.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    speciesSet.set(0, aValue.toString()); // name
                    setSpeciesSetTitle();
                    break;

                case 1:
                    if ((Boolean) aValue) {
                        if (checkCompatibility(speciesSet)) {
                            speciesSet.set(1, aValue.toString());
                        }
                    } else {
                        speciesSet.set(1, aValue.toString());
                    }
                    break;
                default:
                    throw new IllegalArgumentException("unknown column, " + columnIndex);
            }
        }

        public boolean isCellEditable(int row, int col) {
            return true;
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }
    }

    private JPanel createAddRemoveButtonPanel(Action addAction, Icon addIcon, String addToolTip,
                                              Action removeAction, Icon removeIcon, String removeToolTip, int axis) {

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new BoxLayout(buttonPanel, axis));
        buttonPanel.setOpaque(false);
        JButton addButton = new JButton(addAction);
        if (addIcon != null) {
            addButton.setIcon(addIcon);
            addButton.setText(null);
        }
        addButton.setToolTipText(addToolTip);
        addButton.putClientProperty("JButton.buttonType", "roundRect");
        // addButton.putClientProperty("JButton.buttonType", "toolbar");
        addButton.setOpaque(false);
        addAction.setEnabled(false);

        JButton removeButton = new JButton(removeAction);
        if (removeIcon != null) {
            removeButton.setIcon(removeIcon);
            removeButton.setText(null);
        }
        removeButton.setToolTipText(removeToolTip);
        removeButton.putClientProperty("JButton.buttonType", "roundRect");
//        removeButton.putClientProperty("JButton.buttonType", "toolbar");
        removeButton.setOpaque(false);
        removeAction.setEnabled(false);

        buttonPanel.add(addButton);
        buttonPanel.add(new JToolBar.Separator(new Dimension(6, 6)));
        buttonPanel.add(removeButton);

        return buttonPanel;
    }

    private void excludedSpeciesTableSelectionChanged() {
        if (excludedSpeciesTable.getSelectedRowCount() == 0) {
            includeSpeciesAction.setEnabled(false);
        } else {
            includeSpeciesAction.setEnabled(true);
        }
    }

    private void includedSpeciesTableSelectionChanged() {
        if (includedSpeciesTable.getSelectedRowCount() == 0) {
            excludeSpeciesAction.setEnabled(false);
        } else {
            excludeSpeciesAction.setEnabled(true);
        }
    }

    private void includeSelectedSpecies() {
        int[] rows = excludedSpeciesTable.getSelectedRows();

        List<String> transfer = new ArrayList<String>();

        for (int r : rows) {
            transfer.add(excludedSpecies.get(r));
        }

        includedSpecies.addAll(transfer);
        Collections.sort(includedSpecies);

        excludedSpecies.removeAll(includedSpecies);

        includedSpeciesTableModel.fireTableDataChanged();
        excludedSpeciesTableModel.fireTableDataChanged();

        includedSpeciesTable.getSelectionModel().clearSelection();
        for (String sp : transfer) {
            int row = includedSpecies.indexOf(sp);
            includedSpeciesTable.getSelectionModel().addSelectionInterval(row, row);
        }

        speciesSetChanged();
    }

    private void excludeSelectedSpecies() {
        int[] rows = includedSpeciesTable.getSelectedRows();

        List<String> transfer = new ArrayList<String>();

        for (int r : rows) {
            transfer.add(includedSpecies.get(r));
        }

        excludedSpecies.addAll(transfer);
        Collections.sort(excludedSpecies);

        includedSpecies.removeAll(excludedSpecies);

        includedSpeciesTableModel.fireTableDataChanged();
        excludedSpeciesTableModel.fireTableDataChanged();

        excludedSpeciesTable.getSelectionModel().clearSelection();
        for (String sp : transfer) {
            int row = excludedSpecies.indexOf(sp);
            excludedSpeciesTable.getSelectionModel().addSelectionInterval(row, row);
        }

        speciesSetChanged();
    }

    Action includeSpeciesAction = new AbstractAction("->") {
        public void actionPerformed(ActionEvent ae) {
            includeSelectedSpecies();
        }
    };

    Action excludeSpeciesAction = new AbstractAction("<-") {
        public void actionPerformed(ActionEvent ae) {
            excludeSelectedSpecies();
        }
    };

    class SpeciesTableModel extends AbstractTableModel {
        boolean included;

        public SpeciesTableModel(boolean included) {
            this.included = included;
        }

        public int getColumnCount() {
            return 1;
        }

        public int getRowCount() {
            if (currentSpeciesSet == null) return 0;

            if (included) {
                return includedSpecies.size();
            } else {
                return excludedSpecies.size();
            }
        }

        public Object getValueAt(int row, int col) {

            if (included) {
                return includedSpecies.get(row);
            } else {
                return excludedSpecies.get(row);
            }
        }

        public boolean isCellEditable(int row, int col) {
            return false;
        }

        public String getColumnName(int column) {
            if (included) return "Included Species";
            else return "Excluded Species";
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }
    }

}
