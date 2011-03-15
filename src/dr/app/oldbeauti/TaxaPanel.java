/*
 * TaxonSetPanel.java
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

package dr.app.oldbeauti;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: TaxaPanel.java,v 1.1 2006/09/05 13:29:34 rambaut Exp $
 */
public class TaxaPanel extends JPanel implements Exportable {


    /**
     *
     */
    private static final long serialVersionUID = -3138832889782090814L;

    private final String TAXON_SET_DEFAULT = "taxon set...";

    BeautiFrame frame = null;

    BeautiOptions options = null;

    private TaxonList taxa = null;
    private JTable taxonSetsTable = null;
    private TaxonSetsTableModel taxonSetsTableModel = null;

    private JPanel taxonSetEditingPanel = null;

    private Taxa currentTaxonSet = null;

    private List<Taxon> includedTaxa = new ArrayList<Taxon>();
    private List<Taxon> excludedTaxa = new ArrayList<Taxon>();

    private JTable excludedTaxaTable = null;
    private TaxaTableModel excludedTaxaTableModel = null;
    private JComboBox excludedTaxonSetsComboBox = null;
    private boolean excludedSelectionChanging = false;

    private JTable includedTaxaTable = null;
    private TaxaTableModel includedTaxaTableModel = null;
    private JComboBox includedTaxonSetsComboBox = null;
    private boolean includedSelectionChanging = false;

    private static int taxonSetCount = 0;

    public TaxaPanel(BeautiFrame parent) {

        this.frame = parent;

        Icon includeIcon = null, excludeIcon = null;
        try {
            includeIcon = new ImageIcon(IconUtils.getImage(this.getClass(), "images/include.png"));
            excludeIcon = new ImageIcon(IconUtils.getImage(this.getClass(), "images/exclude.png"));
        } catch (Exception e) {
            // do nothing
        }

        // Taxon Sets
        taxonSetsTableModel = new TaxonSetsTableModel();
        taxonSetsTable = new JTable(taxonSetsTableModel);
        final TableColumnModel model = taxonSetsTable.getColumnModel();
        final TableColumn tableColumn0 = model.getColumn(0);
        tableColumn0.setCellRenderer(new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        tableColumn0.setMinWidth(20);

        //final TableColumn tableColumn1 = model.getColumn(1);

        taxonSetsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) { taxonSetsTableSelectionChanged(); }
        });

        JScrollPane scrollPane1 = new JScrollPane(taxonSetsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addTaxonSetAction);
        actionPanel1.setRemoveAction(removeTaxonSetAction);

        addTaxonSetAction.setEnabled(false);
        removeTaxonSetAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.add(actionPanel1);

        // Excluded Taxon List
        excludedTaxaTableModel = new TaxaTableModel(false);
        excludedTaxaTable = new JTable(excludedTaxaTableModel);

        excludedTaxaTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        excludedTaxaTable.getColumnModel().getColumn(0).setMinWidth(20);

        excludedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) { excludedTaxaTableSelectionChanged(); }
        });

        JScrollPane scrollPane2 = new JScrollPane(excludedTaxaTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        Box panel1 = new Box(BoxLayout.X_AXIS);
        panel1.add(new JLabel("Select: "));
        panel1.setOpaque(false);
        excludedTaxonSetsComboBox = new JComboBox(new String[] { TAXON_SET_DEFAULT });
        excludedTaxonSetsComboBox.setOpaque(false);
        panel1.add(excludedTaxonSetsComboBox);

        JPanel buttonPanel = createAddRemoveButtonPanel(includeTaxonAction, includeIcon, "Include selected taxa in the taxon set",
                excludeTaxonAction, excludeIcon, "Exclude selected taxa from the taxon set",
                javax.swing.BoxLayout.Y_AXIS);

        // Included Taxon List
        includedTaxaTableModel = new TaxaTableModel(true);
        includedTaxaTable = new JTable(includedTaxaTableModel);

        includedTaxaTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        includedTaxaTable.getColumnModel().getColumn(0).setMinWidth(20);

        includedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) { includedTaxaTableSelectionChanged(); }
        });
        includedTaxaTable.doLayout();

        JScrollPane scrollPane3 = new JScrollPane(includedTaxaTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        Box panel2 = new Box(BoxLayout.X_AXIS);
        panel2.add(new JLabel("Select: "));
        panel2.setOpaque(false);
        includedTaxonSetsComboBox = new JComboBox(new String[] { TAXON_SET_DEFAULT });
        includedTaxonSetsComboBox.setOpaque(false);
        panel2.add(includedTaxonSetsComboBox);

        taxonSetEditingPanel = new JPanel();
        taxonSetEditingPanel.setBorder(BorderFactory.createTitledBorder("Taxon Set: none selected"));
        taxonSetEditingPanel.setOpaque(false);
        taxonSetEditingPanel.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0.5;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(12,12,4,0);
        taxonSetEditingPanel.add(scrollPane2, c);

        c.gridx = 0;
        c.gridy = 1;
        c.weightx = 0.5;
        c.weighty = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0,12,12,0);
        taxonSetEditingPanel.add(panel1, c);

        c.gridx = 1;
        c.gridy = 0;
        c.weightx = 0;
        c.weighty = 1;
        c.gridheight = 2;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(12,2,12,4);
        taxonSetEditingPanel.add(buttonPanel, c);

        c.gridx = 2;
        c.gridy = 0;
        c.weightx = 0.5;
        c.weighty = 1;
        c.gridheight = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(12,0,4,12);
        taxonSetEditingPanel.add(scrollPane3, c);

        c.gridx = 2;
        c.gridy = 1;
        c.weightx = 0.5;
        c.weighty = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0,0,12,12);
        taxonSetEditingPanel.add(panel2, c);

        JPanel panel3 = new JPanel();
        panel3.setOpaque(false);
        panel3.setLayout(new GridBagLayout());
        c = new GridBagConstraints();

        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0.4;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0,0,2,12);
        panel3.add(scrollPane1, c);

        c.gridx = 0;
        c.gridy = 1;
        c.weightx = 0;
        c.weighty = 0;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.WEST;
        c.insets = new Insets(2,0,0,12);
        panel3.add(actionPanel1, c);

        c.gridx = 1;
        c.gridy = 0;
        c.weightx = 0.6;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0,0,0,0);
        panel3.add(taxonSetEditingPanel, c);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0,0));
        add(panel3, BorderLayout.CENTER);

//		taxonSetsTable.addMouseListener(new MouseAdapter() {
//			public void mouseClicked(MouseEvent e) {
//				if (e.getClickCount() == 2) {
//					JTable target = (JTable)e.getSource();
//					int row = target.getSelectedRow();
//					taxonSetsTableDoubleClicked(row);
//				}
//			}
//		});

        includedTaxaTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    includeSelectedTaxa();
                }
            }
        });
        excludedTaxaTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    excludeSelectedTaxa();
                }
            }
        });

        includedTaxaTable.addFocusListener(new FocusAdapter() {
            public void focusGained(FocusEvent focusEvent) {
                excludedTaxaTable.clearSelection();
            }
        });
        excludedTaxaTable.addFocusListener(new FocusAdapter() {
            public void focusGained(FocusEvent focusEvent) {
                includedTaxaTable.clearSelection();
            }
        });

        includedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (!includedSelectionChanging && includedTaxonSetsComboBox.getSelectedIndex() != 0) {
                    includedTaxonSetsComboBox.setSelectedIndex(0);
                }
            }
        });
        includedTaxonSetsComboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                includedSelectionChanging = true;
                includedTaxaTable.clearSelection();
                if (includedTaxonSetsComboBox.getSelectedIndex() > 0) {
                    Taxa taxa = (Taxa)includedTaxonSetsComboBox.getSelectedItem();
                    for (int i =0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        int index = includedTaxa.indexOf(taxon);
                        includedTaxaTable.getSelectionModel().addSelectionInterval(index, index);

                    }
                }
                includedSelectionChanging = false;
            }
        });

        excludedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (!excludedSelectionChanging && excludedTaxonSetsComboBox.getSelectedIndex() != 0) {
                    excludedTaxonSetsComboBox.setSelectedIndex(0);
                }
            }
        });
        excludedTaxonSetsComboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                excludedSelectionChanging = true;
                excludedTaxaTable.clearSelection();
                if (excludedTaxonSetsComboBox.getSelectedIndex() > 0) {
                    Taxa taxa = (Taxa)excludedTaxonSetsComboBox.getSelectedItem();
                    for (int i =0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        int index = excludedTaxa.indexOf(taxon);
                        excludedTaxaTable.getSelectionModel().addSelectionInterval(index, index);

                    }
                }
                excludedSelectionChanging = false;
            }
        });

        taxonSetsTable.doLayout();
        includedTaxaTable.doLayout();
        excludedTaxaTable.doLayout();
    }

    private void taxonSetChanged() {
        currentTaxonSet.removeAllTaxa();
        for (Taxon anIncludedTaxa : includedTaxa) {
            currentTaxonSet.addTaxon(anIncludedTaxa);
        }

        setupTaxonSetsComboBoxes();

        if (options.taxonSetsMono.get(currentTaxonSet) != null &&
                options.taxonSetsMono.get(currentTaxonSet) &&
                !checkCompatibility(currentTaxonSet)) {
            options.taxonSetsMono.put(currentTaxonSet, Boolean.FALSE);
        }

        frame.taxonSetsChanged();
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        taxa = options.alignment;
        if (taxa == null) {
            addTaxonSetAction.setEnabled(false);
            removeTaxonSetAction.setEnabled(false);
        } else {
            addTaxonSetAction.setEnabled(true);
        }

        taxonSetsTableSelectionChanged();
        taxonSetsTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
//		options.datesUnits = unitsCombo.getSelectedIndex();
//		options.datesDirection = directionCombo.getSelectedIndex();
//		options.translation = translationCombo.getSelectedIndex();
    }

    public JComponent getExportableComponent() {
        return taxonSetsTable;
    }

    private void taxonSetsTableSelectionChanged() {
        int[] rows = taxonSetsTable.getSelectedRows();
        if (rows.length == 0) {
            removeTaxonSetAction.setEnabled(false);
        } else if (rows.length == 1) {
            currentTaxonSet = options.taxonSets.get(rows[0]);
            setCurrentTaxonSet(currentTaxonSet);
            removeTaxonSetAction.setEnabled(true);
        } else {
            setCurrentTaxonSet(null);
            removeTaxonSetAction.setEnabled(true);
        }
    }

//	private void taxonSetsTableDoubleClicked(int row) {
//		currentTaxonSet = (Taxa)taxonSets.get(row);
//
//		Collections.sort(taxonSets);
//		taxonSetsTableModel.fireTableDataChanged();
//
//		setCurrentTaxonSet(currentTaxonSet);
//
//		int sel = taxonSets.indexOf(currentTaxonSet);
//		taxonSetsTable.setRowSelectionInterval(sel, sel);
//	}

    Action addTaxonSetAction = new AbstractAction("+") {

        /**
         *
         */
        private static final long serialVersionUID = 20273987098143413L;

        public void actionPerformed(ActionEvent ae) {
            taxonSetCount ++;
            currentTaxonSet = new Taxa("untitled" + taxonSetCount);

            options.taxonSets.add(currentTaxonSet);
            Collections.sort(options.taxonSets);

            options.taxonSetsMono.put(currentTaxonSet, Boolean.FALSE);

            taxonSetsTableModel.fireTableDataChanged();

            int sel = options.taxonSets.indexOf(currentTaxonSet);
            taxonSetsTable.setRowSelectionInterval(sel, sel);

            taxonSetChanged();
        }
    };

    Action removeTaxonSetAction = new AbstractAction("-") {

        /**
         *
         */
        private static final long serialVersionUID = 6077578872870122265L;

        public void actionPerformed(ActionEvent ae) {
            int row = taxonSetsTable.getSelectedRow();
            if (row != -1) {
                Taxa taxa = options.taxonSets.remove(row);
                options.taxonSetsMono.remove(taxa);
            }
            taxonSetChanged();

            taxonSetsTableModel.fireTableDataChanged();

            if (row >= options.taxonSets.size()) {
                row = options.taxonSets.size() - 1;
            }
            if (row >= 0) {
                taxonSetsTable.setRowSelectionInterval(row, row);
            } else {
                setCurrentTaxonSet(null);
            }
        }
    };

    private void setCurrentTaxonSet(Taxa taxonSet) {

        this.currentTaxonSet = taxonSet;

        includedTaxa.clear();
        excludedTaxa.clear();

        if (currentTaxonSet != null) {
            for (int i = 0; i < taxonSet.getTaxonCount(); i++) {
                includedTaxa.add(taxonSet.getTaxon(i));
            }
            Collections.sort(includedTaxa);

            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                excludedTaxa.add(taxa.getTaxon(i));
            }
            excludedTaxa.removeAll(includedTaxa);
            Collections.sort(excludedTaxa);
        }

        setTaxonSetTitle();

        setupTaxonSetsComboBoxes();

        includedTaxaTableModel.fireTableDataChanged();
        excludedTaxaTableModel.fireTableDataChanged();
    }

    private void setTaxonSetTitle() {

        if (currentTaxonSet == null) {
            taxonSetEditingPanel.setBorder(BorderFactory.createTitledBorder(""));
            taxonSetEditingPanel.setEnabled(false);
        } else {
            taxonSetEditingPanel.setEnabled(true);
            taxonSetEditingPanel.setBorder(BorderFactory.createTitledBorder("Taxon Set: " + currentTaxonSet.getId()));
        }
    }


    private void setupTaxonSetsComboBoxes() {
        setupTaxonSetsComboBox(excludedTaxonSetsComboBox, excludedTaxa);
        excludedTaxonSetsComboBox.setSelectedIndex(0);
        setupTaxonSetsComboBox(includedTaxonSetsComboBox, includedTaxa);
        includedTaxonSetsComboBox.setSelectedIndex(0);
    }

    private void setupTaxonSetsComboBox(JComboBox comboBox, List availableTaxa) {
        comboBox.removeAllItems();

        comboBox.addItem(TAXON_SET_DEFAULT);
        for (Taxa taxa : options.taxonSets) {
            if (taxa != currentTaxonSet) {
                if (isCompatible(taxa, availableTaxa)) {
                    comboBox.addItem(taxa);
                }
            }
        }
    }

    /**
     * Returns true if taxa are all found in availableTaxa
     * @param taxa
     * @param availableTaxa
     * @return true if the taxa are all found in availableTaxa
     */
    private boolean isCompatible(Taxa taxa, List availableTaxa) {

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);
            if (!availableTaxa.contains(taxon)) {
                return false;
            }
        }
        return true;
    }

    private boolean checkCompatibility(Taxa taxa) {
        for (Taxa taxa2 : options.taxonSets) {
            if (taxa2 != taxa && options.taxonSetsMono.get(taxa2)) {
                if (taxa.containsAny(taxa2) && !taxa.containsAll(taxa2) && !taxa2.containsAll(taxa)) {
                    JOptionPane.showMessageDialog(frame,
                            "You cannot enforce monophyly on this taxon set \n" +
                                    "because it is not compatible with another taxon \n" +
                                    "set, " + taxa2.getId() + ", for which monophyly is\n" +
                                    "enforced.",
                            "Warning",
                            JOptionPane.WARNING_MESSAGE);
                    return false;
                }
            }
        }
        return true;
    }

    class TaxonSetsTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = 3318461381525023153L;

        public TaxonSetsTableModel() {
        }

        public int getColumnCount() {
            return 2;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.taxonSets.size();
        }

        public Object getValueAt(int rowIndex, int columnIndex) {
            Taxa taxonSet = options.taxonSets.get(rowIndex);
            switch(columnIndex) {
                case 0: return taxonSet.getId();
                case 1: return options.taxonSetsMono.get(taxonSet);
            }
            return null;
        }

        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            Taxa taxonSet = options.taxonSets.get(rowIndex);
            switch(columnIndex) {
                case 0: {
                    taxonSet.setId(aValue.toString());
                    setTaxonSetTitle();
                    break;
                }
                case 1: {
                    if ((Boolean)aValue) {
                        Taxa taxa = options.taxonSets.get(rowIndex);
                        if (checkCompatibility(taxa)) {
                            options.taxonSetsMono.put(taxonSet, (Boolean)aValue);
                        }
                    } else {
                        options.taxonSetsMono.put(taxonSet, (Boolean)aValue);
                    }
                    break;
                }
            }
        }

        public boolean isCellEditable(int row, int col) {
            return true;
        }

        public String getColumnName(int column) {
            switch(column) {
                case 0: return "Taxon Sets";
                case 1: return "Monophyletic?";
            }
            return null;
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
        addButton.putClientProperty("JButton.buttonType", "toolbar");
        addButton.setOpaque(false);
        addAction.setEnabled(false);

        JButton removeButton = new JButton(removeAction);
        if (removeIcon != null) {
            removeButton.setIcon(removeIcon);
            removeButton.setText(null);
        }
        removeButton.setToolTipText(removeToolTip);
        removeButton.putClientProperty("JButton.buttonType", "toolbar");
        removeButton.setOpaque(false);
        removeAction.setEnabled(false);

        buttonPanel.add(addButton);
        buttonPanel.add(new JToolBar.Separator(new Dimension(6,6)));
        buttonPanel.add(removeButton);

        return buttonPanel;
    }

    private void excludedTaxaTableSelectionChanged() {
        if (excludedTaxaTable.getSelectedRowCount() == 0) {
            includeTaxonAction.setEnabled(false);
        } else {
            includeTaxonAction.setEnabled(true);
        }
    }

    private void includedTaxaTableSelectionChanged() {
        if (includedTaxaTable.getSelectedRowCount() == 0) {
            excludeTaxonAction.setEnabled(false);
        } else {
            excludeTaxonAction.setEnabled(true);
        }
    }

    private void includeSelectedTaxa() {
        int[] rows = excludedTaxaTable.getSelectedRows();

        List<Taxon> transfer = new ArrayList<Taxon>();

        for (int r : rows) {
            transfer.add(excludedTaxa.get(r));
        }

        includedTaxa.addAll(transfer);
        Collections.sort(includedTaxa);

        excludedTaxa.removeAll(includedTaxa);

        includedTaxaTableModel.fireTableDataChanged();
        excludedTaxaTableModel.fireTableDataChanged();

        includedTaxaTable.getSelectionModel().clearSelection();
        for (Taxon taxon : transfer) {
            int row = includedTaxa.indexOf(taxon);
            includedTaxaTable.getSelectionModel().addSelectionInterval(row, row);
        }

        taxonSetChanged();
    }

    private void excludeSelectedTaxa() {
        int[] rows = includedTaxaTable.getSelectedRows();

        List<Taxon> transfer = new ArrayList<Taxon>();

        for (int r : rows) {
            transfer.add(includedTaxa.get(r));
        }

        excludedTaxa.addAll(transfer);
        Collections.sort(excludedTaxa);

        includedTaxa.removeAll(excludedTaxa);

        includedTaxaTableModel.fireTableDataChanged();
        excludedTaxaTableModel.fireTableDataChanged();

        excludedTaxaTable.getSelectionModel().clearSelection();
        for (Taxon taxon : transfer) {
            int row = excludedTaxa.indexOf(taxon);
            excludedTaxaTable.getSelectionModel().addSelectionInterval(row, row);
        }

        taxonSetChanged();
    }

    Action includeTaxonAction = new AbstractAction("->") {
        /**
         *
         */
        private static final long serialVersionUID = 7510299673661594128L;

        public void actionPerformed(ActionEvent ae) {
            includeSelectedTaxa();
        }
    };

    Action excludeTaxonAction = new AbstractAction("<-") {

        /**
         *
         */
        private static final long serialVersionUID = 449692708602410206L;

        public void actionPerformed(ActionEvent ae) {
            excludeSelectedTaxa();
        }
    };

    class TaxaTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -8027482229525938010L;
        boolean included;

        public TaxaTableModel(boolean included) {
            this.included = included;
        }

        public int getColumnCount() {
            return 1;
        }

        public int getRowCount() {
            if (currentTaxonSet == null) return 0;

            if (included) {
                return includedTaxa.size();
            } else {
                return excludedTaxa.size();
            }
        }

        public Object getValueAt(int row, int col) {

            if (included) {
                return includedTaxa.get(row).getId();
            } else {
                return excludedTaxa.get(row).getId();
            }
        }

        public boolean isCellEditable(int row, int col) {
            return false;
        }

        public String getColumnName(int column) {
            if (included) return "Included Taxa";
            else return "Excluded Taxa";
        }

        public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
    }

}
