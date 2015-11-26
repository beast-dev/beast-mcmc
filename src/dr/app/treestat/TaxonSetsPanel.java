/*
 * TaxonSetsPanel.java
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

package dr.app.treestat;

import jam.framework.Exportable;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;


public class TaxonSetsPanel extends JPanel implements Exportable {


	/**
	 *
	 */
	private static final long serialVersionUID = -9013475414423166476L;
	TreeStatFrame frame = null;
	TreeStatData treeStatData = null;
	TreeStatData.TaxonSet selectedTaxonSet = null;

	JScrollPane scrollPane1 = null;
	JTable taxonSetsTable = null;
	TaxonSetsTableModel taxonSetsTableModel = null;

	JScrollPane scrollPane2 = null;
	JTable excludedTaxaTable = null;
	TaxaTableModel excludedTaxaTableModel = null;

	JScrollPane scrollPane3 = null;
	JTable includedTaxaTable = null;
	TaxaTableModel includedTaxaTableModel = null;

	public TaxonSetsPanel(TreeStatFrame frame, TreeStatData treeStatData) {

		this.frame = frame;
		setOpaque(false);

		this.treeStatData = treeStatData;

		Icon addIcon = null, removeIcon = null, includeIcon = null, excludeIcon = null;
 		try {
			addIcon = new ImageIcon(dr.app.util.Utils.getImage(this, "images/add.png"));
			removeIcon = new ImageIcon(dr.app.util.Utils.getImage(this, "images/minus.png"));
			includeIcon = new ImageIcon(dr.app.util.Utils.getImage(this, "images/include.png"));
			excludeIcon = new ImageIcon(dr.app.util.Utils.getImage(this, "images/exclude.png"));
		} catch (Exception e) {
			// do nothing
		}

        JPanel buttonPanel = new JPanel(new BorderLayout());
        buttonPanel.setOpaque(false);
        JButton importButton = new JButton(frame.getImportAction());
        importButton.setFocusable(false);
        importButton.putClientProperty("JButton.buttonType", "textured");
        importButton.setMargin(new Insets(4,4,4,4));
        buttonPanel.add(importButton, BorderLayout.WEST);
        buttonPanel.add(new JLabel("    To define taxon sets, first import a list of taxa (i.e., from the trees to be analysed)"), BorderLayout.SOUTH);

        // Taxon Sets
		taxonSetsTableModel = new TaxonSetsTableModel();
        TableSorter sorter = new TableSorter(taxonSetsTableModel);
		taxonSetsTable = new JTable(sorter);
        sorter.addTableModelListener(taxonSetsTable);

		taxonSetsTable.getColumnModel().getColumn(0).setCellRenderer(
			new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

		taxonSetsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) { taxonSetsTableSelectionChanged(); }
		});

   		scrollPane1 = new JScrollPane(taxonSetsTable,
										JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
										JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

 		JPanel buttonPanel1 = createAddRemoveButtonPanel(addTaxonSetAction, addIcon, "Create a new taxon set",
 															removeTaxonSetAction, removeIcon, "Remove a taxon set",
 															javax.swing.BoxLayout.X_AXIS);

		// Excluded Taxon List
		excludedTaxaTableModel = new TaxaTableModel(false);
        sorter = new TableSorter(excludedTaxaTableModel);
		excludedTaxaTable = new JTable(sorter);
        sorter.addTableModelListener(excludedTaxaTable);

		excludedTaxaTable.getColumnModel().getColumn(0).setCellRenderer(
			new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

		excludedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) { excludedTaxaTableSelectionChanged(); }
		});

   		scrollPane2 = new JScrollPane(excludedTaxaTable,
										JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
										JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

 		JPanel buttonPanel2 = createAddRemoveButtonPanel(includeTaxonAction, includeIcon, "Include selected taxa in the taxon set",
 															excludeTaxonAction, excludeIcon, "Exclude selected taxa from the taxon set",
 															javax.swing.BoxLayout.Y_AXIS);

 		// Included Taxon List
		includedTaxaTableModel = new TaxaTableModel(true);
        sorter = new TableSorter(includedTaxaTableModel);
		includedTaxaTable = new JTable(sorter);
        sorter.addTableModelListener(includedTaxaTable);

		includedTaxaTable.getColumnModel().getColumn(0).setCellRenderer(
			new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

		includedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) { includedTaxaTableSelectionChanged(); }
		});

   		scrollPane3 = new JScrollPane(includedTaxaTable,
										JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
										JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        JPanel panel = new JPanel();
        panel.setOpaque(false);
        panel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

        c.weightx = 0.333333;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(6,6,6,6);
		c.gridx = 0;
		c.gridy = 0;
		panel.add(scrollPane1, c);

		c.weightx = 0.333333;
		c.weighty = 0;
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.LINE_START;
		c.insets = new Insets(0,6,6,6);
		c.gridx = 0;
		c.gridy = 1;
		panel.add(buttonPanel1, c);

		c.weightx = 0.333333;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(6,6,6,6);
		c.gridx = 1;
		c.gridy = 0;
		panel.add(scrollPane2, c);

		c.weightx = 0.0;
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0,0,0,0);
		c.gridx = 2;
		c.gridy = 0;
		panel.add(buttonPanel2, c);

		c.weightx = 0.333333;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(6,6,6,6);
		c.gridx = 3;
		c.gridy = 0;
		panel.add(scrollPane3, c);

        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(6,6,6,6));
        add(buttonPanel, BorderLayout.NORTH);
        add(panel, BorderLayout.CENTER);

	}

	JPanel createAddRemoveButtonPanel(Action addAction, Icon addIcon, String addToolTip,
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
		addButton.putClientProperty("JButton.buttonType", "textured");
		addButton.setOpaque(false);
		addAction.setEnabled(false);

		JButton removeButton = new JButton(removeAction);
		if (removeIcon != null) {
			removeButton.setIcon(removeIcon);
			removeButton.setText(null);
		}
   	  	removeButton.setToolTipText(removeToolTip);
		removeButton.putClientProperty("JButton.buttonType", "textured");
		removeButton.setOpaque(false);
		removeAction.setEnabled(false);

		buttonPanel.add(addButton);
		buttonPanel.add(new JToolBar.Separator(new Dimension(6,6)));
		buttonPanel.add(removeButton);

		return buttonPanel;
	}

	public void dataChanged() {

		addTaxonSetAction.setEnabled(treeStatData.allTaxa.size() > 0);

		taxonSetsTableModel.fireTableDataChanged();
		includedTaxaTableModel.fireTableDataChanged();
		excludedTaxaTableModel.fireTableDataChanged();
	}

	private void taxonSetsTableSelectionChanged() {
		if (taxonSetsTable.getSelectedRowCount() == 0) {
			selectedTaxonSet = null;
			removeTaxonSetAction.setEnabled(false);
		} else {
			selectedTaxonSet = treeStatData.taxonSets.get(taxonSetsTable.getSelectedRow());
			removeTaxonSetAction.setEnabled(true);
		}
		includedTaxaTableModel.fireTableDataChanged();
		excludedTaxaTableModel.fireTableDataChanged();
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

    public JComponent getExportableComponent() {
		return this;
	}

  	Action addTaxonSetAction = new AbstractAction("+") {

  		/**
		 *
		 */
		private static final long serialVersionUID = 1831933175582860833L;

		public void actionPerformed(ActionEvent ae) {
			TreeStatData.TaxonSet taxonSet = new TreeStatData.TaxonSet();
			taxonSet.name = "untitled";
			taxonSet.taxa = new ArrayList();
			treeStatData.taxonSets.add(taxonSet);
			dataChanged();

			int sel = treeStatData.taxonSets.size() - 1;
			taxonSetsTable.setRowSelectionInterval(sel, sel);
  		}
  	};

  	Action removeTaxonSetAction = new AbstractAction("-") {

  		/**
		 *
		 */
		private static final long serialVersionUID = -8662527333546044639L;

		public void actionPerformed(ActionEvent ae) {
			int saved = taxonSetsTable.getSelectedRow();
			int row = taxonSetsTable.getSelectedRow();
			if (row != -1) {
				treeStatData.taxonSets.remove(row);
			}
			dataChanged();
			if (saved >= treeStatData.taxonSets.size()) saved = treeStatData.taxonSets.size() - 1;
			taxonSetsTable.setRowSelectionInterval(saved, saved);
  		}
  	};

  	Action includeTaxonAction = new AbstractAction("->") {
  		/**
		 *
		 */
		private static final long serialVersionUID = -1875904513948242608L;

		public void actionPerformed(ActionEvent ae) {
			int saved = taxonSetsTable.getSelectedRow();
			int[] rows = excludedTaxaTable.getSelectedRows();
			ArrayList<String> exclList = new ArrayList<String>(treeStatData.allTaxa);
			exclList.removeAll(selectedTaxonSet.taxa);
            for (int row : rows) {
                selectedTaxonSet.taxa.add(exclList.get(row));
            }
			dataChanged();
			taxonSetsTable.setRowSelectionInterval(saved, saved);
  		}
  	};

  	Action excludeTaxonAction = new AbstractAction("<-") {

  		/**
		 *
		 */
		private static final long serialVersionUID = 4523480086490780822L;

		public void actionPerformed(ActionEvent ae) {
			int saved = taxonSetsTable.getSelectedRow();
			int[] rows = includedTaxaTable.getSelectedRows();
			for (int i = rows.length - 1; i >= 0 ; i--) {
				selectedTaxonSet.taxa.remove(rows[i]);
			}
			dataChanged();
			taxonSetsTable.setRowSelectionInterval(saved, saved);
  		}
  	};

	class TaxonSetsTableModel extends AbstractTableModel {

		/**
		 *
		 */
		private static final long serialVersionUID = 219223813257870207L;

		public TaxonSetsTableModel() {
		}

		public int getColumnCount() {
			return 1;
		}

		public int getRowCount() {
			return treeStatData.taxonSets.size();
		}

		public Object getValueAt(int row, int col) {
			return (treeStatData.taxonSets.get(row)).name;
		}

		public void setValueAt(Object value, int row, int col) {
			(treeStatData.taxonSets.get(row)).name = (String)value;
		}

        public boolean isCellEditable(int row, int col) {
 			return true;
        }

		public String getColumnName(int column) {
			return "Taxon Sets";
		}

		public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
	}

    class TaxaTableModel extends AbstractTableModel {

		/**
		 *
		 */
		private static final long serialVersionUID = 1559408662356843275L;
		boolean included;

		public TaxaTableModel(boolean included) {
			this.included = included;
		}

		public int getColumnCount() {
			return 1;
		}

		public int getRowCount() {
			if (selectedTaxonSet == null) return 0;

			if (included) {
				return selectedTaxonSet.taxa.size();
			} else {
				return treeStatData.allTaxa.size() - selectedTaxonSet.taxa.size();
			}
		}

		public Object getValueAt(int row, int col) {

			if (included) {
				return selectedTaxonSet.taxa.get(row);
			} else {
				ArrayList<String> exclList = new ArrayList<String>(treeStatData.allTaxa);
				exclList.removeAll(selectedTaxonSet.taxa);
				return exclList.get(row);
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
