/*
 * TaxonSetPanel.java
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

import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.DateCellEditor;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;
import jam.util.IconUtils;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TaxaPanel.java,v 1.1 2006/09/05 13:29:34 rambaut Exp $
 */
public class TaxonSetsPanel extends JPanel implements Exportable {

	private static final long serialVersionUID = -3138832889782090814L;
	
	protected String TAXA;
	protected String TAXON;

	TreeStatFrame frame = null;
	TreeStatData treeStatData = null;
	TreeStatData.TaxonSet selectedTaxonSet = null;

	//    private TaxonList taxa = null;
	protected JTable taxonSetsTable = null;
	private TableColumnModel tableColumnModel;
	protected TaxonSetsTableModel taxonSetsTableModel = new TaxonSetsTableModel();
	ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();

	protected JPanel taxonSetEditingPanel = null;

	protected TreeStatData.TaxonSet currentTaxonSet = null;

	protected final List<Taxon> includedTaxa = new ArrayList<>();
	protected final List<Taxon> excludedTaxa = new ArrayList<>();

	private JTextField excludedTaxaSearchField = new JTextField();

	protected JTable excludedTaxaTable = null;
	protected TaxaTableModel excludedTaxaTableModel = null;
	private JLabel excludedTaxaLabel = new JLabel();
	protected JComboBox excludedTaxonSetsComboBox = null;
	protected boolean excludedSelectionChanging = false;

	private JTextField includedTaxaSearchField = new JTextField();

	protected JTable includedTaxaTable = null;
	protected TaxaTableModel includedTaxaTableModel = null;
	private JLabel includedTaxaLabel = new JLabel();
	protected JComboBox includedTaxonSetsComboBox = null;
	protected boolean includedSelectionChanging = false;

	public TaxonSetsPanel(TreeStatFrame parent, TreeStatData treeStatData) {

		this.frame = parent;

		this.treeStatData = treeStatData;
		setText(false);

		// Taxon Sets
		initTaxonSetsTable(taxonSetsTableModel);

		initTableColumn();

		initPanel(addTaxonSetAction, removeTaxonSetAction);
	}

	protected void setText(boolean useStarBEAST) {
		if (useStarBEAST) {
			TAXA = "Species";
			TAXON = "Species set";
		} else {
			TAXA = "Taxa";
			TAXON = "Taxon set";
		}
	}

	protected void initPanel(Action addTaxonSetAction, Action removeTaxonSetAction) {
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
			public void valueChanged(ListSelectionEvent evt) {
				excludedTaxaTableSelectionChanged();
			}
		});

		JScrollPane scrollPane2 = new JScrollPane(excludedTaxaTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		includedTaxonSetsComboBox = new JComboBox(new String[]{TAXON.toLowerCase() + "..."});
		excludedTaxonSetsComboBox = new JComboBox(new String[]{TAXON.toLowerCase() + "..."});

		includedTaxaLabel.setText("");
		excludedTaxaLabel.setText("");

		Box panel1 = new Box(BoxLayout.X_AXIS);
		panel1.add(new JLabel("Select: "));
		panel1.setOpaque(false);
		excludedTaxonSetsComboBox.setOpaque(false);
		panel1.add(excludedTaxonSetsComboBox);

		// Included Taxon List
		includedTaxaTableModel = new TaxaTableModel(true);
		includedTaxaTable = new JTable(includedTaxaTableModel);

		includedTaxaTable.getColumnModel().getColumn(0).setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
		includedTaxaTable.getColumnModel().getColumn(0).setMinWidth(20);

		includedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) {
				includedTaxaTableSelectionChanged();
			}
		});
		includedTaxaTable.doLayout();

		JScrollPane scrollPane3 = new JScrollPane(includedTaxaTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		Box panel2 = new Box(BoxLayout.X_AXIS);
		panel2.add(new JLabel("Select: "));
		panel2.setOpaque(false);
		includedTaxonSetsComboBox.setOpaque(false);
		panel2.add(includedTaxonSetsComboBox);

		Icon includeIcon = null, excludeIcon = null;
		try {
			includeIcon = new ImageIcon(IconUtils.getImage(TreeStatApp.class, "images/include.png"));
			excludeIcon = new ImageIcon(IconUtils.getImage(TreeStatApp.class, "images/exclude.png"));
		} catch (Exception e) {
			// do nothing
		}

		JPanel buttonPanel = createAddRemoveButtonPanel(includeTaxonAction, includeIcon, "Include selected "
						+ TAXA.toLowerCase() + " in the " + TAXON.toLowerCase(),
				excludeTaxonAction, excludeIcon, "Exclude selected " + TAXA.toLowerCase()
						+ " from the " + TAXON.toLowerCase(), BoxLayout.Y_AXIS);

		taxonSetEditingPanel = new JPanel();
		taxonSetEditingPanel.setBorder(BorderFactory.createTitledBorder(""));
		taxonSetEditingPanel.setOpaque(false);
		taxonSetEditingPanel.setLayout(new GridBagLayout());

		excludedTaxaSearchField.setColumns(12);
//        excludedTaxaSearchField.putClientProperty("JTextField.variant", "search");
		excludedTaxaSearchField.putClientProperty("Quaqua.TextField.style","search");
		excludedTaxaSearchField.putClientProperty("Quaqua.TextField.sizeVariant","small");
		includedTaxaSearchField.setColumns(12);
//        includedTaxaSearchField.putClientProperty("JTextField.variant", "search");
		includedTaxaSearchField.putClientProperty("Quaqua.TextField.style","search");
		includedTaxaSearchField.putClientProperty("Quaqua.TextField.sizeVariant","small");

		GridBagConstraints c = new GridBagConstraints();

		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0.5;
		c.weighty = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(3, 6, 3, 0);
		taxonSetEditingPanel.add(excludedTaxaSearchField, c);

		c.gridx = 0;
		c.gridy = 1;
		c.weightx = 0.5;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 6, 0, 0);
		taxonSetEditingPanel.add(scrollPane2, c);

		c.gridx = 0;
		c.gridy = 2;
		c.weightx = 0.5;
		c.weighty = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 6, 3, 0);
		taxonSetEditingPanel.add(excludedTaxaLabel, c);

		c.gridx = 0;
		c.gridy = 3;
		c.weightx = 0.5;
		c.weighty = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 6, 3, 0);
		taxonSetEditingPanel.add(panel1, c);

		c.gridx = 1;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 1;
		c.gridheight = 4;
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(12, 2, 12, 4);
		taxonSetEditingPanel.add(buttonPanel, c);

		c.gridx = 2;
		c.gridy = 0;
		c.weightx = 0.5;
		c.weighty = 0;
		c.gridheight = 1;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(3, 0, 3, 6);
		taxonSetEditingPanel.add(includedTaxaSearchField, c);

		c.gridx = 2;
		c.gridy = 1;
		c.weightx = 0.5;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 0, 0, 6);
		taxonSetEditingPanel.add(scrollPane3, c);

		c.gridx = 2;
		c.gridy = 2;
		c.weightx = 0.5;
		c.weighty = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 0, 3, 6);
		taxonSetEditingPanel.add(includedTaxaLabel, c);

		c.gridx = 2;
		c.gridy = 3;
		c.weightx = 0.5;
		c.weighty = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 0, 3, 6);
		taxonSetEditingPanel.add(panel2, c);

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
		panel3.add(taxonSetEditingPanel, c);

		setOpaque(false);
		setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
		setLayout(new BorderLayout(0, 0));
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

		includedTaxaSearchField.getDocument().addDocumentListener(new DocumentListener() {
																	  public void changedUpdate(DocumentEvent e) {
																		  selectIncludedTaxa(includedTaxaSearchField.getText());
																	  }

																	  public void removeUpdate(DocumentEvent e) {
																		  selectIncludedTaxa(includedTaxaSearchField.getText());
																	  }

																	  public void insertUpdate(DocumentEvent e) {
																		  selectIncludedTaxa(includedTaxaSearchField.getText());
																	  }
																  }
		);
		excludedTaxaSearchField.getDocument().addDocumentListener(new DocumentListener() {
																	  public void changedUpdate(DocumentEvent e) {
																		  selectExcludedTaxa(excludedTaxaSearchField.getText());
																	  }

																	  public void removeUpdate(DocumentEvent e) {
																		  selectExcludedTaxa(excludedTaxaSearchField.getText());
																	  }

																	  public void insertUpdate(DocumentEvent e) {
																		  selectExcludedTaxa(excludedTaxaSearchField.getText());
																	  }
																  }
		);


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
				if (!includedSelectionChanging) {
					if (includedTaxonSetsComboBox.getSelectedIndex() != 0) {
						includedTaxonSetsComboBox.setSelectedIndex(0);
					}
					includedTaxaSearchField.setText("");
				}
			}
		});
		includedTaxonSetsComboBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				includedSelectionChanging = true;
				includedTaxaTable.clearSelection();
				if (includedTaxonSetsComboBox.getSelectedIndex() > 0) {
					String taxaName = includedTaxonSetsComboBox.getSelectedItem().toString();
					if (!taxaName.endsWith("...")) {
						TreeStatData.TaxonSet taxonSet = treeStatData.taxonSets.get(taxaName);
						if (taxonSet != null) {
							for (int i = 0; i < taxonSet.taxa.getTaxonCount(); i++) {
								Taxon taxon = taxonSet.taxa.getTaxon(i);
								int index = includedTaxa.indexOf(taxon);
								includedTaxaTable.getSelectionModel().addSelectionInterval(index, index);

							}
						}
					}
				}
				includedSelectionChanging = false;
			}
		});

		excludedTaxaTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				if (!excludedSelectionChanging) {
					if (excludedTaxonSetsComboBox.getSelectedIndex() != 0) {
						excludedTaxonSetsComboBox.setSelectedIndex(0);
					}
					excludedTaxaSearchField.setText("");
				}

			}
		});
		excludedTaxonSetsComboBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				excludedSelectionChanging = true;
				excludedTaxaTable.clearSelection();
				if (excludedTaxonSetsComboBox.getSelectedIndex() > 0) {
					String taxaName = excludedTaxonSetsComboBox.getSelectedItem().toString();
					if (!taxaName.endsWith("...")) {
						TreeStatData.TaxonSet taxonSet = treeStatData.taxonSets.get(taxaName);
						if (taxonSet != null) {
							for (int i = 0; i < taxonSet.taxa.getTaxonCount(); i++) {
								Taxon taxon = taxonSet.taxa.getTaxon(i);
								int index = excludedTaxa.indexOf(taxon);
								excludedTaxaTable.getSelectionModel().addSelectionInterval(index, index);

							}
						}
					}
				}
				excludedSelectionChanging = false;
			}
		});

		includedTaxaTable.doLayout();
		excludedTaxaTable.doLayout();
	}

	private void selectIncludedTaxa(String text) {
		includedSelectionChanging = true;
		includedTaxaTable.clearSelection();
		int index = 0;
		for (Taxon taxon : includedTaxa) {
			if (taxon.getId().contains(text)) {
				includedTaxaTable.getSelectionModel().addSelectionInterval(index, index);
			}
			index ++;

		}
		includedSelectionChanging = false;
	}

	private void selectExcludedTaxa(String text) {
		excludedSelectionChanging = true;
		excludedTaxaTable.clearSelection();
		int index = 0;
		for (Taxon taxon : excludedTaxa) {
			if (taxon.getId().contains(text)) {
				excludedTaxaTable.getSelectionModel().addSelectionInterval(index, index);
			}
			index ++;

		}
		excludedSelectionChanging = false;

	}

	protected void initTableColumn() {
		tableColumnModel = taxonSetsTable.getColumnModel();
		TableColumn tableColumn = tableColumnModel.getColumn(0);
		tableColumn.setCellRenderer(new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
		tableColumn.setMinWidth(20);
	}

	protected void initTaxonSetsTable(AbstractTableModel tableModel) {
		taxonSetsTable = new JTable(tableModel);
		taxonSetsTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		taxonSetsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) {
				taxonSetsTableSelectionChanged();
			}
		});
		taxonSetsTable.doLayout();
	}

	protected void taxonSetChanged() {
		currentTaxonSet.taxa.removeAllTaxa();
		for (Taxon anIncludedTaxa : includedTaxa) {
			currentTaxonSet.taxa.addTaxon(anIncludedTaxa);
		}

		setupTaxonSetsComboBoxes();

		includedTaxaLabel.setText("" + includedTaxa.size() + " taxa included");
		excludedTaxaLabel.setText("" + excludedTaxa.size() + " taxa excluded");

		frame.setDirty();
	}

	protected void resetPanel() {
//		if (!treeStatData.hasData() || treeStatData.taxonSets == null || treeStatData.taxonSets.size() < 1) {
//			setCurrentTaxonSet(null);
//		}
	}

	public JComponent getExportableComponent() {
		return taxonSetsTable;
	}

	private void taxonSetsTableSelectionChanged() {
		if (taxonSetsTable.getSelectedRowCount() == 0) {
			selectedTaxonSet = null;
			removeTaxonSetAction.setEnabled(false);
		} else {
			String name = treeStatData.taxonSetNames.get(taxonSetsTable.getSelectedRow());
			selectedTaxonSet = treeStatData.taxonSets.get(name);
			removeTaxonSetAction.setEnabled(true);
		}
		setCurrentTaxonSet(selectedTaxonSet);
		includedTaxaTableModel.fireTableDataChanged();
		excludedTaxaTableModel.fireTableDataChanged();
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

		private static final long serialVersionUID = 20273987098143413L;

		public void actionPerformed(ActionEvent ae) {
			int taxonSetCount = treeStatData.taxonSets.size();
			TreeStatData.TaxonSet taxonSet = new TreeStatData.TaxonSet();
			taxonSet.name = "untitled" + taxonSetCount;
			taxonSet.taxa = new Taxa();
			treeStatData.taxonSetNames.add(taxonSet.name);
			treeStatData.taxonSets.put(taxonSet.name, taxonSet);
			dataChanged();

			int sel = treeStatData.taxonSets.size() - 1;
			taxonSetsTable.setRowSelectionInterval(sel, sel);

			currentTaxonSet = taxonSet;

			taxonSetChanged();

			taxonSetsTableModel.fireTableDataChanged();
		}
	};

	Action removeTaxonSetAction = new AbstractAction("-") {

		private static final long serialVersionUID = 6077578872870122265L;

		public void actionPerformed(ActionEvent ae) {
			int row = taxonSetsTable.getSelectedRow();
			if (row != -1) {
				String name = treeStatData.taxonSetNames.get(row);
				TreeStatData.TaxonSet taxonSet = treeStatData.taxonSets.remove(name);
			}
			taxonSetChanged();

			taxonSetsTableModel.fireTableDataChanged();

			if (row >= treeStatData.taxonSets.size()) {
				row = treeStatData.taxonSets.size() - 1;
			}
			if (row >= 0) {
				taxonSetsTable.setRowSelectionInterval(row, row);
			} else {
				setCurrentTaxonSet(null);
			}
		}
	};

	protected void setCurrentTaxonSet(TreeStatData.TaxonSet taxonSet) {

		this.currentTaxonSet = taxonSet;

		includedTaxa.clear();
		excludedTaxa.clear();

		if (currentTaxonSet != null) {
			for (Taxon taxon : taxonSet.taxa) {
				includedTaxa.add(taxon);
			}
			Collections.sort(includedTaxa);

			for (Taxon taxon : treeStatData.allTaxa) {
				excludedTaxa.add(taxon);
			}
			excludedTaxa.removeAll(includedTaxa);
			Collections.sort(excludedTaxa);
		}

		setTaxonSetTitle();

		setupTaxonSetsComboBoxes();

		includedTaxaTableModel.fireTableDataChanged();
		excludedTaxaTableModel.fireTableDataChanged();
	}

	protected void setTaxonSetTitle() {

		if (currentTaxonSet == null) {
			taxonSetEditingPanel.setBorder(BorderFactory.createTitledBorder(""));
			taxonSetEditingPanel.setEnabled(false);
		} else {
			taxonSetEditingPanel.setEnabled(true);
			taxonSetEditingPanel.setBorder(new TitledBorder(null, TAXON + ": " + currentTaxonSet.name, TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.ABOVE_TOP));
		}
	}


	protected void setupTaxonSetsComboBoxes() {
		setupTaxonSetsComboBox(excludedTaxonSetsComboBox, excludedTaxa);
		excludedTaxonSetsComboBox.setSelectedIndex(0);
		setupTaxonSetsComboBox(includedTaxonSetsComboBox, includedTaxa);
		includedTaxonSetsComboBox.setSelectedIndex(0);
	}

	protected void setupTaxonSetsComboBox(JComboBox comboBox, List<Taxon> availableTaxa) {
		comboBox.removeAllItems();

		comboBox.addItem(TAXON.toLowerCase() + "...");
		for (TreeStatData.TaxonSet taxonSet : treeStatData.taxonSets.values()) {
			// AR - as these comboboxes are just intended to be handy ways of selecting taxa, I have removed
			// these requirements (it was just confusing why they weren't in the lists.
//       if (taxa != currentTaxonSet) {
//                if (isCompatible(taxa, availableTaxa)) {
			comboBox.addItem(taxonSet.name); // have to add String, otherwise it will throw Exception to cast "taxa..." into Taxa
//                }
//            }
		}
	}

	/**
	 * The table on the left side of panel
	 */
	protected class TaxonSetsTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 3318461381525023153L;

		String[] columnNames = {"Taxon Set"};

		public TaxonSetsTableModel() {
		}

		public int getColumnCount() {
			return columnNames.length;
		}

		public String getColumnName(int column) {
			return columnNames[column];
		}

		public int getRowCount() {
			if (treeStatData == null) return 0;
			return treeStatData.taxonSets.size();
		}

		public Object getValueAt(int rowIndex, int columnIndex) {
			String name = treeStatData.taxonSetNames.get(rowIndex);
			TreeStatData.TaxonSet taxonSet = treeStatData.taxonSets.get(name);
			switch (columnIndex) {
				case 0:
					return taxonSet.name;
				default:
					throw new IllegalArgumentException("unknown column, " + columnIndex);
			}
		}

		public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
//			Taxa taxonSet = treeStatData.taxonSets.get(rowIndex);
			TreeStatData.TaxonSet taxonSet = treeStatData.taxonSets.get(rowIndex);
			switch (columnIndex) {
				case 0:
					taxonSet.name=aValue.toString();
					setTaxonSetTitle();
					break;

				default:
					throw new IllegalArgumentException("unknown column, " + columnIndex);
			}
		}

		public boolean isCellEditable(int row, int col) {
			return true;
		}

		public Class getColumnClass(int columnIndex) {
			switch (columnIndex) {
				case 0:
					return String.class;
				default:
					throw new IllegalArgumentException("unknown column, " + columnIndex);
			}
		}
	}

	protected JPanel createAddRemoveButtonPanel(Action addAction, Icon addIcon, String addToolTip,
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

	public void dataChanged() {

		addTaxonSetAction.setEnabled(treeStatData.allTaxa.size() > 0);

		taxonSetsTableModel.fireTableDataChanged();
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
			if (included) return "Included " + TAXA;
			else return "Excluded " + TAXA;
		}

		public Class getColumnClass(int c) {
			return getValueAt(0, c).getClass();
		}
	}

}
