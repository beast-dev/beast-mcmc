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

package dr.app.beauti.taxonsetspanel;

import dr.app.beauti.BeautiApp;
import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.DateCellEditor;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.DataTable;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TaxaPanel.java,v 1.1 2006/09/05 13:29:34 rambaut Exp $
 */
public class TaxonSetPanel extends BeautiPanel implements Exportable {

    private static final long serialVersionUID = -3138832889782090814L;

    private final String[] columnToolTips = {null,
            "Enforce the selected taxon set to be monophyletic on the specified tree",
            "The tmrcaStatistic will represent that age of the parent node of the MRCA, rather than the MRCA itself",
            "Select the tree from which to report the MRCA of the taxa"};

    protected String TAXA;
    protected String TAXON;

    protected BeautiFrame frame = null;
    protected BeautiOptions options = null;

    ImportTaxonSetsAction importTaxonSetsAction = new ImportTaxonSetsAction();

    //    private TaxonList taxa = null;
    protected JTable taxonSetsTable = null;
    private TableColumnModel tableColumnModel;
    protected TaxonSetsTableModel taxonSetsTableModel = new TaxonSetsTableModel();
    ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();

    protected JPanel taxonSetEditingPanel = null;

    protected Taxa currentTaxonSet = null;

    protected final List<Taxon> includedTaxa = new ArrayList<Taxon>();
    protected final List<Taxon> excludedTaxa = new ArrayList<Taxon>();

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

    protected static int taxonSetCount = 0;

    public TaxonSetPanel() {
    }

    public TaxonSetPanel(BeautiFrame parent) {

        this.frame = parent;

        setText(false);

        // Taxon Sets
        initTaxonSetsTable(taxonSetsTableModel, columnToolTips);

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

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setBorder(BorderFactory.createEmptyBorder());

        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        JButton button = new JButton(importTaxonSetsAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

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
            includeIcon = new ImageIcon(IconUtils.getImage(BeautiApp.class, "images/include.png"));
            excludeIcon = new ImageIcon(IconUtils.getImage(BeautiApp.class, "images/exclude.png"));
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
        excludedTaxaSearchField.putClientProperty("JTextField.variant", "search");
        excludedTaxaSearchField.putClientProperty("Quaqua.TextField.style","search");
        excludedTaxaSearchField.putClientProperty("Quaqua.TextField.sizeVariant","small");
        includedTaxaSearchField.setColumns(12);
        includedTaxaSearchField.putClientProperty("JTextField.variant", "search");
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
        add(toolBar1, BorderLayout.NORTH);
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
                        Taxa taxa = options.getTaxa(taxaName);
                        if (taxa != null) {
                            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                                Taxon taxon = taxa.getTaxon(i);
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
                        Taxa taxa = options.getTaxa(taxaName);
                        if (taxa != null) {
                            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                                Taxon taxon = taxa.getTaxon(i);
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

        tableColumn = tableColumnModel.getColumn(1);
        tableColumn.setPreferredWidth(20);
        tableColumn = tableColumnModel.getColumn(2);
        tableColumn.setPreferredWidth(20);

        tableColumn = tableColumnModel.getColumn(3);
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        tableColumn.setCellRenderer(comboBoxRenderer);
        tableColumn.setPreferredWidth(30);

        tableColumn = tableColumnModel.getColumn(4);
        tableColumn.setCellEditor(new DateCellEditor(true));
        tableColumn.setPreferredWidth(30);
    }

    protected void initTaxonSetsTable(AbstractTableModel tableModel, final String[] columnToolTips) {
        taxonSetsTable = new JTable(tableModel) {
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
        taxonSetsTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        taxonSetsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                taxonSetsTableSelectionChanged();
            }
        });
        taxonSetsTable.doLayout();
    }

    protected void taxonSetChanged() {
        currentTaxonSet.removeAllTaxa();
        for (Taxon anIncludedTaxa : includedTaxa) {
            currentTaxonSet.addTaxon(anIncludedTaxa);
        }

        setupTaxonSetsComboBoxes();

        includedTaxaLabel.setText("" + includedTaxa.size() + " taxa included");
        excludedTaxaLabel.setText("" + excludedTaxa.size() + " taxa excluded");

        if (options.taxonSetsMono.get(currentTaxonSet) != null &&
                options.taxonSetsMono.get(currentTaxonSet) &&
                !checkCompatibility(currentTaxonSet)) {
            options.taxonSetsMono.put(currentTaxonSet, Boolean.FALSE);
        }

        frame.setAllOptions();
        frame.setDirty();
    }

    protected void treeModelsChanged() {
        Object[] modelArray = options.getPartitionTreeModels().toArray();
        TableColumn col = tableColumnModel.getColumn(3);
        col.setCellEditor(new DefaultCellEditor(new JComboBox(modelArray)));
    }

    protected void resetPanel() {
        if (!options.hasData() || options.taxonSets == null || options.taxonSets.size() < 1) {
            setCurrentTaxonSet(null);
        }
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        resetPanel();

        if (options.taxonSets == null) {
            addTaxonSetAction.setEnabled(false);
            removeTaxonSetAction.setEnabled(false);
        } else {
            addTaxonSetAction.setEnabled(options.hasData());
        }

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

    protected void taxonSetsTableSelectionChanged() {
        treeModelsChanged();

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

        private static final long serialVersionUID = 20273987098143413L;

        public void actionPerformed(ActionEvent ae) {
            taxonSetCount++;

            String newTaxonSetName = "untitled" + taxonSetCount;
            Taxa newTaxonSet = new Taxa(newTaxonSetName); // cannot use currentTaxonSet

            options.taxonSets.add(newTaxonSet);
            Collections.sort(options.taxonSets);

            options.taxonSetsMono.put(newTaxonSet, Boolean.FALSE);
            options.taxonSetsIncludeStem.put(newTaxonSet, Boolean.FALSE);
            // initialize currentTaxonSet with 1st PartitionTreeModel
            options.taxonSetsTreeModel.put(newTaxonSet, options.getPartitionTreeModels().get(0));

            setCurrentTaxonSet(newTaxonSet);

            taxonSetChanged();

            taxonSetsTableModel.fireTableDataChanged();

            int sel = options.getTaxaIndex(newTaxonSetName);
            if (sel < 0) {
                taxonSetsTable.setRowSelectionInterval(0, 0);
            } else {
                taxonSetsTable.setRowSelectionInterval(sel, sel);
            }
        }
    };

    Action removeTaxonSetAction = new AbstractAction("-") {

        private static final long serialVersionUID = 6077578872870122265L;

        public void actionPerformed(ActionEvent ae) {
            int row = taxonSetsTable.getSelectedRow();
            if (row != -1) {
                Taxa taxa = options.taxonSets.remove(row);
                options.taxonSetsMono.remove(taxa);
                options.taxonSetsIncludeStem.remove(taxa);
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

    protected void setCurrentTaxonSet(Taxa taxonSet) {

        this.currentTaxonSet = taxonSet;

        includedTaxa.clear();
        excludedTaxa.clear();

        if (currentTaxonSet != null) {
            for (int i = 0; i < taxonSet.getTaxonCount(); i++) {
                includedTaxa.add(taxonSet.getTaxon(i));
            }
            Collections.sort(includedTaxa);

            // get taxa associated to each tree
            PartitionTreeModel treeModel = options.taxonSetsTreeModel.get(currentTaxonSet);
            TaxonList alignment = options.getDataPartitions(treeModel).get(0).getTaxonList();
            Taxa taxa = new Taxa(alignment);
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

    protected void setTaxonSetTitle() {

        if (currentTaxonSet == null) {
            taxonSetEditingPanel.setBorder(BorderFactory.createTitledBorder(""));
            taxonSetEditingPanel.setEnabled(false);
        } else {
            taxonSetEditingPanel.setEnabled(true);
            taxonSetEditingPanel.setBorder(new TitledBorder(null, TAXON + ": " + currentTaxonSet.getId(), TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.ABOVE_TOP));
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
        for (Taxa taxa : options.taxonSets) {
            // AR - as these comboboxes are just intended to be handy ways of selecting taxa, I have removed
            // these requirements (it was just confusing why they weren't in the lists.
//       if (taxa != currentTaxonSet) {
//                if (isCompatible(taxa, availableTaxa)) {
            comboBox.addItem(taxa.getId()); // have to add String, otherwise it will throw Exception to cast "taxa..." into Taxa
//                }
//            }
        }
    }

    public void importTaxonSets() {

        File[] files = frame.selectImportFiles("Import Taxon Set File...", false, new FileNameExtensionFilter[]{
                new FileNameExtensionFilter("Tab-delimited text files", "txt", "tab", "dat")});

        DataTable<String[]> dataTable;

        if (files != null && files.length != 0) {
            try {
                // Load the file as a table
                dataTable = DataTable.Text.parse(new FileReader(files[0]), false, true, false, false);

            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
                return;
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe.getMessage(),
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
                return;
            } catch (Exception ex) {
                ex.printStackTrace(System.err);
                JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                        "Error reading file",
                        JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace();
                return;
            }
        } else {
            return;
        }

        String[] taxonNames = dataTable.getRowLabels();

        if (dataTable.getColumnCount() == 0) {
            String name = files[0].getName();
            int lastDot = name.lastIndexOf('.');
            if (lastDot >= 1) {
                name = name.substring(0, lastDot);
            }

            Taxa taxa = new Taxa();

            // Only one column so assume the entire set is defined in this file.
            for (String taxonName : taxonNames) {
                final int index = options.taxonList.getTaxonIndex(taxonName);
                if (index >= 0) {
                    taxa.addTaxon(options.taxonList.getTaxon(index));
                }
            }
            taxa.setId(name);
            options.taxonSets.add(taxa);
            options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
        } else {
            // assume column one is the taxon labels and column two is the set names
            String[] taxonSetNames = dataTable.getColumn(0);

            Map<String, Taxa> taxonSets = new HashMap<String, Taxa>();

            int i = 0;
            for (String taxonName : taxonNames) {
                final int index = options.taxonList.getTaxonIndex(taxonName);
                if (index >= 0) {
                    Taxa taxa = taxonSets.get(taxonSetNames[i]);
                    if (taxa == null) {
                        taxa = new Taxa();
                        taxa.setId(taxonSetNames[i]);
                        taxonSets.put(taxonSetNames[i], taxa);
                    }
                    taxa.addTaxon(options.taxonList.getTaxon(index));
                }
                i++;
            }

            for (Taxa taxa : taxonSets.values()) {
                options.taxonSets.add(taxa);
                options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
                options.taxonSetsMono.put(taxa, false);
                options.taxonSetsIncludeStem.put(taxa, false);
                options.taxonSetsHeights.put(taxa, 0.0);
            }
        }

        taxonSetsTableModel.fireTableDataChanged();
    }


    /**
     * Returns true if taxa are all found in availableTaxa
     *
     * @param taxa          a set of taxa
     * @param availableTaxa a potential superset of taxa
     * @return true if the taxa are all found in availableTaxa
     */
    protected boolean isCompatible(Taxa taxa, List<Taxon> availableTaxa) {

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);
            if (!availableTaxa.contains(taxon)) {
                return false;
            }
        }
        return true;
    }

    protected boolean checkCompatibility(Taxa taxa) {
        for (Taxa taxa2 : options.taxonSets) {
            if (taxa2 != taxa && options.taxonSetsMono.get(taxa2)
                    && options.taxonSetsTreeModel.get(taxa) == options.taxonSetsTreeModel.get(taxa2)) { // no matter if diff tree
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

    /**
     * The table on the left side of panel
     */
    protected class TaxonSetsTableModel extends AbstractTableModel {
        private static final long serialVersionUID = 3318461381525023153L;

        String[] columnNames = {"Taxon Set", "Mono?", "Stem?", "Tree", "Age"};

        public TaxonSetsTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.taxonSets.size();
        }

        public Object getValueAt(int rowIndex, int columnIndex) {
            Taxa taxonSet = options.taxonSets.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    return taxonSet.getId();
                case 1:
                    return options.taxonSetsMono.get(taxonSet);
                case 2:
                    return options.taxonSetsIncludeStem.get(taxonSet);
                case 3:
                    return options.taxonSetsTreeModel.get(taxonSet);
                case 4:
                    return options.taxonSetsHeights.get(taxonSet);
                default:
                    throw new IllegalArgumentException("unknown column, " + columnIndex);
            }
        }

        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            Taxa taxonSet = options.taxonSets.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    taxonSet.setId(aValue.toString());
                    options.renameTMRCAStatistic(taxonSet);
                    setTaxonSetTitle();
                    break;

                case 1:
                    if ((Boolean) aValue) {
                        Taxa taxa = options.taxonSets.get(rowIndex);
                        if (checkCompatibility(taxa)) {
                            options.taxonSetsMono.put(taxonSet, (Boolean) aValue);
                        }
                    } else {
                        options.taxonSetsMono.put(taxonSet, (Boolean) aValue);
                    }
                    break;

                case 2:
                    options.taxonSetsIncludeStem.put(taxonSet, (Boolean) aValue);
                    break;
                case 3:
                    options.taxonSetsTreeModel.put(taxonSet, (PartitionTreeModel) aValue);
                    this.fireTableDataChanged();
                    taxonSetsTable.getSelectionModel().setSelectionInterval(rowIndex, rowIndex);
                    break;
                case 4:
                    options.taxonSetsHeights.put(taxonSet, (Double) aValue);
                    this.fireTableDataChanged();
                    taxonSetsTable.getSelectionModel().setSelectionInterval(rowIndex, rowIndex);
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
                case 1:
                case 2:
                    return Boolean.class;
                case 3:
                    return String.class;
                case 4:
                    return Double.class;
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

    public class ImportTaxonSetsAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = 8514706149822252033L;

        public ImportTaxonSetsAction() {
            super("Import Taxon Sets");
            setToolTipText("Use this tool to import the taxon sets from a file");
        }

        public void actionPerformed(ActionEvent ae) {
            importTaxonSets();
        }
    }

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
