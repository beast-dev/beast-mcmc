/*
 * CharactersPanel.java
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

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;

import dr.app.gui.table.TableSorter;
import jam.table.TableRenderer;
import jam.framework.Exportable;


public class CharactersPanel extends JPanel implements Exportable {

    /**
	 *
	 */
	private static final long serialVersionUID = 1063807543195481382L;
	TreeStatFrame frame = null;
    TreeStatData treeStatData = null;
    TreeStatData.Character selectedCharacter = null;
    TreeStatData.State selectedState = null;

    JScrollPane scrollPane1 = new JScrollPane();
    JTable charactersTable = null;
    CharactersTableModel charactersTableModel = null;

    JScrollPane scrollPane2 = new JScrollPane();
    JTable statesTable = null;
    StatesTableModel statesTableModel = null;

    JScrollPane scrollPane3 = new JScrollPane();
    JTable excludedTaxaTable = null;
    TaxaTableModel excludedTaxaTableModel = null;

    JScrollPane scrollPane4 = new JScrollPane();
    JTable includedTaxaTable = null;
    TaxaTableModel includedTaxaTableModel = null;


    public CharactersPanel(TreeStatFrame frame, TreeStatData treeStatData) {

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

        // Characters
        charactersTableModel = new CharactersTableModel();
        TableSorter sorter = new TableSorter(charactersTableModel);
        charactersTable = new JTable(sorter);
        sorter.addTableModelListener(charactersTable);

        charactersTable.getColumnModel().getColumn(0).setCellRenderer(
            new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        charactersTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) { charactersTableSelectionChanged(); }
        });

           scrollPane1 = new JScrollPane(charactersTable,
                                        JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                                        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

         JPanel buttonPanel1 = createAddRemoveButtonPanel(addCharacterAction, addIcon, "Create a new character",
                                                             removeCharacterAction, removeIcon, "Remove a character",
                                                             javax.swing.BoxLayout.X_AXIS);

        // States
        statesTableModel = new StatesTableModel();
        sorter = new TableSorter(statesTableModel);
        statesTable = new JTable(sorter);
        sorter.addTableModelListener(statesTable);

        statesTable.getColumnModel().getColumn(0).setCellRenderer(
            new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        statesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) { statesTableSelectionChanged(); }
        });

           scrollPane2 = new JScrollPane(statesTable,
                                        JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                                        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

         JPanel buttonPanel2 = createAddRemoveButtonPanel(addStateAction, addIcon, "Create a new state",
                                                             removeStateAction, removeIcon, "Remove a state",
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

           scrollPane3 = new JScrollPane(excludedTaxaTable,
                                        JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                                        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

         JPanel buttonPanel3 = createAddRemoveButtonPanel(includeTaxonAction, includeIcon, "Include selected taxa in the taxon set",
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

           scrollPane4 = new JScrollPane(includedTaxaTable,
                                        JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                                        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);


        setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        c.weightx = 0.333333;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6,6,6,6);
        c.gridx = 0;
        c.gridy = 0;
        c.gridwidth = 1;
        c.gridheight = 3;
        add(scrollPane1, c);

        c.weightx = 0.333333;
        c.weighty = 0;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.LINE_START;
        c.insets = new Insets(0,6,6,6);
        c.gridx = 0;
        c.gridy = 3;
        c.gridwidth = 1;
        c.gridheight = 1;
        add(buttonPanel1, c);

        c.weightx = 0.666664;
        c.weighty = 0.333333;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6,6,6,6);
        c.gridx = 1;
        c.gridy = 0;
        c.gridwidth = 3;
        c.gridheight = 1;
        add(scrollPane2, c);

        c.weightx = 0.666664;
        c.weighty = 0;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.LINE_START;
        c.insets = new Insets(0,6,6,6);
        c.gridx = 1;
        c.gridy = 1;
        c.gridwidth = 3;
        c.gridheight = 1;
        add(buttonPanel2, c);

        c.weightx = 0.333333;
        c.weighty = 0.666664;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6,6,6,6);
        c.gridx = 1;
        c.gridy = 2;
        c.gridwidth = 1;
        c.gridheight = 1;
        add(scrollPane3, c);

        c.weightx = 0.0;
        c.weighty = 0.666664;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(0,0,0,0);
        c.gridx = 2;
        c.gridy = 2;
        c.gridwidth = 1;
        c.gridheight = 1;
        add(buttonPanel3, c);

        c.weightx = 0.333333;
        c.weighty = 0.666664;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6,6,6,6);
        c.gridx = 3;
        c.gridy = 2;
        c.gridwidth = 1;
        c.gridheight = 1;
        add(scrollPane4, c);

        setBorder(BorderFactory.createEmptyBorder(6,6,6,6));
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

    public void dataChanged() {

        addCharacterAction.setEnabled(treeStatData.allTaxa.size() > 0);

        charactersTableModel.fireTableDataChanged();
        statesTableModel.fireTableDataChanged();
        includedTaxaTableModel.fireTableDataChanged();
        excludedTaxaTableModel.fireTableDataChanged();
    }

    private void charactersTableSelectionChanged() {
        if (charactersTable.getSelectedRowCount() == 0) {
            selectedCharacter = null;
            removeCharacterAction.setEnabled(false);
            addStateAction.setEnabled(false);
        } else {
            selectedCharacter = treeStatData.characters.get(charactersTable.getSelectedRow());
            removeCharacterAction.setEnabled(true);
            addStateAction.setEnabled(true);
        }
        includedTaxaTableModel.fireTableDataChanged();
        excludedTaxaTableModel.fireTableDataChanged();
    }

    private void statesTableSelectionChanged() {
        if (statesTable.getSelectedRowCount() == 0) {
            selectedState = null;
            removeStateAction.setEnabled(false);
        } else {
            selectedCharacter = null;
            removeStateAction.setEnabled(true);
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

      Action addCharacterAction = new AbstractAction("Add") {

          /**
		 *
		 */
		private static final long serialVersionUID = -2222139644754040949L;

		public void actionPerformed(ActionEvent ae) {
            TreeStatData.Character character = new TreeStatData.Character();
            character.name = "untitled";
            character.states = new ArrayList<TreeStatData.State>();
            treeStatData.characters.add(character);
            dataChanged();

            int sel = treeStatData.characters.size() - 1;
            charactersTable.setRowSelectionInterval(sel, sel);
          }
      };

      Action removeCharacterAction = new AbstractAction("Remove") {

          /**
		 *
		 */
		private static final long serialVersionUID = -6836455052115579291L;

		public void actionPerformed(ActionEvent ae) {
            int saved = charactersTable.getSelectedRow();
            int row = charactersTable.getSelectedRow();
            if (row != -1) {
                treeStatData.characters.remove(row);
            }
            dataChanged();
            if (saved >= treeStatData.characters.size()) saved = treeStatData.characters.size() - 1;
            charactersTable.setRowSelectionInterval(saved, saved);
          }
      };

      Action addStateAction = new AbstractAction("Add") {

          /**
		 *
		 */
		private static final long serialVersionUID = -2304872597649350237L;

		public void actionPerformed(ActionEvent ae) {
          }
      };

      Action removeStateAction = new AbstractAction("Remove") {

          /**
		 *
		 */
		private static final long serialVersionUID = -6390058458520173491L;

		public void actionPerformed(ActionEvent ae) {
          }
      };

      Action includeTaxonAction = new AbstractAction("->") {

          /**
		 *
		 */
		private static final long serialVersionUID = 4577920870740752531L;

		public void actionPerformed(ActionEvent ae) {
            int saved1 = charactersTable.getSelectedRow();
            int saved2 = statesTable.getSelectedRow();
            int[] rows = excludedTaxaTable.getSelectedRows();
            ArrayList<String> exclList = new ArrayList<String>(treeStatData.allTaxa);
            exclList.removeAll(selectedState.taxa);
            for (int row : rows) {
                selectedState.taxa.add(exclList.get(row));
            }
            dataChanged();
            charactersTable.setRowSelectionInterval(saved1, saved1);
            statesTable.setRowSelectionInterval(saved2, saved2);
          }
      };

      Action excludeTaxonAction = new AbstractAction("<-") {

          /**
		 *
		 */
		private static final long serialVersionUID = 7911132810956409390L;

		public void actionPerformed(ActionEvent ae) {
            int saved1 = charactersTable.getSelectedRow();
            int saved2 = statesTable.getSelectedRow();
            int[] rows = includedTaxaTable.getSelectedRows();
            for (int i = rows.length - 1; i >= 0 ; i--) {
                selectedState.taxa.remove(rows[i]);
            }
            dataChanged();
            charactersTable.setRowSelectionInterval(saved1, saved1);
            statesTable.setRowSelectionInterval(saved2, saved2);
          }
      };

    class CharactersTableModel extends AbstractTableModel {

        /**
		 *
		 */
		private static final long serialVersionUID = -3916166866378281436L;

		public CharactersTableModel() {
        }

        public int getColumnCount() {
            return 1;
        }

        public int getRowCount() {
            return treeStatData.characters.size();
        }

        public Object getValueAt(int row, int col) {
            return treeStatData.characters.get(row).name;
        }

        public void setValueAt(Object value, int row, int col) {
            treeStatData.characters.get(row).name = (String)value;
        }

        public boolean isCellEditable(int row, int col) {
             return true;
        }

        public String getColumnName(int column) {
            return "Characters";
        }

        public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
    }

    class StatesTableModel extends AbstractTableModel {

        /**
		 *
		 */
		private static final long serialVersionUID = -1912262346368463655L;
		String[] columnNames = new String[] { "State", "Description" };

        public StatesTableModel() {
        }

        public int getColumnCount() {
            return 2;
        }

        public int getRowCount() {
            if (selectedCharacter == null) return 0;
            return selectedCharacter.states.size();
        }

        public Object getValueAt(int row, int col) {
            final TreeStatData.State state = selectedCharacter.states.get(row);
            if (col == 0) {
                return state.name;
            } else {
                return state.description;
            }
        }

        public void setValueAt(Object value, int row, int col) {
            final TreeStatData.State state = selectedCharacter.states.get(row);
            if (col == 0) {
                state.name = (String)value;
            } else {
                state.description = (String)value;
            }
        }

        public boolean isCellEditable(int row, int col) {
             return true;
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
    }

    class TaxaTableModel extends AbstractTableModel {

        /**
		 *
		 */
		private static final long serialVersionUID = 2786966293840685962L;
		boolean included;

        public TaxaTableModel(boolean included) {
            this.included = included;
        }

        public int getColumnCount() {
            return 1;
        }

        public int getRowCount() {
            if (selectedState == null || treeStatData.allTaxa == null) return 0;

            if (included) {
                return selectedState.taxa.size();
            } else {
                return treeStatData.allTaxa.size() - selectedState.taxa.size();
            }
        }

        public Object getValueAt(int row, int col) {

            if (included) {
                return selectedState.taxa.get(row);
            } else {
                ArrayList<String> exclList = new ArrayList<String>(treeStatData.allTaxa);
                exclList.removeAll(selectedState.taxa);
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
