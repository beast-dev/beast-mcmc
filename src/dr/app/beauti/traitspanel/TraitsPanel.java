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

package dr.app.beauti.traitspanel;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.*;
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
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class TraitsPanel extends BeautiPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = 5283922195494563924L;

    private final JTable traitsTable;
    private final TraitsTableModel traitsTableModel;

    private final JTable dataTable;
    private final DataTableModel dataTableModel;

    private final ClearTraitAction clearTraitAction = new ClearTraitAction();
    private final GuessLocationsAction guessLocationsAction = new GuessLocationsAction();

    private final BeautiFrame frame;

    private BeautiOptions options = null;

    private String selectedTrait = null;

    private CreateTraitDialog createTraitDialog = null;
    private GuessTraitDialog guessTraitDialog = null;

    public TraitsPanel(BeautiFrame parent, Action importTraitsAction) {

        this.frame = parent;

        traitsTableModel = new TraitsTableModel();
        TableSorter sorter = new TableSorter(traitsTableModel);
        traitsTable = new JTable(sorter);
        sorter.setTableHeader(traitsTable.getTableHeader());

        traitsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        traitsTable.getColumnModel().getColumn(0).setPreferredWidth(80);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(traitsTable);

        traitsTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        traitsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                traitSelectionChanged();
            }
        });

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
                traitSelectionChanged();
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
        JButton button = new JButton(guessLocationsAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);
        button = new JButton(clearTraitAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);
//        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addTraitAction);
        actionPanel1.setRemoveAction(removeTraitAction);

        removeTraitAction.setEnabled(false);

        JToolBar toolBar2 = new JToolBar();
        toolBar2.setFloatable(false);
        toolBar2.setOpaque(false);
        toolBar2.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        button = new JButton(importTraitsAction);
        PanelUtils.setupComponent(button);
        toolBar2.add(button);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(toolBar2, BorderLayout.NORTH);
        panel1.add(scrollPane1, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(toolBar1, BorderLayout.NORTH);
        panel2.add(scrollPane2, BorderLayout.CENTER);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                panel1, panel2);
        splitPane.setDividerLocation(240);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.CENTER);
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        traitsTableModel.fireTableDataChanged();
        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
//        options.datesUnits = unitsCombo.getSelectedIndex();
//        options.datesDirection = directionCombo.getSelectedIndex();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    private void fireTraitsChanged() {
        frame.setDirty();
    }

    private void traitSelectionChanged() {
        int selRow = traitsTable.getSelectedRow();
        if (selRow >= 0) {
            selectedTrait = options.selecetedTraits.get(selRow);
            removeTraitAction.setEnabled(true);
        } else {
            selectedTrait = null;
            removeTraitAction.setEnabled(false);
        }
        dataTableModel.fireTableDataChanged();
    }


    public void clearTrait() {
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

    public void guessTrait() {
    	if( options.taxonList != null ) { // validation of check empty taxonList
	    	TraitGuesser guesser = options.traitGuesser;
	    	
	        if (guessTraitDialog == null) {
	        	guessTraitDialog = new GuessTraitDialog(frame, guesser);
	        }
	
	        int result = guessTraitDialog.showDialog();
	
	        if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
	            return;
	        }
	
	        guesser.guessTrait = true; // ?? no use?
	        guessTraitDialog.setupGuesser();
		        
	        try {
	        	guesser.guessTrait(options);	               
	        	        
				switch (guesser.traitAnalysisType) {
					case TRAIT_SPECIES :
						addTrait(TraitGuesser.Traits.TRAIT_SPECIES.toString(), guesser.traitType);
						break;
					default:
						throw new IllegalArgumentException("unknown trait selected");       	
		        }
	        } catch (IllegalArgumentException iae) {
                JOptionPane.showMessageDialog(this, iae.getMessage(), "Unable to guess trait value", JOptionPane.ERROR_MESSAGE);                
            } 
	
	        dataTableModel.fireTableDataChanged();
	        frame.setDirty();
    	} else {
    		JOptionPane.showMessageDialog(this, "No taxa loaded yet, please import Alignment file!",
                  "No taxa loaded", JOptionPane.ERROR_MESSAGE);
    	}
    }

    private void addTrait() {
        if (createTraitDialog == null) {
            createTraitDialog = new CreateTraitDialog(frame);
        }

        int result = createTraitDialog.showDialog(options);
        if (result != JOptionPane.CANCEL_OPTION) {
            String name = createTraitDialog.getName();
            TraitGuesser.TraitType type = createTraitDialog.getType();

            if (!options.selecetedTraits.contains(name)) {
                // The createTraitDialog will have already checked if the
                // user is overwriting an existing trait
                options.selecetedTraits.add(name);
            }
            options.traitTypes.put(name, type);
            traitsTableModel.fireTableDataChanged();
            int row = options.selecetedTraits.size() - 1;
            traitsTable.getSelectionModel().setSelectionInterval(row, row);
        }

        fireTraitsChanged();

    }
    
    private void addTrait(String traitName, TraitGuesser.TraitType traitType) {
    	
        if (!options.selecetedTraits.contains(traitName)) {
        	options.selecetedTraits.add(traitName);
        }
	    options.traitTypes.put(traitName, traitType);
	    
	    traitsTableModel.fireTableDataChanged();
	    int row = options.selecetedTraits.size() - 1;
	    traitsTable.getSelectionModel().setSelectionInterval(row, row);    	
    }

    private void removeTrait() {
        int selRow = traitsTable.getSelectedRow();
        String trait = options.selecetedTraits.get(selRow);
        options.selecetedTraits.remove(trait);
        options.traitTypes.remove(trait);
        traitsTableModel.fireTableDataChanged();
        fireTraitsChanged();
    }

    public class ClearTraitAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = -7281309694753868635L;

        public ClearTraitAction() {
            super("Clear trait values");
            setToolTipText("Use this tool to remove trait values from each taxon");
        }

        public void actionPerformed(ActionEvent ae) {
            clearTrait();
        }
    }

    public class GuessLocationsAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = 8514706149822252033L;

        public GuessLocationsAction() {
            super("Guess trait values");
            setToolTipText("Use this tool to guess the trait values from the taxon labels");
        }

        public void actionPerformed(ActionEvent ae) {        	
        	guessTrait();        	
        }
    }

    AbstractAction addTraitAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            addTrait();
        }
    };


    AbstractAction removeTraitAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            removeTrait();
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
            if (options.selecetedTraits == null) return 0;

            return options.selecetedTraits.size();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.selecetedTraits.get(row);
                case 1:
                    return options.traitTypes.get(options.selecetedTraits.get(row));
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            if (col == 0) {
                options.selecetedTraits.set(row, aValue.toString());
            } else if (col == 1) {
                options.traitTypes.put(options.selecetedTraits.get(row), (TraitGuesser.TraitType)aValue);
            }
        }

        public boolean isCellEditable(int row, int col) {
//            if (col == 0) return true;
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
                    Object value = null;
                    if (selectedTrait != null) {
                        value = options.taxonList.getTaxon(row).getAttribute(selectedTrait);
                    }
                    if (value != null) {
                        return value;
                    } else {
                        return "-";
                    }
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            if (col == 1) {
                Location location = options.taxonList.getTaxon(row).getLocation();
                if (location != null) {
                    options.taxonList.getTaxon(row).setLocation(location);
                }
            }
        }

        public boolean isCellEditable(int row, int col) {
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