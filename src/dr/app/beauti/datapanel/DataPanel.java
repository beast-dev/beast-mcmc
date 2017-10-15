/*
 * DataPanel.java
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

package dr.app.beauti.datapanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.alignmentviewer.AlignmentViewer;
import dr.app.beauti.alignmentviewer.AminoAcidDecorator;
import dr.app.beauti.alignmentviewer.NucleotideDecorator;
import dr.app.beauti.alignmentviewer.StateCellDecorator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.TableEditorStopper;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.util.Taxa;
import jam.framework.Exportable;
import jam.panels.ActionPanel;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class DataPanel extends BeautiPanel implements Exportable {

    JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    UnlinkModelsAction unlinkModelsAction = new UnlinkModelsAction();
    LinkModelsAction linkModelsAction = new LinkModelsAction();

    UnlinkClocksAction unlinkClocksAction = new UnlinkClocksAction();
    LinkClocksAction linkClocksAction = new LinkClocksAction();

    UnlinkTreesAction unlinkTreesAction = new UnlinkTreesAction();
    LinkTreesAction linkTreesAction = new LinkTreesAction();

    CreateTraitPartitionAction createTraitPartitionAction = new CreateTraitPartitionAction();

    ViewPartitionAction viewPartitionAction = new ViewPartitionAction();

//    ShowAction showAction = new ShowAction();

    public JCheckBox useStarBEASTCheck = new JCheckBox("Use species tree ancestral reconstruction (*BEAST) Heled & Drummond 2010 ");

    SelectModelDialog selectModelDialog = null;
    SelectClockDialog selectClockDialog = null;
    SelectTreeDialog selectTreeDialog = null;

    SelectTraitDialog selectTraitDialog = null;

    BeautiFrame frame = null;

    BeautiOptions options = null;

    public DataPanel(BeautiFrame parent, Action importDataAction, Action removeDataAction/*, Action importTraitsAction*/) {

        this.frame = parent;

        dataTableModel = new DataTableModel();
        dataTable = new JTable(dataTableModel);

        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = dataTable.getColumnModel().getColumn(5);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
        col.setCellRenderer(comboBoxRenderer);

        col = dataTable.getColumnModel().getColumn(6);
        comboBoxRenderer = new ComboBoxRenderer();
        col.setCellRenderer(comboBoxRenderer);

        col = dataTable.getColumnModel().getColumn(7);
        comboBoxRenderer = new ComboBoxRenderer();
        col.setCellRenderer(comboBoxRenderer);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        dataTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    showAlignment();
                }
            }
        });
//        dataTable.setFocusable(false);

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setLayout(new BoxLayout(toolBar1, BoxLayout.X_AXIS));
        toolBar1.setBorder(BorderFactory.createEmptyBorder());

        JButton button = new JButton(unlinkModelsAction);
        unlinkModelsAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(linkModelsAction);
        linkModelsAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        toolBar1.addSeparator();

        button = new JButton(unlinkClocksAction);
        unlinkClocksAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(linkClocksAction);
        linkClocksAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        toolBar1.addSeparator();

        button = new JButton(unlinkTreesAction);
        unlinkTreesAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(linkTreesAction);
        linkTreesAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(importDataAction);
        actionPanel1.setRemoveAction(removeDataAction);
        removeDataAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        button = new JButton(viewPartitionAction);
        controlPanel1.add(new JLabel("   "));
        viewPartitionAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);

        button = new JButton(createTraitPartitionAction);
        controlPanel1.add(new JLabel("   "));
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);

        //JPanel panel1 = new JPanel(new BorderLayout());
        //panel1.setOpaque(false);
        //panel1.add(useStarBEASTCheck, BorderLayout.NORTH);
        //panel1.add(toolBar1, BorderLayout.SOUTH);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
        add(controlPanel1, BorderLayout.SOUTH);

        useStarBEASTCheck.setEnabled(false);
        useStarBEASTCheck.setToolTipText(STARBEASTOptions.CITATION);
        useStarBEASTCheck.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (!frame.setupStarBEAST(useStarBEASTCheck.isSelected())) {
                    useStarBEASTCheck.setSelected(false); // go back to unchecked
                }

                dataTableModel.fireTableDataChanged();
            }
        });

    }

    private void showAlignment() {

        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);
            Alignment alignment = null;

            if (partition instanceof PartitionData) alignment = ((PartitionData) partition).getAlignment();

            // alignment == null if partition is trait or microsat http://code.google.com/p/beast-mcmc/issues/detail?id=343
            if (alignment == null) {
                JOptionPane.showMessageDialog(this, "Cannot display traits or microsatellite data currently.\nUse the traits panel to view and edit traits.",
                        "Illegal Argument Exception",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            JFrame frame = new JFrame();
            frame.setSize(800, 600);

            AlignmentViewer viewer = new AlignmentViewer();
            if (alignment.getDataType().getType() == DataType.NUCLEOTIDES) {
                viewer.setCellDecorator(new StateCellDecorator(new NucleotideDecorator(), false));
            } else if (alignment.getDataType().getType() == DataType.AMINO_ACIDS) {
                viewer.setCellDecorator(new StateCellDecorator(new AminoAcidDecorator(), false));
            } else {
                // no colouring
            }
            viewer.setAlignmentBuffer(new BeautiAlignmentBuffer(alignment));

            JPanel panel = new JPanel(new BorderLayout());
            panel.setOpaque(false);
            panel.add(viewer, BorderLayout.CENTER);

            JPanel infoPanel = new JPanel(new BorderLayout());
            infoPanel.setOpaque(false);
            panel.add(infoPanel, BorderLayout.SOUTH);

            frame.setContentPane(panel);
            frame.setVisible(true);
        }

    }

    private void fireDataChanged() {
        // options.updatePartitionAllLinks();
        frame.setDirty();
    }

    private void modelsChanged() {
        TableColumn col = dataTable.getColumnModel().getColumn(5);
        col.setCellEditor(new ComboBoxCellEditor());

        col = dataTable.getColumnModel().getColumn(6);
        col.setCellEditor(new ComboBoxCellEditor());

        col = dataTable.getColumnModel().getColumn(7);
        col.setCellEditor(new DefaultCellEditor(new JComboBox(options.getPartitionTreeModels().toArray())));
    }

    public class ComboBoxCellEditor extends DefaultCellEditor {
        public ComboBoxCellEditor() {
            super(new JComboBox());
        }

        public Component getTableCellEditorComponent(JTable table, Object value,
                                                     boolean isSelected,
                                                     int row, int column) {

            ((JComboBox) editorComponent).removeAllItems();

            if (column == 5) {
                for (Object ob : options.getPartitionSubstitutionModels()) {
                    ((JComboBox) editorComponent).addItem(ob);
                }
            } else if (column == 6) {
                for (Object ob : options.getPartitionClockModels()) {
                    ((JComboBox) editorComponent).addItem(ob);
                }
            }

//            if (((JComboBox) editorComponent).contains(value)) // todo need validate whether value in the editorComponent

            ((JComboBox) editorComponent).setSelectedItem(value);
            delegate.setValue(value);

            return editorComponent;
        }
    }

    public void selectionChanged() {
        int[] selRows = dataTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
        frame.dataSelectionChanged(hasSelection);

        boolean canUnlink = options.dataPartitions.size() > 1 && hasSelection;
        boolean canLink = options.dataPartitions.size() > 1 && hasSelection && selRows.length > 1;

        unlinkModelsAction.setEnabled(canUnlink);
        linkModelsAction.setEnabled(canLink);

        unlinkClocksAction.setEnabled(canUnlink);
        linkClocksAction.setEnabled(canLink);

        unlinkTreesAction.setEnabled(canUnlink);
        linkTreesAction.setEnabled(canLink);

        viewPartitionAction.setEnabled(options.dataPartitions.size() > 0 && hasSelection);
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        modelsChanged();

        boolean taxaAvailable = options.taxonList != null && options.taxonList.getTaxonCount() > 0;
        boolean traitAvailable = options.traits != null && options.traits.size() > 0 && (!options.useStarBEAST);

        useStarBEASTCheck.setEnabled(taxaAvailable);
        createTraitPartitionAction.setEnabled(traitAvailable);

        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void removeSelection() {
        int[] selRows = dataTable.getSelectedRows();
        Set<AbstractPartitionData> partitionsToRemove = new HashSet<AbstractPartitionData>();
        for (int row : selRows) {
            partitionsToRemove.add(options.dataPartitions.get(row));
        }

        boolean hasIdenticalTaxa = options.hasIdenticalTaxa(); // need to check this before removing partitions

        // TODO: would probably be a good idea to check if the user wants to remove the last partition
        options.dataPartitions.removeAll(partitionsToRemove);

        if (options.dataPartitions.size() == 0) {
            // all data partitions removed so reset the taxa
            options.reset();
            useStarBEASTCheck.setSelected(false);
            frame.setupStarBEAST(false);
            frame.statusLabel.setText("");
            frame.setAllOptions();
            frame.getExportAction().setEnabled(false);
        } else if (!hasIdenticalTaxa) {
            options.updateTaxonList();
        }

        dataTableModel.fireTableDataChanged();

        fireDataChanged();
    }

    public void selectAll() {
        dataTable.selectAll();
    }

    public boolean createFromTraits(List<TraitData> traits) {
        int selRow = -1;

        if (selectTraitDialog == null) {
            selectTraitDialog = new SelectTraitDialog(frame);
        }

        if (traits==null || traits.size() == 0) {
            int result = selectTraitDialog.showDialog(options.traits, null);
            if (result != JOptionPane.CANCEL_OPTION) {
                TraitData trait = selectTraitDialog.getTrait();
                String name = trait.getName();
                if (selectTraitDialog.getMakeCopy()) {
                    name = selectTraitDialog.getName();
                }

                selRow = options.createPartitionForTraits(name, trait);
            } else {
                return false;
            }
        } else {
            if (traits.size() > 1) {
                // a set of traits have been passed to the function
                int result = selectTraitDialog.showDialog(null, null);
                if (result != JOptionPane.CANCEL_OPTION) {
                    String name = selectTraitDialog.getName();
                    selRow = options.createPartitionForTraits(name, traits);
                }  else {
                    return false;
                }
            } else {
                selRow = options.createPartitionForTraits(traits.get(0).getName(), traits);
            }
        }

        modelsChanged();
        dataTableModel.fireTableDataChanged();

        if (selRow != -1) {
            dataTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }
        fireDataChanged();
        repaint();

        return true;
    }

    public void unlinkSubstitutionModels() {
        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);

            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (!model.getName().equals(partition.getName())) {
                PartitionSubstitutionModel newModel = new PartitionSubstitutionModel(options, partition.getName(), model);
                partition.setPartitionSubstitutionModel(newModel);
            }
        }

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void linkSubstitutionModels() {
        int[] selRows = dataTable.getSelectedRows();
        List<AbstractPartitionData> selectedPartitionData = new ArrayList<AbstractPartitionData>();
        DataType dateType = null;
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);
            if (dateType == null) {
                dateType = partition.getDataType();
            } else {
                if (partition.getDataType() != dateType) {
                    JOptionPane.showMessageDialog(this, "Can only link the models for data partitions \n" +
                            "of the same data type (e.g., nucleotides)",
                            "Unable to link models",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }
            }

            if (!selectedPartitionData.contains(partition))
                selectedPartitionData.add(partition);
        }

        Object[] modelArray = options.getPartitionSubstitutionModels(selectedPartitionData).toArray();

        if (selectModelDialog == null) {
            selectModelDialog = new SelectModelDialog(frame);
        }

        int result = selectModelDialog.showDialog(modelArray);
        if (result != JOptionPane.CANCEL_OPTION) {
            PartitionSubstitutionModel model = selectModelDialog.getModel();
            if (selectModelDialog.getMakeCopy()) {
                model.setName(selectModelDialog.getName());
            }

            for (AbstractPartitionData partition : selectedPartitionData) {
                partition.setPartitionSubstitutionModel(model);
            }
        }

        if (options.getPartitionSubstitutionModels(Microsatellite.INSTANCE).size() <= 1) {
            options.shareMicroSat = true;
        }

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void unlinkClockModels() { // reuse previous PartitionTreePrior
        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);

            PartitionClockModel clockModel = partition.getPartitionClockModel();
            if (!clockModel.getName().equals(partition.getName())) {
                clockModel = new PartitionClockModel(options, partition.getName(), clockModel);
                partition.setPartitionClockModel(clockModel);
            }

            // Clock models need to refer to the same tree as the data (many to one relationship).
            clockModel.setPartitionTreeModel(partition.getPartitionTreeModel());
        }

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void linkClockModels() {
        int[] selRows = dataTable.getSelectedRows();

        List<AbstractPartitionData> selectedPartitionData = new ArrayList<AbstractPartitionData>();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);

            if (!selectedPartitionData.contains(partition)) {
                selectedPartitionData.add(partition);
            }
        }

        if (selectedPartitionData.size() > 1) {
            if (!options.hasIdenticalTaxa(selectedPartitionData)) {
                String errMsg = "To share a clock model, partitions need to have identical taxa.";
                JOptionPane.showMessageDialog(this, errMsg, "Unsuppoted Configuration", JOptionPane.ERROR_MESSAGE);
                return;
            }
        }


        Object[] modelArray = options.getPartitionClockModels(selectedPartitionData).toArray();

        if (selectClockDialog == null) {
            selectClockDialog = new SelectClockDialog(frame);
        }

        int result = selectClockDialog.showDialog(modelArray);
        if (result != JOptionPane.CANCEL_OPTION) {
            PartitionClockModel clockModel = selectClockDialog.getModel();
            if (selectClockDialog.getMakeCopy()) {
                clockModel.setName(selectClockDialog.getName());
            }

            for (AbstractPartitionData partition : selectedPartitionData) {
                partition.setPartitionClockModel(clockModel);

                // Clock models need to refer to the same tree as the data (many to one relationship).
                partition.setPartitionTreeModel(clockModel.getPartitionTreeModel());
            }
        }

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void unlinkTreeModels() { // reuse previous PartitionTreePrior
        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);

            PartitionTreeModel treeModel = partition.getPartitionTreeModel();
            if (!treeModel.getName().equals(partition.getName()) && partition.getTraits() == null) {// not a trait
                treeModel = new PartitionTreeModel(options, partition.getName(), treeModel);
                partition.setPartitionTreeModel(treeModel);
            }

            // Clock models need to refer to the same tree as the data (many to one relationship).
            // When unlinking trees, the clocks must also be unlinked. Perhaps a dialog is need
            // to say this? The parameters of the clocks could be linked in the priors panel.
            PartitionClockModel clockModel = partition.getPartitionClockModel();
            if (!clockModel.getName().equals(partition.getName())) {
                clockModel = new PartitionClockModel(options, partition.getName(), clockModel);
                partition.setPartitionClockModel(clockModel);
            }

            clockModel.setPartitionTreeModel(treeModel);
        }


        options.linkTreePriors(frame.getCurrentPartitionTreePrior());

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void linkTreeModels() { // keep previous PartitionTreePrior for reuse
        int[] selRows = dataTable.getSelectedRows();

        List<AbstractPartitionData> selectedPartitionData = new ArrayList<AbstractPartitionData>();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);

            if (!selectedPartitionData.contains(partition))
                selectedPartitionData.add(partition);
        }

        if (selectedPartitionData.size() > 1) {
            if (!options.hasIdenticalTaxa(selectedPartitionData)) {
                String errMsg = "To share a tree, partitions need to have identical taxa.";
                if (selectedPartitionData.get(0).getDataType().getType() == DataType.MICRO_SAT)
                    errMsg += "\nThe data must be all diploid or all haploid when you want to link the tree.";
                JOptionPane.showMessageDialog(this, errMsg, "Unsuppoted Configuration", JOptionPane.ERROR_MESSAGE);
                return;
            }
        }

        Object[] treeArray = options.getPartitionTreeModels(selectedPartitionData).toArray();

        if (selectTreeDialog == null) {
            selectTreeDialog = new SelectTreeDialog(frame);
        }

        int result = selectTreeDialog.showDialog(treeArray);
        if (result != JOptionPane.CANCEL_OPTION) {
            PartitionTreeModel treeModel = selectTreeDialog.getTree();
            if (selectTreeDialog.getMakeCopy()) {
                treeModel.setName(selectTreeDialog.getName());
            }
            PartitionTreePrior prior = treeModel.getPartitionTreePrior();
            options.linkTreePriors(prior);

            for (AbstractPartitionData partition : selectedPartitionData) {
                partition.setPartitionTreeModel(treeModel);

                // Clock models need to refer to the same tree as the data (many to one relationship).
                // Make sure the clock model for this partition refers to the same tree as the partition.
                PartitionClockModel clockModel = partition.getPartitionClockModel();
                clockModel.setPartitionTreeModel(treeModel);

            }

            for (Taxa taxa : options.taxonSets) { // Issue 454: all the taxon sets are deleted when link/unlink tree
                PartitionTreeModel prevModel = options.taxonSetsTreeModel.get(taxa);
                if (prevModel != treeModel) options.taxonSetsTreeModel.put(taxa, treeModel);
            }
        }

        modelsChanged();

        fireDataChanged();
        repaint();
    }

    public void unlinkAll() {
        unlinkSubstitutionModels();
        unlinkClockModels();
        unlinkTreeModels();
    }

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Partition Name", "File Name", "Taxa", "Sites", "Data Type", "Site Model", "Clock Model", "Partition Tree"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
//            return options.getPartitionDataNoSpecies().size();
            return options.dataPartitions.size();
        }

        public Object getValueAt(int row, int col) {
//            PartitionData partition = options.getPartitionDataNoSpecies().get(row);
            AbstractPartitionData partition = options.dataPartitions.get(row);
            switch (col) {
                case 0:
                    return partition.getName();
                case 1:
                    return partition.getFileName();
                case 2:
                    return "" + (partition.getTaxonCount() >= 0 ? partition.getTaxonCount() : "-");
                case 3:
                    return "" + (partition.getSiteCount() >= 0 ? partition.getSiteCount() : "-");
                case 4:
                    return partition.getDataDescription();
                case 5:
//                    return partition.getPloidyType();
//                case 6:
                    return partition.getPartitionSubstitutionModel().getName();
                case 6:
                    return "" + (partition.getPartitionClockModel() != null ? partition.getPartitionClockModel().getName() : "-");
                case 7:
                    return partition.getPartitionTreeModel().getName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
            AbstractPartitionData partition = options.dataPartitions.get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (options.hasPartitionData(name)) {
                        JOptionPane.showMessageDialog(frame, "Duplicate partition name.",
                                "Illegal Argument Exception", JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                    if (name.length() > 0) {
                        options.renamePartition(partition, name);
                    }
                    break;
                case 5:
//                    partition.setPloidyType((PloidyType) aValue);
//                    break;
//                case 6:
                    if (((PartitionSubstitutionModel) aValue).getDataType().equals(partition.getDataType())) {
                        partition.setPartitionSubstitutionModel((PartitionSubstitutionModel) aValue);
                    }
                    break;
                case 6:
                    partition.setPartitionClockModel((PartitionClockModel) aValue);
                    break;
                case 7:
                    partition.setPartitionTreeModel((PartitionTreeModel) aValue);
                    break;
            }
            fireDataChanged();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            AbstractPartitionData partition = options.dataPartitions.get(row);

            switch (col) {
                case 0:// name
                    editable = true;
                    break;
//                case 5:// ploidy type selection menu
//                    editable = true;
//                    break;
                case 5:// substitution model selection menu
                    editable = partition.getDataType().getType() != DataType.CONTINUOUS;
                    break;
                case 6:// clock model selection menu
                    editable = partition.getDataType().getType() != DataType.CONTINUOUS;
                    break;
                case 7:// tree selection menu
                    editable = true;
                    break;
                default:
                    editable = false;
            }

            return editable;
        }


        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            if (getRowCount() == 0) {
                return Object.class;
            }
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

    public class UnlinkModelsAction extends AbstractAction {
        public UnlinkModelsAction() {
            super("Unlink Subst. Models");
            setToolTipText("Use this tool to use a different substitution model for each selected data partition");
        }

        public void actionPerformed(ActionEvent ae) {
            unlinkSubstitutionModels();
        }
    }


    public class LinkModelsAction extends AbstractAction {
        public LinkModelsAction() {
            super("Link Subst. Models");
            setToolTipText("Use this tool to set all the selected partitions to the same substitution model");
        }

        public void actionPerformed(ActionEvent ae) {
            linkSubstitutionModels();
        }
    }

    public class UnlinkClocksAction extends AbstractAction {
        public UnlinkClocksAction() {
            super("Unlink Clock Models");
            setToolTipText("Use this tool to use a different clock model for each selected data partition");
        }

        public void actionPerformed(ActionEvent ae) {
            unlinkClockModels();
        }
    }


    public class LinkClocksAction extends AbstractAction {
        public LinkClocksAction() {
            super("Link Clock Models");
            setToolTipText("Use this tool to set all the selected partitions to the same clock model");
        }

        public void actionPerformed(ActionEvent ae) {
            linkClockModels();
        }
    }

    public class UnlinkTreesAction extends AbstractAction {
        public UnlinkTreesAction() {
            super("Unlink Trees");
            setToolTipText("Use this tool to use a different tree for each selected data partition");
        }

        public void actionPerformed(ActionEvent ae) {
            unlinkTreeModels();
        }
    }


    public class LinkTreesAction extends AbstractAction {
        public LinkTreesAction() {
            super("Link Trees");
            setToolTipText("Use this tool to set all the selected partitions to the same tree");
        }

        public void actionPerformed(ActionEvent ae) {
            linkTreeModels();
        }
    }

    public class ViewPartitionAction extends AbstractAction {
        public ViewPartitionAction() {
            super("View Partition ...");
            setToolTipText("View alignment in a window.");
        }

        public void actionPerformed(ActionEvent ae) {
            showAlignment();
        }
    }

    public class CreateTraitPartitionAction extends AbstractAction {
        public CreateTraitPartitionAction() {
            super("Create partition from trait ...");
            setToolTipText("Create a data partition from a trait. Traits can be defined in the Traits panel.");
        }

        public void actionPerformed(ActionEvent ae) {
            createFromTraits(null);
        }
    }

//    public class ShowAction extends AbstractAction {
//        public ShowAction() {
//            super("Show");
//            setToolTipText("Display the selected alignments");
//        }
//
//        public void actionPerformed(ActionEvent ae) {
//            showAlignment();
//        }
//    }

//    public class UnlinkAllAction extends AbstractAction {
//        public UnlinkAllAction() {
//            super("Unlink All");
//            setToolTipText("Use this tool to use a different substitution model, different clock model and different tree for each selected data partition");
//        }
//
//        public void actionPerformed(ActionEvent ae) {
//            unlinkAll();
//        }
//    }
//
//    public class LinkAllAction extends AbstractAction {
//        public LinkAllAction() {
//            super("Link All");
//            setToolTipText("Use this tool to set all the selected partitions to the same substitution, clock model and tree");
//        }
//
//        public void actionPerformed(ActionEvent ae) {
//            linkAll();
//        }
//    }
}