/*
 * DataPanel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
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
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.DummyDataType;
import dr.evolution.util.Taxa;
import jam.framework.Exportable;

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
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class DataPanel extends BeautiPanel implements Exportable {

    private static final int NAME_COLUMN = 0;
    private static final int FILENAME_COLUMN = 1;
    private static final int TAXA_COLUMN = 2;
    private static final int SITE_COLUMN = 3;
    private static final int PATTERN_COLUMN = 4;
    private static final int DATATYPE_COLUMN = 5;
    private static final int SUBSTITUTION_COLUMN = 6;
    private static final int CLOCK_COLUMN = 7;
    private static final int TREE_COLUMN = 8;

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
    CreateTreePartitionAction createTreePartitionAction = new CreateTreePartitionAction();
    CompressPatternsAction compressPatternsAction = new CompressPatternsAction();


    private final Action removePartitionAction;

    ViewPartitionAction viewPartitionAction = new ViewPartitionAction();

//    ShowAction showAction = new ShowAction();

    LinkSubstitutionModelDialog selectModelDialog = null;
    LinkClockModelDialog selectClockDialog = null;
    LinkTreeModelDialog selectTreeDialog = null;

    CreateTraitPartitionDialog selectTraitDialog = null;
    CreateTreePartitionDialog createTreeDialog = null;
    CompressPatternsDialog compressPatternsDialog = null;
    BeautiFrame frame = null;

    BeautiOptions options = null;

    public DataPanel(BeautiFrame parent, Action importDataAction, Action removePartitionAction) {

        this.frame = parent;

        this.removePartitionAction = removePartitionAction;

        dataTableModel = new DataTableModel();
        dataTable = new JTable(dataTableModel);

        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = dataTable.getColumnModel().getColumn(SUBSTITUTION_COLUMN);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
        col.setCellRenderer(comboBoxRenderer);

        col = dataTable.getColumnModel().getColumn(CLOCK_COLUMN);
        comboBoxRenderer = new ComboBoxRenderer();
        col.setCellRenderer(comboBoxRenderer);

        col = dataTable.getColumnModel().getColumn(TREE_COLUMN);
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

        toolBar1.addSeparator();

//        button = new JButton(compressPatternsAction);
//        //button.setEnabled(false);
//        compressPatternsAction.setEnabled(false);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
//        controlPanel1.add(actionPanel1);

        button = new JButton(importDataAction);
        button.setToolTipText("Import data from a file - sequences, trees or traits");
        controlPanel1.add(new JLabel("   "));
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);

        button = new JButton(removePartitionAction);
        removePartitionAction.setEnabled(false);
        button.setToolTipText("Remove selected data partition");
        controlPanel1.add(new JLabel("   "));
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);

        button = new JButton(createTraitPartitionAction);
        controlPanel1.add(new JLabel("   "));
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);

        button = new JButton(createTreePartitionAction);
        //button.setEnabled(false);
        controlPanel1.add(new JLabel("   "));
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);


        button = new JButton(viewPartitionAction);
        controlPanel1.add(new JLabel("   "));
        viewPartitionAction.setEnabled(false);
        PanelUtils.setupComponent(button);
        controlPanel1.add(button);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
        add(controlPanel1, BorderLayout.SOUTH);

    }

    private void showAlignment() {

        int[] selRows = dataTable.getSelectedRows();
        for (int row : selRows) {
            AbstractPartitionData partition = options.dataPartitions.get(row);
            Alignment alignment = null;

            if (partition instanceof PartitionData) alignment = ((PartitionData) partition).getAlignment();

            if (alignment == null) {
                JOptionPane.showMessageDialog(this, "Cannot display traits or trees data currently.\nUse the traits panel to view and edit traits.",
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
        TableColumn col = dataTable.getColumnModel().getColumn(SUBSTITUTION_COLUMN);
        col.setCellEditor(new ComboBoxCellEditor());

        col = dataTable.getColumnModel().getColumn(CLOCK_COLUMN);
        col.setCellEditor(new ComboBoxCellEditor());

        col = dataTable.getColumnModel().getColumn(TREE_COLUMN);
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

            if (column == SUBSTITUTION_COLUMN) {
                for (Object ob : options.getPartitionSubstitutionModels()) {
                    ((JComboBox) editorComponent).addItem(ob);
                }
            } else if (column == CLOCK_COLUMN) {
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

        boolean canCompress = hasSelection;

        boolean canUnlink = options.dataPartitions.size() > 1 && hasSelection;
        boolean canLink = options.dataPartitions.size() > 1 && hasSelection && selRows.length > 1;

        if (hasSelection) {
            for (int selRow : selRows) {
                if (options.dataPartitions.get(selRow) instanceof TreePartitionData) {
                    canUnlink = false;
                    canLink = false;
                    break;
                }
            }
            for (int selRow : selRows) {
                if (!(options.dataPartitions.get(selRow) instanceof PartitionData)) {
                    canCompress = false;
                    break;
                }
            }
        }

        unlinkModelsAction.setEnabled(canUnlink);
        linkModelsAction.setEnabled(canLink);

        unlinkClocksAction.setEnabled(canUnlink);
        linkClocksAction.setEnabled(canLink);

        unlinkTreesAction.setEnabled(canUnlink);
        linkTreesAction.setEnabled(canLink);
        compressPatternsAction.setEnabled(canCompress);

        removePartitionAction.setEnabled(!options.dataPartitions.isEmpty() && hasSelection);
        viewPartitionAction.setEnabled(!options.dataPartitions.isEmpty() && hasSelection);
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        modelsChanged();

        boolean traitAvailable = options.traits != null && options.traits.size() > 0;
        boolean treePartition = !options.userTrees.isEmpty();

        createTraitPartitionAction.setEnabled(traitAvailable);
        createTreePartitionAction.setEnabled(treePartition);

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

        removePartitions(partitionsToRemove);

    }

    private void removePartitions(Set<AbstractPartitionData> partitionsToRemove) {
        boolean hasIdenticalTaxa = options.hasIdenticalTaxa(); // need to check this before removing partitions

        // TODO: would probably be a good idea to check if the user wants to remove the last partition
        options.dataPartitions.removeAll(partitionsToRemove);

        if (options.dataPartitions.size() == 0) {
            // all data partitions removed so reset the taxa
            options.reset();
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

    public boolean createPartitionFromTraits(List<TraitData> selectedTraits, String partitionName, Component parent) {

        if (selectTraitDialog == null) {
            selectTraitDialog = new CreateTraitPartitionDialog(frame);
        }

        List<TraitData> traits = (selectedTraits == null ? options.traits : selectedTraits);
        boolean selectTraits = selectedTraits == null;

        String name = partitionName;

        boolean done = false;
        while (!done) {
            int result = selectTraitDialog.showDialog(traits, name, this, selectTraits);
            if (result == JOptionPane.CANCEL_OPTION) {
                return false;
            }

            if (!checkSelectedTraits(selectTraitDialog.getTraits())) {
                JOptionPane.showMessageDialog(parent, "You can't mix discrete and continuous traits when creating partition(s).", "Mixed Trait Types", JOptionPane.ERROR_MESSAGE);
                continue;
            }

            if (selectTraitDialog.getRenamePartition()) {
                name = selectTraitDialog.getName();
                if (name.trim().isEmpty()) {
                    continue;
                }
            }
            if (options.hasPartitionData(name)) {
                JOptionPane.showMessageDialog(this, "A partition with this name already exists.\nProvide a new name.",
                        "Duplicate partition name",
                        JOptionPane.ERROR_MESSAGE);
                continue;
            }

            if (selectTraits) {
                selectedTraits = selectTraitDialog.getTraits();
            }
            done = true;
        }

//        if (selectedTraits.size() > 1 && !selectTraitDialog.getForceIndependent()) {
//            // a set of traits have been passed to the function
//            int result;
//            if (alreadySelected && selectTraitDialog.getMakeCopy()) { //selectTraitDialog should not allow an empty name
//                name = selectTraitDialog.getName();
//            } else {
//                result = selectTraitDialog.showDialog(selectedTraits, partitionName, this, false);
//                name = selectTraitDialog.getName();
//
//                if (result == JOptionPane.CANCEL_OPTION) {
//                    return false;
//                }
//            }
//
//        } else {
//            name = selectedTraits.get(0).getName();
//            if (selectTraitDialog.getMakeCopy()) {
//                name = selectTraitDialog.getName();
//            }
//        }

        int minRow = -1;
        int maxRow;

        if (selectedTraits.get(0).getTraitType() == TraitData.TraitType.DISCRETE || selectTraitDialog.getForceIndependent()) {

            for (int i = 0; i < selectedTraits.size(); i++) {

                TraitData trait = selectedTraits.get(i);

                int selRow = options.createPartitionForTraits(trait.getName(), trait);

                if (i == 0) {
                    minRow = selRow;
                }
            }

            maxRow = minRow + selectedTraits.size() - 1;
        } else {
            minRow = options.createPartitionForTraits(name, selectedTraits);
            maxRow = minRow;

        }

        modelsChanged();
        dataTableModel.fireTableDataChanged();

        if (minRow != -1) {
            dataTable.getSelectionModel().setSelectionInterval(minRow, maxRow);
        }

        selectTraitDialog.reset();
        fireDataChanged();
        repaint();

        return true;
    }

    private boolean checkSelectedTraits(List<TraitData> traits) {
        int discreteCount = 0;
        int continuousCount = 0;

        for (TraitData trait : traits) {

            if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                discreteCount++;
            }
            if (trait.getTraitType() == TraitData.TraitType.CONTINUOUS) {
                continuousCount++;
            }
        }

        if (discreteCount > 0 && continuousCount > 0) {
            return false;
        }

        return true;
    }

    public void createPartitionFromTrees() {

        if (options.userTrees.isEmpty()) {
            throw new UnsupportedOperationException("No trees loaded - this button should be disabled");
        }

        if (createTreeDialog == null) {
            createTreeDialog = new CreateTreePartitionDialog(frame);
        }

        int selRow = -1;

        boolean done = false;
        while (!done) {
            int result = createTreeDialog.showDialog(options.userTrees, null, this);
            if (result != JOptionPane.CANCEL_OPTION) {
                TreeHolder tree = createTreeDialog.getTree();
                String name = tree.getName();
                if (createTreeDialog.getRenamePartition()) {
                    name = createTreeDialog.getName();
                    if (name.trim().isEmpty()) {
                        continue;
                    }
                }
                if (options.hasPartitionData(name)) {
                    JOptionPane.showMessageDialog(this, "A partition with this name already exists.\nProvide a new name.",
                            "Duplicate partition name",
                            JOptionPane.ERROR_MESSAGE);
                    continue;
                }

                selRow = options.createPartitionForTree(tree, name);
                done = true;
            } else {
                return;
            }
        }

        modelsChanged();
        dataTableModel.fireTableDataChanged();

        if (selRow != -1) {
            dataTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        frame.setAllOptions();

        fireDataChanged();
        repaint();

    }

    public void compressSitePatterns() {
        int[] selRows = dataTable.getSelectedRows();

        if (compressPatternsDialog == null) {
            compressPatternsDialog = new CompressPatternsDialog(frame);
        }

        int result = compressPatternsDialog.showDialog(this);
        if (result == JOptionPane.CANCEL_OPTION) {
            return;
        }
        SitePatterns.CompressionType compressionType = compressPatternsDialog.getCompressionType();

//        compressPatternsAction.setEnabled(false);
        setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

        for (int row : selRows) {
            PartitionData partition = (PartitionData)options.dataPartitions.get(row);
            partition.compressPatterns(compressionType);
        }
        setCursor(null);

//        SwingWorker task = new SwingWorker() {
//            @Override
//            protected Object doInBackground() throws Exception {
////                for (int row : selRows) {
////                    PartitionData partition = (PartitionData)options.dataPartitions.get(row);
////                    partition.compressPatterns();
////                }
//                Random random = new Random();
//                int progress = 0;
//                //Initialize progress property.
//                setProgress(0);
//                while (progress < 100) {
//                    //Sleep for up to one second.
//                    try {
//                        Thread.sleep(random.nextInt(1000));
//                    } catch (InterruptedException ignore) {}
//                    //Make random progress.
//                    progress += random.nextInt(10);
//                    setProgress(Math.min(progress, 100));
//                }
//                return null;
//            }
//
//            @Override
//            protected void done() {
//                setCursor(null);
//            }
//
//        };

//        ProgressMonitor progressMonitor = new ProgressMonitor(frame,
//                "Running a Long Task",
//                "", 0, 100);
//        int progress = task.getProgress();
//        String message = String.format("Completed %d%%.\n", progress);
//        progressMonitor.setNote(message);
//        progressMonitor.setProgress(progress);
//        taskOutput.append(message);

        //task.addPropertyChangeListener(this);
//        task.execute();

//        if (progressMonitor.isCanceled() || task.isDone()) {
//            progressMonitor.close();
//            Toolkit.getDefaultToolkit().beep();
//            if (progressMonitor.isCanceled()) {
//                task.cancel(true);
////                taskOutput.append("Task canceled.\n");
//            } else {
////                taskOutput.append("Task completed.\n");
//            }
////            startButton.setEnabled(true);
//        }
//
//        compressPatternsAction.setEnabled(true);


        modelsChanged();

        fireDataChanged();
        repaint();

        modelsChanged();
        dataTableModel.fireTableDataChanged();

        frame.setAllOptions();

        fireDataChanged();
        repaint();

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
            selectModelDialog = new LinkSubstitutionModelDialog(frame);
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
            selectClockDialog = new LinkClockModelDialog(frame);
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
            if (clockModel != null) {
                if (!clockModel.getName().equals(partition.getName())) {
                    clockModel = new PartitionClockModel(options, partition.getName(), clockModel);
                    partition.setPartitionClockModel(clockModel);
                }

                clockModel.setPartitionTreeModel(treeModel);
            }
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
                JOptionPane.showMessageDialog(this, errMsg, "Unsuppoted Configuration", JOptionPane.ERROR_MESSAGE);
                return;
            }
        }

        Object[] treeArray = options.getPartitionTreeModels(selectedPartitionData).toArray();

        if (selectTreeDialog == null) {
            selectTreeDialog = new LinkTreeModelDialog(frame);
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
        String[] columnNames = {"Partition Name", "File Name", "Taxa", "Sites", "Patterns", "Data Type", "Site Model", "Clock Model", "Partition Tree"};

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
            boolean isTreePartition = partition.getDataType().getType() == DataType.TREE;
            boolean isTraitPartition = partition.getTraits() != null;
            switch (col) {
                case NAME_COLUMN:
                    return partition.getName();
                case FILENAME_COLUMN:
                    return partition.getFileName();
                case TAXA_COLUMN:
                    final String nTaxa;
                    if (!isTraitPartition) {
                        nTaxa = partition.getTaxonCount() >= 0 ? String.valueOf(partition.getTaxonCount()) : "-";
                    } else {
                        nTaxa = String.valueOf(options.taxonList.getTaxonCount());
                    }
                    return nTaxa;
                case SITE_COLUMN:
                    return "" + (partition.getSiteCount() >= 0 ? partition.getSiteCount() : "-");
                case PATTERN_COLUMN:
                    return "" + (partition.getPatternCount() >= 0 ? partition.getPatternCount() : "-");
                case DATATYPE_COLUMN:
                    return partition.getDataDescription();
                case SUBSTITUTION_COLUMN:
                    return "" + (!isTreePartition ? partition.getPartitionSubstitutionModel().getName() : "-");
                case CLOCK_COLUMN:
                    return "" + (partition.getPartitionClockModel() != null ? partition.getPartitionClockModel().getName() : "-");
                case TREE_COLUMN:
                    return "" + (!isTreePartition ? partition.getPartitionTreeModel().getName() : "-");
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
            AbstractPartitionData partition = options.dataPartitions.get(row);
            switch (col) {
                case NAME_COLUMN:
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
                case SUBSTITUTION_COLUMN:
                    if (((PartitionSubstitutionModel) aValue).getDataType().equals(partition.getDataType())) {
                        partition.setPartitionSubstitutionModel((PartitionSubstitutionModel) aValue);
                    }
                    break;
                case CLOCK_COLUMN:
                    partition.setPartitionClockModel((PartitionClockModel) aValue);
                    break;
                case TREE_COLUMN:
                    partition.setPartitionTreeModel((PartitionTreeModel) aValue);
                    break;
            }
            fireDataChanged();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            AbstractPartitionData partition = options.dataPartitions.get(row);

            switch (col) {
                case NAME_COLUMN:// name
                    editable = true;
                    break;
                case SUBSTITUTION_COLUMN:// substitution model selection menu
                    editable = partition.getDataType().getType() != DataType.CONTINUOUS && partition.getDataType().getType() != DataType.TREE;
                    break;
                case CLOCK_COLUMN:// clock model selection menu
                    editable = partition.getDataType().getType() != DataType.CONTINUOUS && partition.getDataType().getType() != DataType.TREE;
                    break;
                case TREE_COLUMN:// tree selection menu
                    editable = partition.getDataType().getType() != DataType.TREE;
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
            createPartitionFromTraits(null, null, DataPanel.this);

            for (AbstractPartitionData partition : options.dataPartitions) {

                if (partition.getDataType() instanceof DummyDataType) {

                    Set<AbstractPartitionData> partitionsToRemove = new HashSet<AbstractPartitionData>();
                    partitionsToRemove.add(partition);
                    removePartitions(partitionsToRemove);

                    break; //there should only be one dummy partion

                }

            }


            selectDataPanel();
            //TODO: fix error when deleting partition

        }


    }

    public class CreateTreePartitionAction extends AbstractAction {
        public CreateTreePartitionAction() {
            super("Create partition from a tree ...");
            setToolTipText("Create a data partition from a tree. The tree should be imported first.");
        }

        public void actionPerformed(ActionEvent ae) {
            createPartitionFromTrees();

            selectDataPanel();
        }

    }

    public class CompressPatternsAction extends AbstractAction {
        public CompressPatternsAction() {
            super("Compress patterns...");
            setToolTipText("Compress an alignment into a weighted set of unique site patterns.");
        }

        public void actionPerformed(ActionEvent ae) {
            compressSitePatterns();

            selectDataPanel();
        }

    }

    private void selectDataPanel() {
        frame.tabbedPane.setSelectedComponent(this);
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