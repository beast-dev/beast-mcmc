/*
 * ModelPanel.java
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

package dr.app.beauti.treespanel;

import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.modelsPanel.CreateModelDialog;
import dr.app.beauti.*;
import dr.app.beauti.options.*;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.components.RealNumberField;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class TreesPanel extends BeautiPanel implements Exportable {

    public final static boolean DEBUG = false;

    private static final long serialVersionUID = 2778103564318492601L;

    private OptionsPanel treePriorPanel = new OptionsPanel();
    public JComboBox treePriorCombo;
    private JComboBox parameterizationCombo = new JComboBox(new String[]{
            "Growth Rate", "Doubling Time"});
    private JComboBox bayesianSkylineCombo = new JComboBox(new String[]{
            "Piecewise-constant", "Piecewise-linear"});
    private WholeNumberField groupCountField = new WholeNumberField(2, Integer.MAX_VALUE);

    JComboBox extendedBayesianSkylineCombo = new JComboBox(new String[]{
            "Single-Locus", "Multi-Loci"});

    JComboBox gmrfBayesianSkyrideCombo = new JComboBox(new String[]{
            "Uniform", "Time-aware"});

    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

    JScrollPane scrollPane = new JScrollPane();
    JTable treeTable = null;
    TreeTableModel treeTableModel = null;
    BeautiOptions options = null;

    JPanel treePanelParent;
    PartitionTree currentTree = null;

    Map<PartitionTree, PartitionTreePanel> treePanels = new HashMap<PartitionTree, PartitionTreePanel>();
    TitledBorder treeBorder;

    // General settings ////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////

    BeautiFrame frame = null;
    CreateModelDialog createModelDialog = null;
    boolean settingOptions = false;
    boolean hasAlignment = false;

    public TreesPanel(BeautiFrame parent, Action removeTreeAction) {

        super();

        this.frame = parent;

        treeTableModel = new TreeTableModel();
        treeTable = new JTable(treeTableModel);

        treeTable.getTableHeader().setReorderingAllowed(false);
        treeTable.getTableHeader().setResizingAllowed(false);
        treeTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        final TableColumnModel model = treeTable.getColumnModel();
        final TableColumn tableColumn0 = model.getColumn(0);
        tableColumn0.setCellRenderer(new TreesTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(treeTable);

        treeTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        treeTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        scrollPane = new JScrollPane(treeTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addTreeAction);
        actionPanel1.setRemoveAction(removeTreeAction);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        setCurrentTree(null);

        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                fireTreePriorsChanged();
            }
        };

        treePriorCombo = new JComboBox(EnumSet.range(TreePrior.CONSTANT, TreePrior.BIRTH_DEATH).toArray());

        PanelUtils.setupComponent(treePriorCombo);
        treePriorCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        KeyListener keyListener = new KeyAdapter() {
            public void keyTyped(KeyEvent ev) {
                if (ev.getKeyCode() == KeyEvent.VK_ENTER) {
                    fireTreePriorsChanged();
                }
            }
        };

        groupCountField.addKeyListener(keyListener);
        samplingProportionField.addKeyListener(keyListener);

        FocusListener focusListener = new FocusAdapter() {
            public void focusLost(FocusEvent focusEvent) {
                fireTreePriorsChanged();
            }
        };
        groupCountField.addFocusListener(focusListener);
//        samplingProportionField.addFocusListener(focusListener);

        PanelUtils.setupComponent(parameterizationCombo);
        parameterizationCombo.addItemListener(listener);

        PanelUtils.setupComponent(bayesianSkylineCombo);
        bayesianSkylineCombo.addItemListener(listener);

        PanelUtils.setupComponent(extendedBayesianSkylineCombo);
        extendedBayesianSkylineCombo.addItemListener(listener);

        PanelUtils.setupComponent(gmrfBayesianSkyrideCombo);
        gmrfBayesianSkyrideCombo.addItemListener(listener);

        treePriorPanel.setBorder(null);
        add(treePriorPanel, BorderLayout.CENTER);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);

        OptionsPanel panel = new OptionsPanel(10, 10);
        panel.addSeparator();
        panel.add(treePriorPanel);

        treePanelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        treePanelParent.setOpaque(false);
        treeBorder = new TitledBorder("Partition Tree");
        treePanelParent.setBorder(treeBorder);
        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, treePanelParent);
        splitPane.setDividerLocation(180);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));

        add(panel, BorderLayout.NORTH);
        add(splitPane, BorderLayout.CENTER);

        setupPanel();
    }

    private void setupPanel() {

        treePriorPanel.removeAll();

        treePriorPanel.addComponentWithLabel("Tree Prior:", treePriorCombo);
        if (treePriorCombo.getSelectedItem() == TreePrior.EXPONENTIAL ||
                treePriorCombo.getSelectedItem() == TreePrior.LOGISTIC ||
                treePriorCombo.getSelectedItem() == TreePrior.EXPANSION) {
            treePriorPanel.addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.SKYLINE) {
            groupCountField.setColumns(6);
            treePriorPanel.addComponentWithLabel("Number of groups:", groupCountField);
            treePriorPanel.addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.BIRTH_DEATH) {
//            samplingProportionField.setColumns(8);
//            treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.EXTENDED_SKYLINE) {
            treePriorPanel.addComponentWithLabel("Type:", extendedBayesianSkylineCombo);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.GMRF_SKYRIDE) {
            treePriorPanel.addComponentWithLabel("Smoothing:", gmrfBayesianSkyrideCombo);        
        }

        treePriorPanel.addSeparator();

//		generateTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

        validate();
        repaint();
    }

    public void setOptions(BeautiOptions options) {

        if (DEBUG) {
            Logger.getLogger("dr.app.beauti").info("ModelsPanel.setOptions");
        }

        this.options = options;

        settingOptions = true;

        int selRow = treeTable.getSelectedRow();
        treeTableModel.fireTableDataChanged();
        if (options.getPartitionTrees().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            treeTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        // set overall settings here
        treePriorCombo.setSelectedItem(options.nodeHeightPrior);

        groupCountField.setValue(options.skylineGroupCount);
        //samplingProportionField.setValue(options.birthDeathSamplingProportion);

        parameterizationCombo.setSelectedIndex(options.parameterization);
        bayesianSkylineCombo.setSelectedIndex(options.skylineModel);

        extendedBayesianSkylineCombo.setSelectedIndex(options.multiLoci ? 1 : 0);

        gmrfBayesianSkyrideCombo.setSelectedIndex(options.skyrideSmoothing);

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {

        // This prevents options be overwritten due to listeners calling
        // this function (indirectly through modelChanged()) whilst in the
        // middle of the setOptions() method.
        if (settingOptions) return;

        // get overall settings here
        options.nodeHeightPrior = (TreePrior) treePriorCombo.getSelectedItem();

        if (options.nodeHeightPrior == TreePrior.SKYLINE) {
            Integer groupCount = groupCountField.getValue();
            if (groupCount != null) {
                options.skylineGroupCount = groupCount;
            } else {
                options.skylineGroupCount = 5;
            }
        } else if (options.nodeHeightPrior == TreePrior.BIRTH_DEATH) {
//            Double samplingProportion = samplingProportionField.getValue();
//            if (samplingProportion != null) {
//                options.birthDeathSamplingProportion = samplingProportion;
//            } else {
//                options.birthDeathSamplingProportion = 1.0;
//            }
        }

        options.parameterization = parameterizationCombo.getSelectedIndex();
        options.skylineModel = bayesianSkylineCombo.getSelectedIndex();
        options.multiLoci = extendedBayesianSkylineCombo.getSelectedIndex() == 1;

        options.skyrideSmoothing = gmrfBayesianSkyrideCombo.getSelectedIndex();
        // the taxon list may not exist yet... this should be set when generating...
//        options.skyrideIntervalCount = options.taxonList.getTaxonCount() - 1;
    }

    private void fireTreesChanged() {
        frame.setDirty();
    }

    private void fireTreePriorsChanged() {
        frame.setDirty();
    }

    private void createPartitionTree() {
//        if (createModelDialog == null) {
//            createModelDialog = new CreateModelDialog(frame);
//        }
//
//        int result = createModelDialog.showDialog();
//        if (result != JOptionPane.CANCEL_OPTION) {
//            PartitionModel model = new PartitionModel(options, createModelDialog.getName(), createModelDialog.getDataType());
//            options.addPartitionModel(model);
//            modelTableModel.fireTableDataChanged();
//            int row = options.getPartitionModels().size() - 1;
//            modelTable.getSelectionModel().setSelectionInterval(row, row);
//        }

        fireTreesChanged();
    }

    public void removeSelection() {
        int selRow = treeTable.getSelectedRow();
        if (!isUsed(selRow)) {
            PartitionModel model = options.getPartitionModels().get(selRow);
            options.getPartitionModels().remove(model);
        }

        treeTableModel.fireTableDataChanged();
        int n = options.getPartitionModels().size();
        if (selRow >= n) {
            selRow--;
        }
        treeTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        if (n == 0) {
            setCurrentTree(null);
        }

        fireTreesChanged();
    }

    private void selectionChanged() {

        int selRow = treeTable.getSelectedRow();
        if (selRow >= 0) {
            setCurrentTree(options.getPartitionTrees().get(selRow));
            frame.modelSelectionChanged(!isUsed(selRow));
        }
    }

    /**
     * Sets the current tree that this trees panel is displaying
     *
     * @param tree the new tree to display
     */
    private void setCurrentTree(PartitionTree tree) {

        if (tree != null) {
            if (currentTree != null) treePanelParent.removeAll();

            PartitionTreePanel panel = treePanels.get(tree);
            if (panel == null) {
                panel = new PartitionTreePanel(tree);
                treePanels.put(tree, panel);
            }

            currentTree = tree;
            treePanelParent.add(panel);

            updateBorder();
        }
    }

    private void updateBorder() {

        String title = "Partition Tree: ";

        treeBorder.setTitle(title + currentTree.getName());
        repaint();
    }

    private boolean isUsed(int row) {
        PartitionModel model = options.getPartitionModels().get(row);
        for (DataPartition partition : options.dataPartitions) {
            if (partition.getPartitionModel() == model) {
                return true;
            }
        }
        return false;
    }

    public JComponent getExportableComponent() {
        return this;
    }

    class TreeTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Partition Trees"};

        public TreeTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.getPartitionTrees().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionTree tree = options.getPartitionTrees().get(row);
            switch (col) {
                case 0:
                    return tree.getName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public boolean isCellEditable(int row, int col) {
            return true;
        }

        public void setValueAt(Object value, int row, int col) {
            String name = ((String) value).trim();
            if (name.length() > 0) {
                PartitionTree tree = options.getPartitionTrees().get(row);
                tree.setName(name);
                updateBorder();
                fireTreesChanged();
            }
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

    class TreesTableCellRenderer extends TableRenderer {

        public TreesTableCellRenderer(int alignment, Insets insets) {
            super(alignment, insets);
        }

        public Component getTableCellRendererComponent(JTable aTable,
                                                       Object value,
                                                       boolean aIsSelected,
                                                       boolean aHasFocus,
                                                       int aRow, int aColumn) {

            if (value == null) return this;

            Component renderer = super.getTableCellRendererComponent(aTable,
                    value,
                    aIsSelected,
                    aHasFocus,
                    aRow, aColumn);

            if (!isUsed(aRow))
                renderer.setForeground(Color.gray);
            else
                renderer.setForeground(Color.black);
            return this;
        }

    }

    Action addTreeAction = new AbstractAction("+") {
        public void actionPerformed(ActionEvent ae) {
            createPartitionTree();
        }
    };
}