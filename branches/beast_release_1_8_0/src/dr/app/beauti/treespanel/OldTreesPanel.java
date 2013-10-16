/*
 * PriorsPanel.java
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

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.WholeNumberField;
import dr.app.gui.table.TableEditorStopper;
import dr.app.pathogen.TemporalRooting;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 * @deprecated
 */
public class OldTreesPanel extends BeautiPanel {

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

//    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

    private JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
    private JComboBox userTreeCombo = new JComboBox();
    private JButton button;

    private CreateTreeAction createTreeAction = new CreateTreeAction();
    private TreeDisplayPanel treeDisplayPanel;

    private BeautiFrame frame = null;
    private BeautiOptions options = null;

    private JScrollPane scrollPane = new JScrollPane();
    private JTable treesTable = null;
    private TreesTableModel treesTableModel = null;

    private GenerateTreeDialog generateTreeDialog = null;

    public OldTreesPanel(BeautiFrame parent) {

        this.frame = parent;

        treesTableModel = new TreesTableModel();
        treesTable = new JTable(treesTableModel);

        treesTable.getTableHeader().setReorderingAllowed(false);
        treesTable.getTableHeader().setResizingAllowed(false);
//        treesTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        final TableColumnModel model = treesTable.getColumnModel();

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(treesTable);

        treesTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        treesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        scrollPane = new JScrollPane(treesTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);
        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        button = new JButton(createTreeAction);
        createTreeAction.setEnabled(true);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

//        button = new JButton(linkModelAction);
//        linkModelAction.setEnabled(false);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);

//        ActionPanel actionPanel1 = new ActionPanel(false);
//        actionPanel1.setAddAction(addTreeAction);
//        actionPanel1.setRemoveAction(removeTreeAction);
//
//        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
//        controlPanel1.setOpaque(false);
//        controlPanel1.add(actionPanel1);

        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                fireTreePriorsChanged();
            }
        };

        treePriorCombo = new JComboBox(EnumSet.range(TreePriorType.CONSTANT, TreePriorType.BIRTH_DEATH).toArray());

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
//        samplingProportionField.addKeyListener(keyListener);

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

        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        PanelUtils.setupComponent(userTreeCombo);
        userTreeCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        fireTreePriorsChanged();
                    }
                }
        );

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
//        panel1.add(controlPanel1, BorderLayout.SOUTH);

        treeDisplayPanel = new TreeDisplayPanel(parent);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, treeDisplayPanel);
        splitPane.setDividerLocation(180);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(toolBar1, BorderLayout.NORTH);
        panel2.add(splitPane, BorderLayout.CENTER);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));

        treePriorPanel.setBorder(null);
        add(treePriorPanel, BorderLayout.NORTH);

        add(panel2, BorderLayout.CENTER);

        setupPanel();
    }

    private void fireTreePriorsChanged() {
        if (!settingOptions) {
            frame.setDirty();
        }
    }

    private void selectionChanged() {
        int selRow = treesTable.getSelectedRow();
        if (selRow >= 0) {
            treeDisplayPanel.setTree(options.userTrees.get(selRow));
        } else {
            treeDisplayPanel.setTree(null);
        }
    }

    private void createTree() {
        if (generateTreeDialog == null) {
            generateTreeDialog = new GenerateTreeDialog(frame);
        }

        int result = generateTreeDialog.showDialog(options);
        if (result != JOptionPane.CANCEL_OPTION) {
            GenerateTreeDialog.MethodTypes methodType = generateTreeDialog.getMethodType();
            PartitionData partition = generateTreeDialog.getDataPartition();

            Patterns patterns = new Patterns(partition.getAlignment());
            DistanceMatrix distances = new F84DistanceMatrix(patterns);
            Tree tree;
            TemporalRooting temporalRooting;

            switch (methodType) {
                case NJ:
                    tree = new NeighborJoiningTree(distances);
                    temporalRooting = new TemporalRooting(tree);
                    tree = temporalRooting.findRoot(tree, TemporalRooting.RootingFunction.CORRELATION);
                    break;
                case UPGMA:
                    tree = new UPGMATree(distances);
                    temporalRooting = new TemporalRooting(tree);
                    break;
                default:
                    throw new IllegalArgumentException("unknown method type");
            }

            tree.setId(generateTreeDialog.getName());
            options.userTrees.add(tree);
            treesTableModel.fireTableDataChanged();
            int row = options.userTrees.size() - 1;
            treesTable.getSelectionModel().setSelectionInterval(row, row);
        }

        fireTreePriorsChanged();

    }

    private void setupPanel() {

        treePriorPanel.removeAll();

        treePriorPanel.addComponentWithLabel("Tree Prior:", treePriorCombo);
        if (treePriorCombo.getSelectedItem() == TreePriorType.EXPONENTIAL ||
                treePriorCombo.getSelectedItem() == TreePriorType.LOGISTIC ||
                treePriorCombo.getSelectedItem() == TreePriorType.EXPANSION) {
            treePriorPanel.addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
            button.setEnabled(true);
        } else if (treePriorCombo.getSelectedItem() == TreePriorType.SKYLINE) {
            groupCountField.setColumns(6);
            treePriorPanel.addComponentWithLabel("Number of groups:", groupCountField);
            treePriorPanel.addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
            button.setEnabled(true);
        } else if (treePriorCombo.getSelectedItem() == TreePriorType.BIRTH_DEATH) {
        	button.setEnabled(true);
//            samplingProportionField.setColumns(8);
//            treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
        } else if (treePriorCombo.getSelectedItem() == TreePriorType.EXTENDED_SKYLINE) {
            treePriorPanel.addComponentWithLabel("Type:", extendedBayesianSkylineCombo);
            button.setEnabled(true);
        } else if (treePriorCombo.getSelectedItem() == TreePriorType.GMRF_SKYRIDE) {
            treePriorPanel.addComponentWithLabel("Smoothing:", gmrfBayesianSkyrideCombo);
            button.setEnabled(true);
        } else {
        	button.setEnabled(true);
        }

        treePriorPanel.addSeparator();

        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        panel.setOpaque(false);
        panel.add(startingTreeCombo);
        if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
            panel.add(new JLabel("  Select Tree:"));
            panel.add(userTreeCombo);
        }
        treePriorPanel.addComponentWithLabel("                          Starting Tree:", panel);

        treePriorPanel.addSeparator();

        createTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

        treesTableModel.fireTableDataChanged();

        validate();
        repaint();
    }

    private boolean settingOptions = false;

    public void setOptions(BeautiOptions options) {
        this.options = options;

        settingOptions = true;

//        treePriorCombo.setSelectedItem(options.nodeHeightPrior);
//
//        groupCountField.setValue(options.skylineGroupCount);
//        //samplingProportionField.setValue(options.birthDeathSamplingProportion);
//
//        parameterizationCombo.setSelectedIndex(options.parameterization);
//        bayesianSkylineCombo.setSelectedIndex(options.skylineModel);
//
//        extendedBayesianSkylineCombo.setSelectedIndex(options.multiLoci ? 1 : 0);
//
//        gmrfBayesianSkyrideCombo.setSelectedIndex(options.skyrideSmoothing);
//
//        startingTreeCombo.setSelectedItem(options.startingTreeType);

        userTreeCombo.removeAllItems();
        if (options.userTrees.size() == 0) {
            userTreeCombo.addItem("no trees loaded");
            userTreeCombo.setEnabled(false);
        } else {
            for (Tree tree : options.userTrees) {
                userTreeCombo.addItem(tree.getId());
            }
            userTreeCombo.setEnabled(true);
        }

        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
//        options.nodeHeightPrior = (TreePriorType) treePriorCombo.getSelectedItem();
//
//        if (options.nodeHeightPrior == TreePriorType.SKYLINE) {
//            Integer groupCount = groupCountField.getValue();
//            if (groupCount != null) {
//                options.skylineGroupCount = groupCount;
//            } else {
//                options.skylineGroupCount = 5;
//            }
//        } else if (options.nodeHeightPrior == TreePriorType.BIRTH_DEATH) {
////            Double samplingProportion = samplingProportionField.getValue();
////            if (samplingProportion != null) {
////                options.birthDeathSamplingProportion = samplingProportion;
////            } else {
////                options.birthDeathSamplingProportion = 1.0;
////            }
//        }
//
//        options.parameterization = parameterizationCombo.getSelectedIndex();
//        options.skylineModel = bayesianSkylineCombo.getSelectedIndex();
//        options.multiLoci = extendedBayesianSkylineCombo.getSelectedIndex() == 1;
//
//        options.skyrideSmoothing = gmrfBayesianSkyrideCombo.getSelectedIndex();
//        // the taxon list may not exist yet... this should be set when generating...
////        options.skyrideIntervalCount = taxonList.taxonList.getTaxonCount() - 1;
//
//        options.startingTreeType = (StartingTreeType) startingTreeCombo.getSelectedItem();
//        options.userStartingTree = getSelectedUserTree();
    }

    private Tree getSelectedUserTree() {
        String treeId = (String) userTreeCombo.getSelectedItem();
        for (Tree tree : options.userTrees) {
            if (tree.getId().equals(treeId)) {
                return tree;
            }
        }
        return null;
    }

    public JComponent getExportableComponent() {
        return treeDisplayPanel;
    }

    class TreesTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Trees"};

        public TreesTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.userTrees.size();
        }

        public Object getValueAt(int row, int col) {
            Tree tree = options.userTrees.get(row);
            switch (col) {
                case 0:
                    return tree.getId();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
            Tree tree = options.userTrees.get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        tree.setId(name);
                    }
                    break;
            }
            fireTreePriorsChanged();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            switch (col) {
                case 0:// name
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

    public class CreateTreeAction extends AbstractAction {
        public CreateTreeAction() {
            super("Create Tree");
            setToolTipText("Create a NJ or UPGMA tree using a data partition");
        }

        public void actionPerformed(ActionEvent ae) {
            createTree();
        }
    }


}