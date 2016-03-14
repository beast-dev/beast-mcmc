/*
 * PartitionTreeModelPanel.java
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

package dr.app.beauti.treespanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PartitionTreeModelPanel extends OptionsPanel {

    private static final long serialVersionUID = 8096349200725353543L;

    private final String NO_TREE = "no tree loaded";
    private JComboBox ploidyTypeCombo = new JComboBox(PloidyType.values());

    private ButtonGroup startingTreeGroup = new ButtonGroup();
    private JRadioButton randomTreeRadio = new JRadioButton("Random starting tree");
    private JRadioButton upgmaTreeRadio = new JRadioButton("UPGMA starting tree");
    private JRadioButton userTreeRadio = new JRadioButton("User-specified starting tree");

    private JLabel userTreeLabel = new JLabel("Select user-specified tree:");
    private JComboBox userTreeCombo = new JComboBox();

    private JLabel treeFormatLabel = new JLabel("Export format for tree:");
    private JComboBox treeFormatCombo = new JComboBox(new String[] {"Newick", "XML"});

    private JLabel userTreeInfo = new JLabel("<html>" +
            "Import user-specified starting trees from <b>NEXUS</b><br>" +
            "format  data files using the 'Import Data' menu option.<br>" +
            "Trees must be rooted and strictly bifurcating (binary).</html>");

//    private JButton treeDisplayButton = new JButton("Display selected tree");
//    private JButton correctBranchLengthButton = new JButton("Correct branch lengths to get ultrametric tree");

    private RealNumberField initRootHeightField = new RealNumberField(Double.MIN_VALUE, Double.POSITIVE_INFINITY, "Init root height");

    private BeautiOptions options = null;
    private final BeautiFrame parent;
    private boolean settingOptions = false;

    PartitionTreeModel partitionTreeModel;

    public PartitionTreeModelPanel(final BeautiFrame parent, PartitionTreeModel parTreeModel, final BeautiOptions options) {
        super(12, (OSType.isMac() ? 6 : 24));

        this.partitionTreeModel = parTreeModel;
        this.options = options;
        this.parent = parent;

        PanelUtils.setupComponent(initRootHeightField);
        initRootHeightField.setColumns(10);
        initRootHeightField.setEnabled(false);

        PanelUtils.setupComponent(ploidyTypeCombo);
        ploidyTypeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                partitionTreeModel.setPloidyType((PloidyType) ploidyTypeCombo.getSelectedItem());
            }
        });
        if (options.isEBSPSharingSamePrior() || options.useStarBEAST) {
            ploidyTypeCombo.setSelectedItem(partitionTreeModel.getPloidyType());
        }

        PanelUtils.setupComponent(randomTreeRadio);
        PanelUtils.setupComponent(upgmaTreeRadio);
        PanelUtils.setupComponent(userTreeRadio);

        startingTreeGroup.add(randomTreeRadio);
        startingTreeGroup.add(upgmaTreeRadio);
        startingTreeGroup.add(userTreeRadio);

        randomTreeRadio.setSelected(partitionTreeModel.getStartingTreeType() == StartingTreeType.RANDOM);
        upgmaTreeRadio.setSelected(partitionTreeModel.getStartingTreeType() == StartingTreeType.UPGMA);
        userTreeRadio.setSelected(partitionTreeModel.getStartingTreeType() == StartingTreeType.USER);
        userTreeRadio.setEnabled(options.userTrees.size() > 0);

        boolean enabled = partitionTreeModel.getStartingTreeType() == StartingTreeType.USER;
        userTreeLabel.setEnabled(enabled);
        userTreeCombo.setEnabled(enabled);
        treeFormatLabel.setEnabled(enabled);
        treeFormatCombo.setEnabled(enabled);
        userTreeInfo.setEnabled(enabled);

        PanelUtils.setupComponent(treeFormatCombo);

        ActionListener listener = new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                if (randomTreeRadio.isSelected()) {
                    partitionTreeModel.setStartingTreeType(StartingTreeType.RANDOM);
                } else if (upgmaTreeRadio.isSelected()) {
                    partitionTreeModel.setStartingTreeType(StartingTreeType.UPGMA);
                } else if (userTreeRadio.isSelected()) {
                    partitionTreeModel.setStartingTreeType(StartingTreeType.USER);
                }
                boolean enabled = partitionTreeModel.getStartingTreeType() == StartingTreeType.USER;
                userTreeLabel.setEnabled(enabled);
                userTreeCombo.setEnabled(enabled);
                treeFormatLabel.setEnabled(enabled);
                treeFormatCombo.setEnabled(enabled);
                userTreeInfo.setEnabled(enabled);
            }
        };
        randomTreeRadio.addActionListener(listener);
        upgmaTreeRadio.addActionListener(listener);
        userTreeRadio.addActionListener(listener);

        treeFormatCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent itemEvent) {
                partitionTreeModel.setNewick(treeFormatCombo.getSelectedItem().equals("Newick"));
            }
        });

        PanelUtils.setupComponent(userTreeCombo);

        userTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                setUserSpecifiedStartingTree();
            }
        });

//        PanelUtils.setupComponent(treeDisplayButton);
//        treeDisplayButton.addActionListener(new ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                SquareTreePainter treePainter = new SquareTreePainter();
//                treePainter.setColorAttribute("color");
//                treePainter.setLineAttribute("line");
////        treePainter.setShapeAttribute("shape");
//                treePainter.setLabelAttribute("label");
//                Tree tree = getSelectedUserTree();
//                JTreeDisplay treeDisplay = new JTreeDisplay(treePainter, tree);
//
//                JTreePanel treePanel = new JTreePanel(treeDisplay);
//
//                JOptionPane optionPane = new JOptionPane(treePanel,
//                        JOptionPane.PLAIN_MESSAGE,
//                        JOptionPane.OK_OPTION,
//                        null,
//                        null,
//                        null);
//                optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));
//
//                final JDialog dialog = optionPane.createDialog(parent, "Display the selected starting tree - " + tree.getId());
//                dialog.setSize(600, 400);
//                dialog.setResizable(true);
//                dialog.setVisible(true);
//            }
//        });

//        PanelUtils.setupComponent(correctBranchLengthButton);
//        correctBranchLengthButton.addActionListener(new ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                Tree tree = getSelectedUserTree();
//                Tree.Utils.correctBranchLengthToGetUltrametricTree(tree);
//                partitionTreeModel.setUserStartingTree(tree);
//            }
//        });

        setupPanel();
    }

    private void setUserSpecifiedStartingTree() {
        if (userTreeCombo.getSelectedItem() != null && (!userTreeCombo.getSelectedItem().toString().equalsIgnoreCase(NO_TREE))) {
            Tree seleTree = getSelectedUserTree();
            if (seleTree == null || isBifurcatingTree(seleTree, seleTree.getRoot())) {
                partitionTreeModel.setUserStartingTree(seleTree);
            } else {
                JOptionPane.showMessageDialog(parent, "The selected user-specified starting tree " +
                        "is not fully bifurcating.\nBEAST requires rooted, bifurcating (binary) trees.",
                        "Illegal user-specified starting tree",
                        JOptionPane.ERROR_MESSAGE);

                userTreeCombo.setSelectedItem(NO_TREE);
                partitionTreeModel.setUserStartingTree(null);
            }
        }
    }

    public void setupPanel() {

        removeAll();

        if (options.isEBSPSharingSamePrior() || options.useStarBEAST) {
            addComponentWithLabel("Ploidy type:", ploidyTypeCombo);
        }

        if (partitionTreeModel.getDataType().getType() != DataType.MICRO_SAT) {
            addSpanningComponent(randomTreeRadio);
            addSpanningComponent(upgmaTreeRadio);
            addSpanningComponent(userTreeRadio);

            addComponents(userTreeLabel, userTreeCombo);
            userTreeCombo.removeAllItems();
            if (options.userTrees.size() < 1) {
                userTreeCombo.addItem(NO_TREE);
            } else {
                Object selectedItem = userTreeCombo.getSelectedItem();
                for (Tree tree : options.userTrees) {
                    userTreeCombo.addItem(tree.getId());
                }
                if (selectedItem != null) {
                    userTreeCombo.setSelectedItem(selectedItem);
                } else {
                    userTreeCombo.setSelectedIndex(0);
                }
            }

//            addComponent(treeDisplayButton);  // todo JTreeDisplay not work properly

            addComponents(treeFormatLabel, treeFormatCombo);
            addComponent(userTreeInfo);

        }
//		generateTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

        setOptions();
        validate();
        repaint();
    }

    public void setOptions() {

        if (partitionTreeModel == null) {
            return;
        }

//        setupPanel();

        settingOptions = true;
        initRootHeightField.setValue(partitionTreeModel.getInitialRootHeight());
//        if (options.isEBSPSharingSamePrior() || options.useStarBEAST) {
//            ploidyTypeCombo.setSelectedItem(partitionTreeModel.getPloidyType());
//        }

//        startingTreeCombo.setSelectedItem(partitionTreeModel.getStartingTreeType());
//
//        if (partitionTreeModel.getUserStartingTree() == null) {
//            userTreeCombo.setSelectedItem(NO_TREE);
//        } else {
//            userTreeCombo.setSelectedItem(partitionTreeModel.getUserStartingTree().getId());
//        }

        userTreeRadio.setEnabled(options.userTrees.size() > 0);

        settingOptions = false;

    }

    public void getOptions(BeautiOptions options) {
        if (settingOptions) return;
        // all moved to event
//    	if (options.isEBSPSharingSamePrior() || options.starBEASTOptions.isSpeciesAnalysis()) {
//
//        	partitionTreeModel.setPloidyType( (PloidyType) ploidyTypeCombo.getSelectedItem());
//        }
//
//    	partitionTreeModel.setStartingTreeType( (StartingTreeType) startingTreeCombo.getSelectedItem());
//    	partitionTreeModel.setUserStartingTree(getSelectedUserTree(options));
    }

    public boolean isBifurcatingTree(Tree tree, NodeRef node) {
        if (tree.getChildCount(node) > 2) return false;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            if (!isBifurcatingTree(tree, tree.getChild(node, i))) return false;
        }
        return true;
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

}