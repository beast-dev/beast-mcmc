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
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockModelGroup;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.types.FixRateType;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.util.TextUtil;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.tree.JTreeDisplay;
import dr.app.gui.tree.JTreePanel;
import dr.app.gui.tree.SquareTreePainter;
import dr.app.util.OSType;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
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

    private JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
    private JComboBox userTreeCombo = new JComboBox();

    private ButtonGroup treeFormatButtonGroup = new ButtonGroup();
    private JRadioButton newickJRadioButton = new JRadioButton("Generate Newick Starting Tree");
    private JRadioButton simpleJRadioButton = new JRadioButton("Generate XML Starting Tree");

    private JButton treeDisplayButton = new JButton("Display Selected Tree");
    private JButton correctBranchLengthButton = new JButton("Correct Branch Length to Get Ultrametric Tree");
    private JButton exampleButton = new JButton("Introduce how to load starting tree by user");

    private RealNumberField initRootHeightField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

    private BeautiOptions options = null;

    private boolean settingOptions = false;

    PartitionTreeModel partitionTreeModel;

    public PartitionTreeModelPanel(final BeautiFrame parent, PartitionTreeModel parTreeModel, BeautiOptions options) {
        super(12, (OSType.isMac() ? 6 : 24));

        this.partitionTreeModel = parTreeModel;
        this.options = options;

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
        
        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                partitionTreeModel.setStartingTreeType((StartingTreeType) startingTreeCombo.getSelectedItem());
                setupPanel();
            }
        });
        startingTreeCombo.setSelectedItem(partitionTreeModel.getStartingTreeType());

        PanelUtils.setupComponent(userTreeCombo);
        userTreeCombo.addItem(NO_TREE);
        if (options.userTrees.size() < 1) {
            userTreeCombo.setEnabled(false);
        } else {
            for (Tree tree : options.userTrees) {
                userTreeCombo.addItem(tree.getId());
            }
            userTreeCombo.setSelectedIndex(1);
            userTreeCombo.setEnabled(true);
        }
        userTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
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
        });

        PanelUtils.setupComponent(treeDisplayButton);
        treeDisplayButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                SquareTreePainter treePainter = new SquareTreePainter();
                treePainter.setColorAttribute("color");
                treePainter.setLineAttribute("line");
//        treePainter.setShapeAttribute("shape");
                treePainter.setLabelAttribute("label");
                Tree tree = getSelectedUserTree();
                JTreeDisplay treeDisplay = new JTreeDisplay(treePainter, tree);

                JTreePanel treePanel = new JTreePanel(treeDisplay);

                JOptionPane optionPane = new JOptionPane(treePanel,
                        JOptionPane.PLAIN_MESSAGE,
                        JOptionPane.OK_OPTION,
                        null,
                        null,
                        null);
                optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

                final JDialog dialog = optionPane.createDialog(parent, "Display the selected starting tree - " + tree.getId());
                dialog.setSize(600, 400);
                dialog.setResizable(true);
                dialog.setVisible(true);
            }
        });

        treeFormatButtonGroup.add(newickJRadioButton);
        treeFormatButtonGroup.add(simpleJRadioButton);
        ItemListener itemListener = new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                partitionTreeModel.setNewick(newickJRadioButton.isSelected());
            }
        };
        newickJRadioButton.addItemListener(itemListener);
        simpleJRadioButton.addItemListener(itemListener);

        PanelUtils.setupComponent(correctBranchLengthButton);
        correctBranchLengthButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
//                Tree tree = getSelectedUserTree();
//                Tree.Utils.correctBranchLengthToGetUltrametricTree(tree);
//                partitionTreeModel.setUserStartingTree(tree);
            }
        });

        exampleButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                JScrollPane scrollPane = TextUtil.createHTMLScrollPane(PartitionTreeModel.USER_SPEC_TREE_FORMAT, new Dimension(400, 200));

                JOptionPane.showMessageDialog(parent, scrollPane,
                        "Introduction of loading starting tree by user",
                        JOptionPane.PLAIN_MESSAGE);
            }
        });

        setupPanel();
    }

    public void setupPanel() {

        removeAll();

        ClockModelGroup group = null;
        if (options.getAllPartitionData(partitionTreeModel).size() > 0)
            group = options.getAllPartitionData(partitionTreeModel).get(0).getPartitionClockModel().getClockModelGroup();

        if (group != null && (group.getRateTypeOption() == FixRateType.FIX_MEAN
                || group.getRateTypeOption() == FixRateType.RELATIVE_TO)) {
            addComponentWithLabel("The Estimated Initial Root Height:", initRootHeightField);
        }

        if (options.isEBSPSharingSamePrior() || options.useStarBEAST) {
            addComponentWithLabel("Ploidy Type:", ploidyTypeCombo);
        }

        if (partitionTreeModel.getDataType().getType() != Microsatellite.INSTANCE.getType())
            addComponentWithLabel("Starting Tree:", startingTreeCombo);

        if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
            addComponentWithLabel("Select User-specified Tree:", userTreeCombo);
//            userTreeCombo.removeAllItems();

//            addComponent(treeDisplayButton);  // todo JTreeDisplay not work properly
            addComponent(newickJRadioButton);
            addComponent(simpleJRadioButton);
            addComponent(exampleButton);
            newickJRadioButton.setSelected(partitionTreeModel.isNewick());
//            simpleJRadioButton.setSelected(!partitionTreeModel.isNewick());

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