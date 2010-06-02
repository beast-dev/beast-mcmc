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
import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.StartingTreeType;
import dr.app.beauti.options.*;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.Tree;

import dr.gui.tree.JTreeDisplay;
import dr.gui.tree.JTreePanel;
import dr.gui.tree.SquareTreePainter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PartitionTreeModelPanel extends OptionsPanel {

    private static final long serialVersionUID = 8096349200725353543L;

    private JComboBox ploidyTypeCombo = new JComboBox(PloidyType.values());

    private JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
    private JComboBox userTreeCombo = new JComboBox();

    private ButtonGroup treeFormatButtonGroup = new ButtonGroup();
    private JRadioButton newickJRadioButton = new JRadioButton("Generate Newick Starting Tree");
    private JRadioButton simpleJRadioButton = new JRadioButton("Generate XML Starting Tree");

    private JButton treeDisplayButton = new JButton("Display Selected Tree");

    private RealNumberField initRootHeightField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

//	private BeautiFrame frame = null;
    private BeautiOptions options = null;

    private boolean settingOptions = false;

    PartitionTreeModel partitionTreeModel;

    public PartitionTreeModelPanel(final BeautiFrame parent, PartitionTreeModel parTreeModel, BeautiOptions options) {
        super(12, 18);

        this.partitionTreeModel = parTreeModel;
        this.options = options;

        PanelUtils.setupComponent(initRootHeightField);

        PanelUtils.setupComponent(ploidyTypeCombo);
        ploidyTypeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                partitionTreeModel.setPloidyType((PloidyType) ploidyTypeCombo.getSelectedItem());
            }
        });

        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                partitionTreeModel.setStartingTreeType((StartingTreeType) startingTreeCombo.getSelectedItem());
                setupPanel();
            }
        });

        PanelUtils.setupComponent(userTreeCombo);
        userTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                fireUserTreeChanged();
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
                        JOptionPane.OK_CANCEL_OPTION,
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

        setupPanel();
    }

    private void fireUserTreeChanged() {
        partitionTreeModel.setUserStartingTree(getSelectedUserTree());
    }

    private void setupPanel() {

        removeAll();

        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
                || options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO) {
            initRootHeightField.setValue(partitionTreeModel.getInitialRootHeight());
            initRootHeightField.setColumns(10);
            initRootHeightField.setEnabled(false);
            addComponentWithLabel("The Estimated Initial Root Height:", initRootHeightField);
        }

        if (options.isEBSPSharingSamePrior() || options.starBEASTOptions.isSpeciesAnalysis()) {

            addComponentWithLabel("Ploidy Type:", ploidyTypeCombo);
        }

        addComponentWithLabel("Starting Tree:", startingTreeCombo);

        if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
            addComponentWithLabel("Select User-specified Tree:", userTreeCombo);
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

            addComponent(treeDisplayButton);
            addComponent(newickJRadioButton);
            addComponent(simpleJRadioButton);
            newickJRadioButton.setSelected(partitionTreeModel.isNewick());
//            simpleJRadioButton.setSelected(!partitionTreeModel.isNewick());

        }

//		generateTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

        validate();
        repaint();
    }

    public void setOptions() {

        if (partitionTreeModel == null) {
            return;
        }

        setupPanel();

        settingOptions = true;

        if (options.isEBSPSharingSamePrior() || options.starBEASTOptions.isSpeciesAnalysis()) {

            ploidyTypeCombo.setSelectedItem(partitionTreeModel.getPloidyType());
        }

        startingTreeCombo.setSelectedItem(partitionTreeModel.getStartingTreeType());

        if (partitionTreeModel.getUserStartingTree() != null) {
            userTreeCombo.setSelectedItem(partitionTreeModel.getUserStartingTree().getId());
        }

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