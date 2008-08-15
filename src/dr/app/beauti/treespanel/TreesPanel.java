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
import dr.app.beauti.PanelUtils;
import dr.app.beauti.options.*;
import dr.evolution.tree.Tree;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreesPanel extends JPanel {

    OptionsPanel treePriorPanel = new OptionsPanel();
    JComboBox treePriorCombo;
    JComboBox parameterizationCombo = new JComboBox(new String[]{
            "Growth Rate", "Doubling Time"});
    JComboBox bayesianSkylineCombo = new JComboBox(new String[]{
            "Piecewise-constant", "Piecewise-linear"});
    WholeNumberField groupCountField = new WholeNumberField(2, Integer.MAX_VALUE);

//    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

    JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
    JComboBox userTreeCombo = new JComboBox();

    UserTreePanel userTreePanel;

    BeautiFrame frame = null;

    public TreesPanel(BeautiFrame parent) {

        this.frame = parent;

        java.awt.event.ItemListener listener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                fireTreePriorsChanged();
            }
        };

        treePriorCombo = new JComboBox(TreePrior.values());

        PanelUtils.setupComponent(treePriorCombo);
        treePriorCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        KeyListener keyListener = new KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
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

        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        PanelUtils.setupComponent(userTreeCombo);
        userTreeCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        fireUserTreeChanged();
                    }
                }
        );

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        treePriorPanel.setBorder(null);
        add(treePriorPanel, BorderLayout.NORTH);

        userTreePanel = new UserTreePanel(parent);
        add(userTreePanel, BorderLayout.CENTER);

        setupPanel();
    }

    private void fireTreePriorsChanged() {
        if (!settingOptions) {
            frame.treePriorsChanged();
        }
    }

    private void fireUserTreeChanged() {
        if (!settingOptions) {
            frame.treePriorsChanged();
            userTreePanel.setTree((Tree)userTreeCombo.getSelectedItem());
        }
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
        }

        treePriorPanel.addSeparator();
        
        treePriorPanel.addComponentWithLabel("                          Starting Tree:", startingTreeCombo);
        if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
            treePriorPanel.addComponentWithLabel("User Tree:", userTreeCombo);
            userTreePanel.setTree((Tree)userTreeCombo.getSelectedItem());
        }

        validate();
        repaint();
    }

    private boolean settingOptions = false;

    public void setOptions(BeautiOptions options) {
        settingOptions = true;

        treePriorCombo.setSelectedItem(options.nodeHeightPrior);

        groupCountField.setValue(options.skylineGroupCount);
        //samplingProportionField.setValue(options.birthDeathSamplingProportion);

        parameterizationCombo.setSelectedIndex(options.parameterization);
        bayesianSkylineCombo.setSelectedIndex(options.skylineModel);

        startingTreeCombo.setSelectedItem(options.startingTreeType);

        userTreeCombo.removeAllItems();
        if (options.trees.size() == 0) {
            userTreeCombo.addItem("no trees loaded");
            userTreeCombo.setEnabled(false);
        } else {
            for (Tree tree : options.trees) {
                userTreeCombo.addItem(tree);
            }
            userTreeCombo.setEnabled(true);
        }

        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
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

        options.startingTreeType = (StartingTreeType)startingTreeCombo.getSelectedItem();
    }

}
