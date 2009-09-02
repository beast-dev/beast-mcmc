/*
 * PartitionTreePriorPanel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.treespanel;

import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.util.PanelUtils;
import dr.evomodel.coalescent.VariableDemographicModel;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.*;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PartitionTreePriorPanel extends OptionsPanel {

    private static final long serialVersionUID = 5016996360264782252L;

    private JComboBox treePriorCombo = new JComboBox(EnumSet.range(TreePriorType.CONSTANT, TreePriorType.BIRTH_DEATH).toArray());

    private JComboBox parameterizationCombo = new JComboBox(new String[]{
            "Growth Rate", "Doubling Time"});
    private JComboBox parameterizationCombo1 = new JComboBox(new String[]{
            "Doubling Time"});
    private JComboBox bayesianSkylineCombo = new JComboBox(new String[]{
            "Piecewise-constant", "Piecewise-linear"});
    private WholeNumberField groupCountField = new WholeNumberField(2, Integer.MAX_VALUE);

    private JComboBox extendedBayesianSkylineCombo = new JComboBox(VariableDemographicModel.Type.values());

    JComboBox gmrfBayesianSkyrideCombo = new JComboBox(new String[]{
            "Uniform", "Time-aware"});

//    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

//	private BeautiFrame frame = null;
//	private BeautiOptions options = null;

    PartitionTreePrior partitionTreePrior;
    private final TreesPanel treesPanel;

    private boolean settingOptions = false;


    public PartitionTreePriorPanel(PartitionTreePrior parTreePrior, TreesPanel parent) {    	
    	super(12, 8);

        this.partitionTreePrior = parTreePrior;
        this.treesPanel = parent;

        PanelUtils.setupComponent(treePriorCombo);
        treePriorCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	partitionTreePrior.setNodeHeightPrior((TreePriorType) treePriorCombo.getSelectedItem());
                setupPanel();
            }
        });
        
        PanelUtils.setupComponent(parameterizationCombo);
        parameterizationCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	partitionTreePrior.setParameterization(parameterizationCombo.getSelectedIndex());                
            }
        });
        
        PanelUtils.setupComponent(parameterizationCombo1);
        parameterizationCombo1.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	partitionTreePrior.setParameterization(parameterizationCombo1.getSelectedIndex());       
            }
        });
        
        PanelUtils.setupComponent(groupCountField);
        groupCountField.addKeyListener(new java.awt.event.KeyAdapter() {
			public void keyTyped(java.awt.event.KeyEvent ev) {
				// move to here?
			}
		});
        
        
        PanelUtils.setupComponent(bayesianSkylineCombo);
        bayesianSkylineCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	partitionTreePrior.setSkylineModel(bayesianSkylineCombo.getSelectedIndex());
            }
        });
        
        PanelUtils.setupComponent(extendedBayesianSkylineCombo);
        extendedBayesianSkylineCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	partitionTreePrior.setExtendedSkylineModel(((VariableDemographicModel.Type) 
            			extendedBayesianSkylineCombo.getSelectedItem()).toString());
            }
        });

        PanelUtils.setupComponent(gmrfBayesianSkyrideCombo);
        gmrfBayesianSkyrideCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	partitionTreePrior.setSkyrideSmoothing(gmrfBayesianSkyrideCombo.getSelectedIndex());
            }
        });
        
        
        
//	        samplingProportionField.addKeyListener(keyListener);

        setupPanel();
    }

    private void setupPanel() {

        removeAll();

        addComponentWithLabel("Tree Prior:", treePriorCombo);

        if (treePriorCombo.getSelectedItem() == TreePriorType.EXPONENTIAL ||
//                treePriorCombo.getSelectedItem() == TreePriorType.LOGISTIC || //TODO Issue 93
                treePriorCombo.getSelectedItem() == TreePriorType.EXPANSION) {
            addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
            
        } else if (treePriorCombo.getSelectedItem() == TreePriorType.LOGISTIC) {//TODO Issue 93
        	addComponentWithLabel("Parameterization for growth:", parameterizationCombo1);
        	partitionTreePrior.setParameterization(parameterizationCombo1.getSelectedIndex());

        } else if (treePriorCombo.getSelectedItem() == TreePriorType.SKYLINE) {
            groupCountField.setColumns(6);
            addComponentWithLabel("Number of groups:", groupCountField);
            addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);

        } else if (treePriorCombo.getSelectedItem() == TreePriorType.BIRTH_DEATH) {

//            samplingProportionField.setColumns(8);
//            treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
        } else if (treePriorCombo.getSelectedItem() == TreePriorType.EXTENDED_SKYLINE) {
            addComponentWithLabel("Model Type:", extendedBayesianSkylineCombo);
            treesPanel.shareSameTreePriorCheck.setSelected(true);
            treesPanel.updateShareSameTreePriorChanged();
            
//            treesPanel.getFrame().setupEBSP(); TODO

        } else if (treePriorCombo.getSelectedItem() == TreePriorType.GMRF_SKYRIDE) {
            addComponentWithLabel("Smoothing:", gmrfBayesianSkyrideCombo);

        }

        if (treePriorCombo.getSelectedItem() == TreePriorType.GMRF_SKYRIDE) {
            //For GMRF, one tree prior has to be associated to one tree model. The validation is in BeastGenerator.checkOptions()
            addLabel("For GMRF, tree model/tree prior combination not implemented by BEAST yet!" +
                    "\nThe shareSameTreePrior has to be unchecked using GMRF.");

            treesPanel.shareSameTreePriorCheck.setSelected(false);
            treesPanel.shareSameTreePriorCheck.setEnabled(false);
            treesPanel.updateShareSameTreePriorChanged();
        } else {
            treesPanel.shareSameTreePriorCheck.setEnabled(true);
//            treesPanel.shareSameTreePriorCheck.setSelected(true);
//            treesPanel.updateShareSameTreePriorChanged();
        }
        
//        getOptions();
//        
//        treesPanel.treeModelPanels.get(treesPanel.currentTreeModel).setOptions();
        for (PartitionTreeModel model : treesPanel.treeModelPanels.keySet()) {
        	if (model != null) {
        		treesPanel.treeModelPanels.get(model).setOptions();
        	}
     	}
        
//        createTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

//        fireTableDataChanged();

        validate();
        repaint();
    }

    public void setOptions() {

        if (partitionTreePrior == null) {
            return;
        }

        settingOptions = true;

        treePriorCombo.setSelectedItem(partitionTreePrior.getNodeHeightPrior());

        groupCountField.setValue(partitionTreePrior.getSkylineGroupCount());
        //samplingProportionField.setValue(partitionTreePrior.birthDeathSamplingProportion);

        parameterizationCombo.setSelectedIndex(partitionTreePrior.getParameterization());
        bayesianSkylineCombo.setSelectedIndex(partitionTreePrior.getSkylineModel());

        extendedBayesianSkylineCombo.setSelectedItem(partitionTreePrior.getExtendedSkylineModel());

        gmrfBayesianSkyrideCombo.setSelectedIndex(partitionTreePrior.getSkyrideSmoothing());

//        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions() {
        if (settingOptions) return;

//        partitionTreePrior.setNodeHeightPrior((TreePriorType) treePriorCombo.getSelectedItem());

        if (partitionTreePrior.getNodeHeightPrior() == TreePriorType.SKYLINE) {
            Integer groupCount = groupCountField.getValue();
            if (groupCount != null) {
                partitionTreePrior.setSkylineGroupCount(groupCount);
            } else {
                partitionTreePrior.setSkylineGroupCount(5);
            }
        } else if (partitionTreePrior.getNodeHeightPrior() == TreePriorType.BIRTH_DEATH) {
//            Double samplingProportion = samplingProportionField.getValue();
//            if (samplingProportion != null) {
//                partitionTreePrior.birthDeathSamplingProportion = samplingProportion;
//            } else {
//                partitionTreePrior.birthDeathSamplingProportion = 1.0;
//            }
        }

//        partitionTreePrior.setParameterization(parameterizationCombo.getSelectedIndex());
//        partitionTreePrior.setSkylineModel(bayesianSkylineCombo.getSelectedIndex());
//        partitionTreePrior.setExtendedSkylineModel(((VariableDemographicModel.Type) extendedBayesianSkylineCombo.getSelectedItem()).toString());
//
//        partitionTreePrior.setSkyrideSmoothing(gmrfBayesianSkyrideCombo.getSelectedIndex());
        // the taxon list may not exist yet... this should be set when generating...
//        partitionTreePrior.skyrideIntervalCount = partitionTreePrior.taxonList.getTaxonCount() - 1;

    }

    public void removeCertainPriorFromTreePriorCombo() {
        treePriorCombo.removeItem(TreePriorType.YULE);
        treePriorCombo.removeItem(TreePriorType.BIRTH_DEATH);
    }

    public void recoveryTreePriorCombo() {
        if (treePriorCombo.getItemCount() < EnumSet.range(TreePriorType.CONSTANT, TreePriorType.BIRTH_DEATH).size()) {
            treePriorCombo.addItem(TreePriorType.YULE);
            treePriorCombo.addItem(TreePriorType.BIRTH_DEATH);
        }
    }
}