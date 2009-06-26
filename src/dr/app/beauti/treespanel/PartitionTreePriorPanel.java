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
import dr.app.beauti.options.*;
import dr.evolution.tree.Tree;

import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;

import java.awt.*;
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
	
    private JComboBox treePriorCombo = new JComboBox(EnumSet.range(TreePrior.CONSTANT, TreePrior.BIRTH_DEATH).toArray());

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

//	private BeautiFrame frame = null;
//	private BeautiOptions options = null;
    
    private final PartitionTreePrior partitionTreePrior;
    
    private boolean settingOptions = false;


    public PartitionTreePriorPanel(PartitionTreePrior partitionTreePrior) {

		this.partitionTreePrior = partitionTreePrior;        

        PanelUtils.setupComponent(treePriorCombo);
        treePriorCombo.addItemListener(new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
//	                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        KeyListener keyListener = new KeyAdapter() {
            public void keyTyped(KeyEvent ev) {
                if (ev.getKeyCode() == KeyEvent.VK_ENTER) {
//	                    fireTreePriorsChanged();
                }
            }
        };
        groupCountField.addKeyListener(keyListener);
//	        samplingProportionField.addKeyListener(keyListener);

        FocusListener focusListener = new FocusAdapter() {
            public void focusLost(FocusEvent focusEvent) {
//	                fireTreePriorsChanged();
            }
        };
        groupCountField.addFocusListener(focusListener);
//	        samplingProportionField.addFocusListener(focusListener);


       ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
//	                fireTreePriorsChanged();
            }
        };
        PanelUtils.setupComponent(parameterizationCombo);
        parameterizationCombo.addItemListener(listener);

        PanelUtils.setupComponent(bayesianSkylineCombo);
        bayesianSkylineCombo.addItemListener(listener);

        PanelUtils.setupComponent(extendedBayesianSkylineCombo);
        extendedBayesianSkylineCombo.addItemListener(listener);

        PanelUtils.setupComponent(gmrfBayesianSkyrideCombo);
        gmrfBayesianSkyrideCombo.addItemListener(listener);

        setupPanel();
	}

    private void setupPanel() {

        removeAll();
        
        addComponentWithLabel("Tree Prior:", treePriorCombo);
        
        if (treePriorCombo.getSelectedItem() == TreePrior.EXPONENTIAL ||
                treePriorCombo.getSelectedItem() == TreePrior.LOGISTIC ||
                treePriorCombo.getSelectedItem() == TreePrior.EXPANSION) {
            addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
           
        } else if (treePriorCombo.getSelectedItem() == TreePrior.SKYLINE) {
            groupCountField.setColumns(6);
            addComponentWithLabel("Number of groups:", groupCountField);
            addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
           
        } else if (treePriorCombo.getSelectedItem() == TreePrior.BIRTH_DEATH) {
        	
//            samplingProportionField.setColumns(8);
//            treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.EXTENDED_SKYLINE) {
            addComponentWithLabel("Type:", extendedBayesianSkylineCombo);
            
        } else if (treePriorCombo.getSelectedItem() == TreePrior.GMRF_SKYRIDE) {
            addComponentWithLabel("Smoothing:", gmrfBayesianSkyrideCombo);
                
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

        extendedBayesianSkylineCombo.setSelectedIndex(partitionTreePrior.isMultiLoci() ? 1 : 0);

        gmrfBayesianSkyrideCombo.setSelectedIndex(partitionTreePrior.getSkyrideSmoothing());

        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions() {
    	if (settingOptions) return;
    	
    	partitionTreePrior.setNodeHeightPrior( (TreePrior) treePriorCombo.getSelectedItem());

        if (partitionTreePrior.getNodeHeightPrior() == TreePrior.SKYLINE) {
            Integer groupCount = groupCountField.getValue();
            if (groupCount != null) {
                partitionTreePrior.setSkylineGroupCount(groupCount);
            } else {
                partitionTreePrior.setSkylineGroupCount(5);
            }
        } else if (partitionTreePrior.getNodeHeightPrior() == TreePrior.BIRTH_DEATH) {
//            Double samplingProportion = samplingProportionField.getValue();
//            if (samplingProportion != null) {
//                partitionTreePrior.birthDeathSamplingProportion = samplingProportion;
//            } else {
//                partitionTreePrior.birthDeathSamplingProportion = 1.0;
//            }
        }

        partitionTreePrior.setParameterization(parameterizationCombo.getSelectedIndex());
        partitionTreePrior.setSkylineModel(bayesianSkylineCombo.getSelectedIndex());
        partitionTreePrior.setMultiLoci(extendedBayesianSkylineCombo.getSelectedIndex() == 1);

        partitionTreePrior.setSkyrideSmoothing(gmrfBayesianSkyrideCombo.getSelectedIndex());
        // the taxon list may not exist yet... this should be set when generating...
//        partitionTreePrior.skyrideIntervalCount = partitionTreePrior.taxonList.getTaxonCount() - 1;
        
    }
    
    public void removeCertainPriorFromTreePriorCombo() {
    	treePriorCombo.removeItem(TreePrior.YULE);
    	treePriorCombo.removeItem(TreePrior.BIRTH_DEATH);
	}

	public void recoveryTreePriorCombo() {
		if (treePriorCombo.getItemCount() < EnumSet.range(TreePrior.CONSTANT, TreePrior.BIRTH_DEATH).size()) {
			treePriorCombo.addItem(TreePrior.YULE);
	    	treePriorCombo.addItem(TreePrior.BIRTH_DEATH);  
		}
	}
}