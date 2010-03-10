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

import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.StartingTreeType;
import dr.app.beauti.options.*;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.Tree;

import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
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
	
	private RealNumberField initRootHeightField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

//	private BeautiFrame frame = null;
	private BeautiOptions options = null;

	private boolean settingOptions = false;
	
    PartitionTreeModel partitionTreeModel;

    public PartitionTreeModelPanel(PartitionTreeModel parTreeModel, BeautiOptions options) {
    	super(12, 18);

		this.partitionTreeModel = parTreeModel;
		this.options = options;
		
		PanelUtils.setupComponent(initRootHeightField);
		
		PanelUtils.setupComponent(ploidyTypeCombo);
		ploidyTypeCombo.addItemListener( new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	partitionTreeModel.setPloidyType( (PloidyType) ploidyTypeCombo.getSelectedItem());
                    }
                }
        );

        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener( new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	partitionTreeModel.setStartingTreeType( (StartingTreeType) startingTreeCombo.getSelectedItem());
                    	setupPanel();
                    }
                }
        );

        PanelUtils.setupComponent(userTreeCombo);
        userTreeCombo.addItemListener( new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	fireUserTreeChanged();
                    }
                }
        );        

		setupPanel();
	}
    
	 private void fireUserTreeChanged() {
		 partitionTreeModel.setUserStartingTree(getSelectedUserTree(options));
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
        	addComponentWithLabel("Select Tree:", userTreeCombo);        	
        } 
        
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
		
//		generateTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

		validate();
		repaint();
	}

    public void setOptions() {     

        if (partitionTreeModel == null) {
            return;
        }

        settingOptions = true;
        
        if (options.isEBSPSharingSamePrior() || options.starBEASTOptions.isSpeciesAnalysis()) {
        	        	
        	ploidyTypeCombo.setSelectedItem(partitionTreeModel.getPloidyType());
        }
        
        startingTreeCombo.setSelectedItem(partitionTreeModel.getStartingTreeType());
        
        if (partitionTreeModel.getUserStartingTree() != null) {
        	userTreeCombo.setSelectedItem(partitionTreeModel.getUserStartingTree().getId());
        }
        
        setupPanel();

        settingOptions = false;

    }

    public void getOptions(BeautiOptions options) {
    	if (settingOptions) return;
    	
//    	if (options.isEBSPSharingSamePrior() || options.isSpeciesAnalysis()) {
//    		
//        	partitionTreeModel.setPloidyType( (PloidyType) ploidyTypeCombo.getSelectedItem());
//        }

//    	partitionTreeModel.setStartingTreeType( (StartingTreeType) startingTreeCombo.getSelectedItem());
//    	partitionTreeModel.setUserStartingTree(getSelectedUserTree(options));
    }

    private Tree getSelectedUserTree(BeautiOptions options) {
        String treeId = (String) userTreeCombo.getSelectedItem();
        for (Tree tree : options.userTrees) {
            if (tree.getId().equals(treeId)) {
                return tree;
            }
        }
        return null;
    }

}