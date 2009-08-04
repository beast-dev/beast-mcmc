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
import dr.app.tools.TemporalRooting;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.PloidyType;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
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

//	private BeautiFrame frame = null;
	private BeautiOptions options = null;

	private boolean settingOptions = false;
	
    private final PartitionTreeModel partitionTreeModel;

    public PartitionTreeModelPanel(PartitionTreeModel partitionTreeModel, BeautiOptions options) {
    	super(12, 18);

		this.partitionTreeModel = partitionTreeModel;
		this.options = options;
		
		PanelUtils.setupComponent(ploidyTypeCombo);
		ploidyTypeCombo.addItemListener( new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
//                        setupPanel();
                    }
                }
        );

        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener( new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
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
		
		if ((options.shareSameTreePrior && options.activedSameTreePrior != null 
        		&& options.activedSameTreePrior.getNodeHeightPrior().equals(TreePrior.EXTENDED_SKYLINE))
     	   || (!options.shareSameTreePrior && partitionTreeModel != null
     			&& partitionTreeModel.getPartitionTreePrior().getNodeHeightPrior().equals(TreePrior.EXTENDED_SKYLINE))
     	   || options.isSpeciesAnalysis()) {
			
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
        
        if ((options.shareSameTreePrior && options.activedSameTreePrior != null 
        		&& options.activedSameTreePrior.getNodeHeightPrior().equals(TreePrior.EXTENDED_SKYLINE))
     	   || (!options.shareSameTreePrior && partitionTreeModel != null
     			&& partitionTreeModel.getPartitionTreePrior().getNodeHeightPrior().equals(TreePrior.EXTENDED_SKYLINE))
     	   || options.isSpeciesAnalysis()) {
        	        	
        	ploidyTypeCombo.setSelectedItem(partitionTreeModel.getPloidyType());
        }
        
        startingTreeCombo.setSelectedItem(partitionTreeModel.getStartingTreeType());
        
        if (partitionTreeModel.getUserStartingTree() != null) {
        	userTreeCombo.setSelectedItem(partitionTreeModel.getUserStartingTree().getId());
        }
        
        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
    	if (settingOptions) return;
    	
    	if ((options.shareSameTreePrior && options.activedSameTreePrior != null 
        		&& options.activedSameTreePrior.getNodeHeightPrior().equals(TreePrior.EXTENDED_SKYLINE))
     	   || (!options.shareSameTreePrior && partitionTreeModel != null
     			&& partitionTreeModel.getPartitionTreePrior().getNodeHeightPrior().equals(TreePrior.EXTENDED_SKYLINE))
     	   || options.isSpeciesAnalysis()) {
    		
        	partitionTreeModel.setPloidyType( (PloidyType) ploidyTypeCombo.getSelectedItem());
        }

    	partitionTreeModel.setStartingTreeType( (StartingTreeType) startingTreeCombo.getSelectedItem());
    	partitionTreeModel.setUserStartingTree(getSelectedUserTree(options));
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