/*
 * SpeciesTreesPanel.java
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
import dr.app.beauti.options.TreePrior;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.*;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: SpeciesTreesPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class SpeciesTreesPanel extends BeautiPanel {

    private OptionsPanel optionsPanel = new OptionsPanel(30, 50);
    public JComboBox treePriorCombo;

    private BeautiFrame frame = null;
    private BeautiOptions options = null;

    public SpeciesTreesPanel(BeautiFrame parent) {
    	
        this.frame = parent;

        setOpaque(false);
        optionsPanel.setOpaque(false);
        setLayout(new BorderLayout());
    	setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
    	

    	treePriorCombo = new JComboBox(EnumSet.range(TreePrior.SPECIES_YULE, TreePrior.SPECIES_BIRTH_DEATH).toArray());
        optionsPanel.addComponentWithLabel("Species Tree Prior:", treePriorCombo);
        
        treePriorCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        fireTreePriorsChanged();
                    }
                }
        );

        optionsPanel.addSeparator();

        optionsPanel.addLabel("Note: the Multispecies Coalescent analysis only needs to select the prior for species tree.");
        
        add(optionsPanel, BorderLayout.NORTH);

    }

    private void fireTreePriorsChanged() {
        if (!settingOptions) {
            frame.setDirty();
        }
    }
 
    private boolean settingOptions = false;

    public void setOptions(BeautiOptions options) {
        this.options = options;

        settingOptions = true;

        treePriorCombo.setSelectedItem(options.nodeHeightPrior);

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
        options.nodeHeightPrior = (TreePrior) treePriorCombo.getSelectedItem();

//        if (options.nodeHeightPrior == ) {
// 
//        } else if (options.nodeHeightPrior == ) {
//
//        }
    }

	@Override
	public JComponent getExportableComponent() {
		// TODO Auto-generated method stub
		return null;
	}

 


}