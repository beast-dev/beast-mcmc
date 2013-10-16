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

import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.types.PopulationSizeModelType;
import dr.app.beauti.util.PanelUtils;

import dr.app.util.OSType;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.*;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: SpeciesTreesPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 * @deprecated
 */
public class SpeciesTreesPanel extends OptionsPanel {

	private static final long serialVersionUID = -2768091530149898538L;

	private JComboBox treePriorCombo = new JComboBox(EnumSet.range(TreePriorType.SPECIES_YULE, TreePriorType.SPECIES_BIRTH_DEATH).toArray());
    private JComboBox populationSizeCombo = new JComboBox(PopulationSizeModelType.values());

    private final PartitionTreePrior partitionTreePrior;
//    private boolean settingOptions = false;

    public SpeciesTreesPanel(final PartitionTreePrior partitionTreePrior) {
    	super(12, (OSType.isMac() ? 6 : 24));

    	this.partitionTreePrior = partitionTreePrior;

    	PanelUtils.setupComponent(treePriorCombo);
    	addComponentWithLabel("Species Tree Prior:", treePriorCombo);

        treePriorCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	partitionTreePrior.setNodeHeightPrior( (TreePriorType) treePriorCombo.getSelectedItem());
                    }
                }
        );

        PanelUtils.setupComponent(populationSizeCombo);
    	addComponentWithLabel("Population Size Model:", populationSizeCombo);

        populationSizeCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	partitionTreePrior.setPopulationSizeModel((PopulationSizeModelType) populationSizeCombo.getSelectedItem());
                    }
                }
        );
//        addSeparator();

        addLabel("Note: *BEAST only needs to select the prior for species tree.");

        validate();
        repaint();

    }

    public void setOptions() {
//        settingOptions = true;

        treePriorCombo.setSelectedItem(partitionTreePrior.getNodeHeightPrior());
        populationSizeCombo.setSelectedItem(partitionTreePrior.getPopulationSizeModel());

//        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions() {

    }

	// @Override
	public JComponent getExportableComponent() {
		// TODO Auto-generated method stub
		return null;
	}
}