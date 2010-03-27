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
package dr.app.beauti.siteModelsPanel;

import dr.app.beauti.options.PartitionDiscreteTraitSubstModel;
import dr.app.beauti.util.PanelUtils;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 *@author Walter Xie
 */
public class DiscreteTraitSiteModelPanel extends PartitionModelPanel {

    private JComboBox discreteTraitSiteModelCombo = new JComboBox(PartitionDiscreteTraitSubstModel.LocationSubstModelType.values());
    private JCheckBox activateBSSVS = new JCheckBox("Activate BSSVS");


    final PartitionDiscreteTraitSubstModel discreteTraitModel;

    private boolean settingOptions = false;


    public DiscreteTraitSiteModelPanel(final PartitionDiscreteTraitSubstModel discreteTraitModel) {
    	super(discreteTraitModel);

        this.discreteTraitModel = discreteTraitModel;

        PanelUtils.setupComponent(discreteTraitSiteModelCombo);
        discreteTraitSiteModelCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                discreteTraitModel.setLocationSubstType(
                        (PartitionDiscreteTraitSubstModel.LocationSubstModelType) discreteTraitSiteModelCombo.getSelectedItem());
            }
        });

        PanelUtils.setupComponent(activateBSSVS);
        activateBSSVS.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                discreteTraitModel.setActivateBSSVS(activateBSSVS.isSelected());
            }
        });

        setupPanel();
    }

    private void setupPanel() {

        removeAll();

        addComponentWithLabel("Discrete Trait Substitution Model:", discreteTraitSiteModelCombo);

        addComponent(activateBSSVS);

        validate();
        repaint();
    }

    public void setOptions() {
        settingOptions = true;

        discreteTraitSiteModelCombo.setSelectedItem(discreteTraitModel.getLocationSubstType());
        activateBSSVS.setSelected(discreteTraitModel.isActivateBSSVS());

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions() {
        if (settingOptions) return;

    }

}