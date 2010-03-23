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

import dr.app.beauti.options.PhylogeographicOptions;
import dr.app.beauti.util.PanelUtils;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 *@author Walter Xie
 */
public class LocationSiteModelPanel extends OptionsPanel {

    private JComboBox locationSiteModelCombo = new JComboBox(PhylogeographicOptions.LocationSubstModelType.values());
    private JCheckBox activeBSSVS = new JCheckBox("Active BSSVS");


    final PhylogeographicOptions phylogeographicOptions;

    private boolean settingOptions = false;


    public LocationSiteModelPanel(final PhylogeographicOptions phylogeographicOptions) {
    	super(12, 30);

        this.phylogeographicOptions = phylogeographicOptions;

        PanelUtils.setupComponent(locationSiteModelCombo);
        locationSiteModelCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                phylogeographicOptions.setLocationSubstType(
                        (PhylogeographicOptions.LocationSubstModelType) locationSiteModelCombo.getSelectedItem());
            }
        });

        PanelUtils.setupComponent(activeBSSVS);
        activeBSSVS.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                phylogeographicOptions.setActiveBSSVS(activeBSSVS.isSelected());
            }
        });

        setupPanel();
    }

    private void setupPanel() {

        removeAll();

        addComponentWithLabel("Location Substitution Model:", locationSiteModelCombo);

        addComponent(activeBSSVS);

        validate();
        repaint();
    }

    public void setOptions() {
        settingOptions = true;

        locationSiteModelCombo.setSelectedItem(phylogeographicOptions.getLocationSubstType());
        activeBSSVS.setSelected(phylogeographicOptions.isActiveBSSVS());

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions() {
        if (settingOptions) return;

    }

}