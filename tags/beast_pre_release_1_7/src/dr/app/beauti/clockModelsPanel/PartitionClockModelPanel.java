/*
 * PartitionModelPanel.java
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

package dr.app.beauti.clockModelsPanel;

import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.types.ClockDistributionType;
import dr.app.beauti.util.PanelUtils;
import dr.app.util.OSType;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @deprecated
 */
public class PartitionClockModelPanel extends OptionsPanel {

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

//    private JComboBox clockTypeCombo = new JComboBox(ClockType.values());
    private JComboBox clockDistributionCombo = new JComboBox(ClockDistributionType.values());

    protected final PartitionClockModel model;

    public PartitionClockModelPanel(final PartitionClockModel partitionModel) {

        super(12, (OSType.isMac() ? 6 : 24));

        this.model = partitionModel;

//        PanelUtils.setupComponent(clockTypeCombo);
//        clockTypeCombo.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent ev) {
//                model.setClockType((ClockType) clockTypeCombo.getSelectedItem());
//                setupPanel();
//            }
//        });
//        clockTypeCombo.setToolTipText("<html>Select the type of molecular clock model.</html>");
//
//        clockTypeCombo.setSelectedItem(model.getClockType());

        PanelUtils.setupComponent(clockDistributionCombo);
        clockDistributionCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setClockDistributionType((ClockDistributionType) clockDistributionCombo.getSelectedItem());
                setupPanel();
            }
        });
        clockDistributionCombo.setToolTipText("<html>Select the distribution that describes the variation in rate.</html>");

        clockDistributionCombo.setSelectedItem(model.getClockDistributionType());

        setupPanel();
        setOpaque(false);
    }


    /**
     * Lays out the appropriate components in the panel for this partition model.
     */
    public void setupPanel() {
        removeAll();
//        addComponentWithLabel("Clock Type:", clockTypeCombo);

        switch (model.getClockType()) {
            case STRICT_CLOCK:
                break;

            case UNCORRELATED:
            case AUTOCORRELATED:
                addComponentWithLabel("Relaxed Distribution:", clockDistributionCombo);
                break;

            case RANDOM_LOCAL_CLOCK:
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");

        }
    }

}
