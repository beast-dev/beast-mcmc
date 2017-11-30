/*
 * PartitionClockModelPanel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.clockmodelspanel;

import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.types.ClockDistributionType;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.util.PanelUtils;
import dr.app.util.OSType;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 */
public class PartitionClockModelPanel extends OptionsPanel {

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

    private JComboBox clockTypeCombo = new JComboBox(EnumSet.range(
            ClockType.STRICT_CLOCK, ClockType.FIXED_LOCAL_CLOCK).toArray());
    private JComboBox clockDistributionCombo = new JComboBox (new ClockDistributionType[] {
            ClockDistributionType.LOGNORMAL,
            ClockDistributionType.GAMMA,
//            ClockDistributionType.CAUCHY,
            ClockDistributionType.EXPONENTIAL,
            //ClockDistributionType.MODEL_AVERAGING
    });
    private JCheckBox continuousQuantileCheck = new JCheckBox("Use continuous quantile parameterization.");

    private JLabel modelAveragingInfo = new JLabel(
            "<html>Using the Bayesian Model Averaging (BMA) approach on the<br>" +
                    "available relaxed clock models as described by Li & Drummond (2012)<br>" +
                    " Mol. Biol. Evol. 29:751-761.</html>");

    protected final PartitionClockModel model;

    public PartitionClockModelPanel(final PartitionClockModel partitionModel) {

        super(12, (OSType.isMac() ? 6 : 24));

        this.model = partitionModel;

        PanelUtils.setupComponent(clockTypeCombo);
        clockTypeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setClockType((ClockType) clockTypeCombo.getSelectedItem());
                setupPanel();
            }
        });
        clockTypeCombo.setToolTipText("<html>Select the type of molecular clock model.</html>");

        clockTypeCombo.setSelectedItem(model.getClockType());

        PanelUtils.setupComponent(clockDistributionCombo);
        clockDistributionCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setClockDistributionType((ClockDistributionType) clockDistributionCombo.getSelectedItem());
                if (clockDistributionCombo.getSelectedItem() == ClockDistributionType.MODEL_AVERAGING) {
                    model.setPerformModelAveraging(true);
                    continuousQuantileCheck.setSelected(true);
                    continuousQuantileCheck.setEnabled(false);
                    addComponent(modelAveragingInfo);
                } else {
                    model.setPerformModelAveraging(false);
                    continuousQuantileCheck.setEnabled(true);
                    remove(modelAveragingInfo);
                }
                model.setContinuousQuantile(continuousQuantileCheck.isSelected());
            }
        });
        clockDistributionCombo.setToolTipText("<html>Select the distribution that describes the variation in rate.</html>");

        clockDistributionCombo.setSelectedItem(model.getClockDistributionType());

        PanelUtils.setupComponent(continuousQuantileCheck);
        continuousQuantileCheck.setToolTipText("<html>" +
                "Select this option to use the continuous quantile form of the relaxed<br>" +
                "clock model described by Li & Drummond (2012) MBE 29:751-61 instead of<br>" +
                "the discretized categorical form.<html>");
        continuousQuantileCheck.setSelected(model.isContinuousQuantile());
        continuousQuantileCheck.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        model.setContinuousQuantile(continuousQuantileCheck.isSelected());
                    }
                });

        setupPanel();
        setOpaque(false);

    }


    /**
     * Lays out the appropriate components in the panel for this partition model.
     */
    public void setupPanel() {
        removeAll();
        addComponentWithLabel("Clock Type:", clockTypeCombo);

        switch (model.getClockType()) {
            case STRICT_CLOCK:
                break;

            case UNCORRELATED:
                addComponent(new JLabel(
                        "<html>" +
                                "Using the uncorrelated relaxed clock model of Drummond, Ho, Phillips & <br>" +
                                "Rambaut (2006) PLoS Biology 4, e88.<html>"));
                addComponentWithLabel("Relaxed Distribution:", clockDistributionCombo);
                addComponent(continuousQuantileCheck);
                break;

            case AUTOCORRELATED:
                addComponentWithLabel("Relaxed Distribution:", clockDistributionCombo);
                break;

            case RANDOM_LOCAL_CLOCK:
            case FIXED_LOCAL_CLOCK:
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model type");

        }

    }

    /**
     * Sets the components up according to the partition model - but does not
     * layout the top level options panel.
     */
    public void setOptions() {

        if (model == null) {
            return;
        }
        clockTypeCombo.setSelectedItem(model.getClockType());
        clockDistributionCombo.setSelectedItem(model.getClockDistributionType());

        setupPanel();
        setOpaque(false);

    }

}
