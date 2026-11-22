/*
 * PartitionClockModelPanel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.clockmodelspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.types.ClockDistributionType;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.util.PanelUtils;
import dr.app.util.OSType;
import dr.evolution.util.Taxa;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

/**
 * @author Andrew Rambaut
 */
public class PartitionClockModelPanel extends OptionsPanel {

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

    private JComboBox<ClockType> clockTypeCombo = new JComboBox<>();

    private JComboBox clockDistributionCombo = new JComboBox (new ClockDistributionType[] {
            ClockDistributionType.LOGNORMAL,
            ClockDistributionType.GAMMA,
//            ClockDistributionType.CAUCHY,
            ClockDistributionType.EXPONENTIAL,
            //ClockDistributionType.MODEL_AVERAGING
    });
    private JComboBox clockHMCDistributionCombo = new JComboBox (new ClockDistributionType[] {
            ClockDistributionType.LOGNORMAL
    });
    private JCheckBox continuousQuantileCheck = new JCheckBox("Use continuous quantile parameterization.");

    private final JTextArea localClocksText = new JTextArea();

    private JLabel modelAveragingInfo = new JLabel(
            "<html>Using the Bayesian Model Averaging (BMA) approach on the<br>" +
                    "available relaxed clock models as described by Li & Drummond (2012)<br>" +
                    " Mol. Biol. Evol. 29:751-761.</html>");

    protected final PartitionClockModel model;

    public PartitionClockModelPanel(final PartitionClockModel partitionModel) {

        super(12, (OSType.isMac() ? 6 : 24));

        this.model = partitionModel;

        for (ClockType clockType : EnumSet.range(ClockType.STRICT_CLOCK, ClockType.MIXED_EFFECTS_CLOCK)) {
            clockTypeCombo.addItem(clockType);
            /*if (clockType == ClockType.STRICT_CLOCK || clockType == ClockType.HMC_CLOCK) {
                clockTypeCombo.addItem(new JSeparator(JSeparator.HORIZONTAL));
            }*/
        }

        //PanelUtils.setupComponent(clockTypeCombo);
        clockTypeCombo.putClientProperty("JButton.buttonType", "textured");
        clockTypeCombo.setRenderer(new ClockComboBoxRenderer());

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
                "the discretized categorical form.</html>");
        continuousQuantileCheck.setSelected(model.isContinuousQuantile());
        continuousQuantileCheck.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        model.setContinuousQuantile(continuousQuantileCheck.isSelected());
                    }
                });

        localClocksText.setEditable(false);
        localClocksText.setLineWrap(true);
        localClocksText.setWrapStyleWord(true);
        localClocksText.setOpaque(false);
        localClocksText.setFont(getFont());
        localClocksText.setColumns(40);

        setupPanel();
        setOpaque(false);

    }

    /**
     * Describes which clades the fixed local clock will give their own rate: every taxon set
     * marked as monophyletic in the Taxa panel (the same rule ClockModelGenerator uses).
     */
    private void updateLocalClocksText() {
        BeautiOptions options = model.getOptions();
        List<String> clades = new ArrayList<>();
        boolean severalTrees = options.getPartitionTreeModels().size() > 1;
        for (Taxa taxonSet : options.taxonSets) {
            if (Boolean.TRUE.equals(options.taxonSetsMono.get(taxonSet))) {
                StringBuilder clade = new StringBuilder("\u2022 ").append(taxonSet.getId());
                if (Boolean.TRUE.equals(options.taxonSetsIncludeStem.get(taxonSet))) {
                    clade.append(" (including stem)");
                }
                if (severalTrees && options.taxonSetsTreeModel.get(taxonSet) != null) {
                    clade.append(" on tree ").append(options.taxonSetsTreeModel.get(taxonSet).getName());
                }
                clades.add(clade.toString());
            }
        }

        String text;
        if (clades.isEmpty()) {
            text = "No monophyletic taxon sets defined. \n" +
                    "In the '" + BeautiFrame.TAXON_SETS + "' panel, define a taxon set for each clade that should " +
                    "have its own rate and select its 'Mono?' box.";
        } else {
            text = "A separate rate will be estimated for " +
                    (clades.size() == 1 ? "taxon set:" : "taxon sets:") +
                    "\n" + String.join("\n", clades);
        }
        localClocksText.setText(text);
        localClocksText.setRows(Math.max(2, clades.size() + 1));
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
                                "Rambaut (2006) PLoS Biology 4, e88.</html>"));
                addComponentWithLabel("Relaxed Distribution:", clockDistributionCombo);
                addComponent(continuousQuantileCheck);
                break;

            case HMC_CLOCK:
                addComponent(new JLabel(
                        "<html>" +
                                "Using the Hamiltonian Monte Carlo relaxed clock model of Ji, Zhang, Holbrook,<br>" +
                                "Nishimura, Baele, Rambaut, Lemey & Suchard (2020) Mol Biol Evol 37, 3047–3060.</html>"));
                addComponentWithLabel("Relaxed Distribution:", clockHMCDistributionCombo);
                break;

            case AUTOCORRELATED:
                addComponentWithLabel("Relaxed Distribution:", clockDistributionCombo);
                break;

            case FIXED_LOCAL_CLOCK:
                addComponent(new JLabel(
                        "<html>" +
                                "Using the a priori local clock model of Yoder & Yang (2000)<br>" +
                                "Mol Biol Evol 17, 1081\u20131090.</html>"));
                updateLocalClocksText();
                addComponentWithLabel("Local clocks:", localClocksText);
                break;

            case SHRINKAGE_LOCAL_CLOCK:
            case RANDOM_LOCAL_CLOCK:
            case MIXED_EFFECTS_CLOCK:
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

    /**
     *
     */
    static class ClockComboBoxRenderer extends JLabel implements ListCellRenderer<Object> {

        public Component getListCellRendererComponent(JList<? extends Object> list,
                                                      Object value, int index, boolean isSelected, boolean cellHasFocus) {

            if (value == ClockType.UNCORRELATED || value == ClockType.HMC_CLOCK){
                setBackground(new Color(0, 100, 255, 30));
                setOpaque(true);
            } else {
                setOpaque(false);
            }

            setText("  " + value.toString() + " ");

            return this;
        }
    }

}
