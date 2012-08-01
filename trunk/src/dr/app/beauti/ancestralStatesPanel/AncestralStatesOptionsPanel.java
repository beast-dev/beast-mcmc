/*
 * PartitionModelPanel.java
 *
 * Copyright (c) 2002-2011 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.ancestralStatesPanel;

import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.types.SequenceErrorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxa;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class AncestralStatesOptionsPanel extends OptionsPanel {

    private static final String ROBUST_COUNTING_TOOL_TIP = "<html>"
            + "Enable counting of reconstructed number of substitutions as described in<br>"
            + "Minin & Suchard (2008). These will be annotated directly in the<br>"
            + "logged trees.</html>";

    private static final String DNDS_ROBUST_COUNTING_TOOL_TIP = "<html>"
            + "Enable counting of synonymous and non-synonymous substitution as described in<br>"
            + "O'Brien, Minin & Suchard (2009) and Lemey, Minin, Bielejec, Kosakovsky-Pond &<br>"
            + "Suchard (in preparation). This model requires a 3-partition codon model to be<br>"
            + "selected in the Site model for this partition.</html>";

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

    private final AbstractPartitionData partition;

    private JCheckBox ancestralReconstructionCheck = new JCheckBox(
            "Reconstruct states at all ancestors");
    private JCheckBox mrcaReconstructionCheck = new JCheckBox(
            "Reconstruct states at ancestor:");
    private JComboBox mrcaReconstructionCombo = new JComboBox();
    private JCheckBox countingCheck = new JCheckBox(
            "Reconstruct state change counts");
    private JCheckBox dNdSRobustCountingCheck = new JCheckBox(
            "Reconstruct synonymous/non-synonymous change counts");


    private JTextArea dNnSText = new JTextArea(
            "This model requires a 3-partition codon model to be selected in " +
                    "the Site model for this partition before it can be selected.");

    // dNdS robust counting is automatic if RC is turned on for a codon
    // partitioned data set.
//    private JCheckBox dNdSCountingCheck = new JCheckBox(
//            "Reconstruct synonymous/non-synonymous counts");

    final BeautiOptions options;

    //    JComboBox errorModelCombo = new JComboBox(EnumSet.range(SequenceErrorType.NO_ERROR, SequenceErrorType.BASE_ALL).toArray());
    JComboBox errorModelCombo = new JComboBox(SequenceErrorType.values());

    AncestralStatesComponentOptions ancestralStatesComponent;
    SequenceErrorModelComponentOptions sequenceErrorComponent;

    public AncestralStatesOptionsPanel(final AncestralStatesPanel ancestralStatesPanel, final BeautiOptions options, final AbstractPartitionData partition) {

        super(12, (OSType.isMac() ? 6 : 24));
        setOpaque(false);

        this.partition = partition;
        this.options = options;

        PanelUtils.setupComponent(ancestralReconstructionCheck);
        ancestralReconstructionCheck
                .setToolTipText("<html>"
                        + "Reconstruct posterior realizations of the states at ancestral nodes.<br>" +
                        "These will be annotated directly in the logged trees.</html>");

        PanelUtils.setupComponent(mrcaReconstructionCheck);
        mrcaReconstructionCheck
                .setToolTipText("<html>"
                        + "Reconstruct posterior realizations of the states at a specific common<br>" +
                        "ancestor defined by a taxon set. This will be recorded in the log file.</html>");
        PanelUtils.setupComponent(mrcaReconstructionCombo);
        mrcaReconstructionCombo
                .setToolTipText("<html>"
                        + "Reconstruct posterior realizations of the states at a specific common.<br>" +
                        "ancestor defined by a taxon set. This will be recorded in the log file.</html>");

        PanelUtils.setupComponent(countingCheck);
        countingCheck.setToolTipText(ROBUST_COUNTING_TOOL_TIP);

        PanelUtils.setupComponent(dNdSRobustCountingCheck);
        dNdSRobustCountingCheck.setToolTipText(DNDS_ROBUST_COUNTING_TOOL_TIP);

        // ////////////////////////
        PanelUtils.setupComponent(errorModelCombo);
        errorModelCombo.setToolTipText("<html>Select how to model sequence error or<br>"
                + "post-mortem DNA damage.</html>");


        // Set the initial options
        ancestralStatesComponent = (AncestralStatesComponentOptions)options.getComponentOptions(AncestralStatesComponentOptions.class);
//        ancestralStatesComponent.createParameters(options);
        ancestralReconstructionCheck.setSelected(ancestralStatesComponent.reconstructAtNodes(partition));
        mrcaReconstructionCheck.setSelected(ancestralStatesComponent.reconstructAtMRCA(partition));
        mrcaReconstructionCombo.setSelectedItem(ancestralStatesComponent.getMRCATaxonSet(partition));
        countingCheck.setSelected(ancestralStatesComponent.isCountingStates(partition));
        dNdSRobustCountingCheck.setSelected(ancestralStatesComponent.dNdSRobustCounting(partition));

        sequenceErrorComponent = (SequenceErrorModelComponentOptions)options.getComponentOptions(SequenceErrorModelComponentOptions.class);
//        sequenceErrorComponent.createParameters(options); // this cannot create correct param here, because of improper design
        errorModelCombo.setSelectedItem(sequenceErrorComponent.getSequenceErrorType(partition));

        setupPanel();

        ItemListener listener = new ItemListener() {
            public void itemStateChanged(final ItemEvent itemEvent) {
                optionsChanged();
                ancestralStatesPanel.fireModelChanged();

                // The following is only necessary is simpleCounting XOR robustCounting
//                if (itemEvent.getItem() == countingCheck) {
//                    boolean enableRC = !countingCheck.isSelected() && ancestralStatesComponent.dNdSRobustCountingAvailable(partition);
//                    dNdSRobustCountingCheck.setEnabled(enableRC);
//                    dNnSText.setEnabled(enableRC);
//                }
//
//                if (itemEvent.getItem() == dNdSRobustCountingCheck) {
//                    boolean enableSimpleCounting = !dNdSRobustCountingCheck.isSelected();
//                    countingCheck.setEnabled(enableSimpleCounting);
//                }
            }
        };

        ancestralReconstructionCheck.addItemListener(listener);
        mrcaReconstructionCheck.addItemListener(listener);
        mrcaReconstructionCombo.addItemListener(listener);
        countingCheck.addItemListener(listener);
        dNdSRobustCountingCheck.addItemListener(listener);

        errorModelCombo.addItemListener(listener);
    }

    private void optionsChanged() {
        if (isUpdating) return;

        ancestralStatesComponent.setReconstructAtNodes(partition, ancestralReconstructionCheck.isSelected());
        ancestralStatesComponent.setReconstructAtMRCA(partition, mrcaReconstructionCheck.isSelected());
        mrcaReconstructionCombo.setEnabled(mrcaReconstructionCheck.isSelected());
        if (mrcaReconstructionCombo.getSelectedIndex() == 0) {
            // root node
            ancestralStatesComponent.setMRCATaxonSet(partition, null);
        } else {
            String text = (String) mrcaReconstructionCombo.getSelectedItem();
            String taxonSetId = text.substring(5,text.length() - 1);

            ancestralStatesComponent.setMRCATaxonSet(partition, taxonSetId);
        }
        ancestralStatesComponent.setCountingStates(partition, countingCheck.isSelected());
//        ancestralStatesComponent.setDNdSRobustCounting(partition, robustCountingCheck.isSelected());
        ancestralStatesComponent.setDNdSRobustCounting(partition, dNdSRobustCountingCheck.isSelected());

        sequenceErrorComponent.setSequenceErrorType(partition, (SequenceErrorType)errorModelCombo.getSelectedItem());
        sequenceErrorComponent.createParameters(options);
    }

    /**
     * Lays out the appropriate components in the panel for this partition
     * model.
     */
    void setupPanel() {

        isUpdating = true;

        String selectedItem = (String)mrcaReconstructionCombo.getSelectedItem();

        if (mrcaReconstructionCombo.getItemCount() > 0) {
            mrcaReconstructionCombo.removeAllItems();
        }
        mrcaReconstructionCombo.addItem("Tree Root");
        if (options.taxonSets.size() > 0) {
            for (Taxa taxonSet : options.taxonSets) {
                mrcaReconstructionCombo.addItem("MRCA("+ taxonSet.getId() + ")");
            }
            if (selectedItem != null) {
                mrcaReconstructionCombo.setSelectedItem(selectedItem);
            }
        }
        mrcaReconstructionCombo.setEnabled(mrcaReconstructionCheck.isSelected());

        boolean ancestralReconstructionAvailable = true;
        boolean countingAvailable = true;
        boolean dNdSRobustCountingAvailable = false;
        boolean errorModelAvailable = false;

        switch (partition.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                errorModelAvailable = true;
                dNdSRobustCountingAvailable = true; // but will be disabled if not codon partitioned
                break;
            case DataType.AMINO_ACIDS:
            case DataType.GENERAL:
            case DataType.TWO_STATES:
                break;
            case DataType.CONTINUOUS:
                countingAvailable = false;
                break;
            case DataType.MICRO_SAT:
                ancestralReconstructionAvailable = false;
                countingAvailable = false;
                break;
            default:
                throw new IllegalArgumentException("Unsupported data type");

        }

        removeAll();

        if (ancestralReconstructionAvailable) {
            addSpanningComponent(new JLabel("Ancestral State Reconstruction:"));

            addComponent(ancestralReconstructionCheck);

            FlowLayout layout = new FlowLayout(FlowLayout.LEFT);
            layout.setHgap(0);
            JPanel panel = new JPanel(layout);
            panel.setOpaque(false);
            panel.setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
            panel.add(mrcaReconstructionCheck);
            panel.add(mrcaReconstructionCombo);
            addComponent(panel);
        }

        if (countingAvailable) {
            if (ancestralReconstructionAvailable) {
                addSeparator();
            }
            addSpanningComponent(new JLabel("State Change Count Reconstruction:"));

            JTextArea text = new JTextArea(
                    "Select this option to reconstruct counts of state changes using " +
                            "Markov Jumps. This approached is described in Minin & Suchard (2008).");
            text.setColumns(40);
            PanelUtils.setupComponent(text);
            addComponent(text);

            addComponent(countingCheck);

            boolean enableSimpleCounting = true;

            // TODO Simple counting is currently not available for codon partitioned models due to BEAUti limitation
            if (ancestralStatesComponent.dNdSRobustCountingAvailable(partition)) {
                enableSimpleCounting = false;
                countingCheck.setSelected(false);
            }
            countingCheck.setEnabled(enableSimpleCounting);

            if (dNdSRobustCountingAvailable) {
                addSeparator();
                text = new JTextArea(
                        "Select this option to reconstruct counts of synonymous and nonsynonymous " +
                                "changes using Robust Counting. This approached is described in O'Brien, Minin " +
                                "& Suchard (2009) and Lemey, Minin, Bielejec, Kosakovsky-Pond & Suchard " +
                                "(in preparation):");
                text.setColumns(40);
                PanelUtils.setupComponent(text);
                addComponent(text);

                addComponent(dNdSRobustCountingCheck);

                dNnSText.setColumns(40);
                dNnSText.setBorder(BorderFactory.createEmptyBorder(0, 32, 0, 0));
                PanelUtils.setupComponent(dNnSText);
                addComponent(dNnSText);

                boolean enableRC = ancestralStatesComponent.dNdSRobustCountingAvailable(partition);
                // && !ancestralStatesComponent.isCountingStates(partition);
                dNdSRobustCountingCheck.setEnabled(enableRC);
                dNnSText.setEnabled(enableRC);
                if (!enableRC) {
                    dNdSRobustCountingCheck.setSelected(false);
                }
            }
        }

        if (errorModelAvailable) {
            if (ancestralReconstructionAvailable || countingAvailable) {
                addSeparator();
            }
            addSpanningComponent(new JLabel("Sequence error model:"));
            addComponentWithLabel("Error Model:", errorModelCombo);
        }
        isUpdating = false;

    }

    private boolean isUpdating = false;

}
