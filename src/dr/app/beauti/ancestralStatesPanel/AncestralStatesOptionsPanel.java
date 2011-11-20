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
import dr.app.beauti.options.*;
import dr.app.beauti.types.*;
import dr.app.beauti.util.PanelUtils;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxa;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class AncestralStatesOptionsPanel extends OptionsPanel {

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

    private final AbstractPartitionData partition;

    private JCheckBox ancestralReconstructionCheck = new JCheckBox(
            "Reconstruct states at all ancestors");
    private JCheckBox mrcaReconstructionCheck = new JCheckBox(
            "Reconstruct states at MRCA for Taxon Set:");
    private JComboBox mrcaReconstructionCombo = new JComboBox();
    private JCheckBox robustCountingCheck = new JCheckBox(
            "Reconstruct state change counts");
    private JCheckBox dNdSCountingCheck = new JCheckBox(
            "Reconstruct synonymous/non-synonymous counts");

    final BeautiOptions options;

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

        PanelUtils.setupComponent(robustCountingCheck);
        robustCountingCheck
                .setToolTipText("<html>"
                        + "Enable counting of reconstructed number of substitutions as described in<br>" +
                        "Minin & Suchard (in preparation). These will be annotated directly in the<br>"
                        + "logged trees.</html>");

        PanelUtils.setupComponent(dNdSCountingCheck);
        dNdSCountingCheck
                .setToolTipText("<html>"
                        + "Enable counting of synonymous and non-synonymous substitution as described in<br>" +
                        "Lemey, Minin, Bielejec, Kosakovsky-Pond & Suchard (in preparation). This model<br>"
                        + "requires a 3-partition codon model to be selected, above.</html>");

        robustCountingCheck.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                dNdSCountingCheck.setEnabled(robustCountingCheck.isSelected());
//                if (robustCountingCheck.isSelected()) {
//                    if (checkRobustCounting()) {
//                        setRobustCountingModel();
//                    }
//
//                } else {
//                    removeRobustCountingModel();
//                }

            }// END: actionPerformed
        });

        // ////////////////////////
        PanelUtils.setupComponent(errorModelCombo);
        errorModelCombo.setToolTipText("<html>Select how to model sequence error or<br>"
                + "post-mortem DNA damage.</html>");


        // Set the initial options
        ancestralStatesComponent = (AncestralStatesComponentOptions)options.getComponentOptions(AncestralStatesComponentOptions.class);
        ancestralReconstructionCheck.setSelected(ancestralStatesComponent.reconstructAtNodes(partition));
        mrcaReconstructionCheck.setSelected(ancestralStatesComponent.reconstructAtMRCA(partition));
        mrcaReconstructionCombo.setSelectedItem(ancestralStatesComponent.getMRCATaxonSet(partition));
        robustCountingCheck.setSelected(ancestralStatesComponent.robustCounting(partition));
        dNdSCountingCheck.setSelected(ancestralStatesComponent.dNdSRobustCounting(partition));

        sequenceErrorComponent = (SequenceErrorModelComponentOptions)options.getComponentOptions(SequenceErrorModelComponentOptions.class);
        errorModelCombo.setSelectedItem(sequenceErrorComponent.getSequenceErrorType(partition));

        setupPanel();

        ItemListener listener = new ItemListener() {
            public void itemStateChanged(final ItemEvent itemEvent) {
                optionsChanged();
                ancestralStatesPanel.fireModelChanged();
            }
        };

        ancestralReconstructionCheck.addItemListener(listener);
        mrcaReconstructionCheck.addItemListener(listener);
        mrcaReconstructionCombo.addItemListener(listener);
        robustCountingCheck.addItemListener(listener);
        dNdSCountingCheck.addItemListener(listener);

        errorModelCombo.addItemListener(listener);
    }

    private void optionsChanged() {
        ancestralStatesComponent.setReconstructAtNodes(partition, ancestralReconstructionCheck.isSelected());
        ancestralStatesComponent.setReconstructAtMRCA(partition, mrcaReconstructionCheck.isSelected());
        ancestralStatesComponent.setMRCATaxonSet(partition, (Taxa) mrcaReconstructionCombo.getSelectedItem());
        ancestralStatesComponent.setRobustCounting(partition, robustCountingCheck.isSelected());
        ancestralStatesComponent.setDNdSRobustCounting(partition, dNdSCountingCheck.isSelected());

        sequenceErrorComponent.setSequenceErrorType(partition, (SequenceErrorType)errorModelCombo.getSelectedItem());
    }

    /**
     * Lays out the appropriate components in the panel for this partition
     * model.
     */
    void setupPanel() {

        String selectedItem = (String)mrcaReconstructionCombo.getSelectedItem();

        if (mrcaReconstructionCombo.getItemCount() > 0) {
            mrcaReconstructionCombo.removeAllItems();
        }
        if (options.taxonSets.size() > 0) {
            for (Taxa taxonSet : options.taxonSets) {
                mrcaReconstructionCombo.addItem(taxonSet.getId());
            }
            if (selectedItem != null) {
                mrcaReconstructionCombo.setSelectedItem(selectedItem);
            }
            mrcaReconstructionCheck.setEnabled(true);
            mrcaReconstructionCombo.setEnabled(true);
        } else {
            mrcaReconstructionCheck.setEnabled(false);
            mrcaReconstructionCombo.setEnabled(false);
        }

        boolean ancestralReconstruction = true;
        boolean robustCounting = true;
        boolean errorModel = false;

        switch (partition.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                errorModel = true;
                break;
            case DataType.CONTINUOUS:
                robustCounting = false;
                break;
            case DataType.MICRO_SAT:
                ancestralReconstruction = false;
                robustCounting = false;
                break;
            default:
                throw new IllegalArgumentException("Unsupported data type");

        }

        if (ancestralReconstruction || robustCounting) {
            addSpanningComponent(new JLabel("Ancestral State Reconstruction:"));
        }
        if (ancestralReconstruction) {
            addComponent(ancestralReconstructionCheck);
            JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
            panel.setOpaque(false);
            panel.setBorder(BorderFactory.createEmptyBorder());
            panel.add(mrcaReconstructionCheck);
            panel.add(mrcaReconstructionCombo);
            addComponent(panel);
        }

        if (robustCounting) {
            addComponent(robustCountingCheck);
            addComponent(dNdSCountingCheck);
        }

        if (errorModel) {
            if (ancestralReconstruction || robustCounting) {
                addSeparator();
            }
            addComponentWithLabel("Sequence Error Model:", errorModelCombo);
        }
    }


//    private boolean checkRobustCounting() {
//
//        if (heteroCombo.getSelectedIndex() == 0  && codingCombo.getSelectedIndex() == 2) {
//            return true;
//        } else {
//            SwingUtilities.invokeLater(new Runnable() {
//
//                public void run() {
//
//                    String msg = String.format("Wrong settings. \n"
//                            + "Set site heterogeneity model to none \n"
//                            + "and partition into 3 codon position.");
//
//                    JOptionPane.showMessageDialog(PanelUtils
//                            .getActiveFrame(), msg, "Error",
//                            JOptionPane.ERROR_MESSAGE);
//
//                    robustCountingCheck.setSelected(false);
//
//                }
//            });
//
//            return false;
//        }
//    }// END: checkRobustCounting
//
//    private void setRobustCountingModel() {
//        DnDsComponentOptions comp = (DnDsComponentOptions) model
//                .getOptions().getComponentOptions(
//                        DnDsComponentOptions.class);
//
//        // Add model to ComponentOptions
//        comp.addPartition(model);
//    }
//
//    private void removeRobustCountingModel() {
//        DnDsComponentOptions comp = (DnDsComponentOptions) model
//                .getOptions().getComponentOptions(
//                        DnDsComponentOptions.class);
//
//        // Remove model from ComponentOptions
//        comp.removePartition(model);
//    }
//
}
