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

import dr.app.beauti.components.continuous.ContinuousSubstModelType;
import dr.app.beauti.components.discrete.DiscreteSubstModelType;
import dr.app.beauti.components.dnds.DnDsComponentOptions;
import dr.app.beauti.components.dollo.DolloComponentOptions;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
import dr.app.beauti.options.*;
import dr.app.beauti.siteModelsPanel.SiteModelsPanel;
import dr.app.beauti.types.*;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import dr.evomodel.substmodel.AminoAcidModelType;
import dr.evomodel.substmodel.NucModelType;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.EnumSet;
import java.util.logging.Logger;

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
            "Reconstruct states at MRCA:");
    private JComboBox mrcaReconstructionCombo = new JComboBox();
    private JCheckBox robustCountingCheck = new JCheckBox(
            "Reconstruct state change counts");
    private JCheckBox dNdSCountingCheck = new JCheckBox(
            "Reconstruct synonymous/non-synonymous counts");

    JComboBox errorModelCombo = new JComboBox(SequenceErrorType.values());

    public AncestralStatesOptionsPanel(final BeautiOptions options, final AbstractPartitionData partition) {

        super(12, (OSType.isMac() ? 6 : 24));

        this.partition = partition;

        PanelUtils.setupComponent(ancestralReconstructionCheck);
        ancestralReconstructionCheck
                .setToolTipText("<html>"
                        + "Reconstruct posterior realizations of the states at ancestral nodes.<br>" +
                        "These will be annotated directly in the logged trees.</html>");

        PanelUtils.setupComponent(mrcaReconstructionCheck);
        mrcaReconstructionCheck
                .setToolTipText("<html>"
                        + "Reconstruct posterior realizations of the states at a specific common.<br>" +
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

        errorModelCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
//                fireModelsChanged();
            }
        });
        errorModelCombo.setToolTipText("<html>Select how to model sequence error or<br>"
                + "post-mortem DNA damage.</html>");

        SequenceErrorModelComponentOptions comp = (SequenceErrorModelComponentOptions)options.getComponentOptions(SequenceErrorModelComponentOptions.class);

        errorModelCombo.setSelectedItem(comp.getSequenceErrorType(partition));

        setupPanel();
        setOpaque(false);
    }

    /**
     * Sets the components up according to the partition model - but does not
     * layout the top level options panel.
     */
    public void setOptions() {
//        ancestralReconstructionCheck.setSelected();
    }

    /**
     * Lays out the appropriate components in the panel for this partition
     * model.
     */
    private void setupPanel() {

        boolean ancestralReconstruction = true;
        boolean robustCounting = true;
        boolean errorModel = false;

        switch (partition.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                errorModel = true;
                break;
            case DataType.CONTINUOUS:
                robustCounting = false;
                errorModel = false;
                break;
            case DataType.MICRO_SAT:
                ancestralReconstruction = false;
                robustCounting = false;
                errorModel = false;
                break;
            default:
                throw new IllegalArgumentException("Unsupported data type");

        }

        if (ancestralReconstruction) {
            addComponent(ancestralReconstructionCheck);
            JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
            panel.add(mrcaReconstructionCheck);
            panel.add(mrcaReconstructionCombo);
            addComponent(panel);
        }

        if (robustCounting) {
            addComponent(robustCountingCheck);
            addComponent(dNdSCountingCheck);
        }

        if (errorModel) {
            addSeparator();
            addComponentWithLabel("Sequence Error Model:", errorModelCombo);
        }

        setOptions();
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
