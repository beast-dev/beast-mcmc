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

package dr.app.beauti.siteModelsPanel;

import dr.app.beauti.components.continuous.ContinuousComponentOptions;
import dr.app.beauti.components.continuous.ContinuousSubstModelType;
import dr.app.beauti.components.discrete.DiscreteSubstModelType;
import dr.app.beauti.components.dollo.DolloComponentOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.BinaryModelType;
import dr.app.beauti.types.FrequencyPolicyType;
import dr.app.beauti.types.MicroSatModelType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import dr.evomodel.substmodel.AminoAcidModelType;
import dr.evomodel.substmodel.NucModelType;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.*;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class ClonePartitionModelPanel extends OptionsPanel {

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

    private JComboBox sourceModelCombo;
    private JButton cloneModelButton;

    private List<PartitionSubstitutionModel> models;

    public ClonePartitionModelPanel() {

        super(12, (OSType.isMac() ? 6 : 24));
        setOpaque(false);

        cloneModelButton = new JButton();
        cloneModelButton = new JButton("Clone source to selected models");
        cloneModelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                cloneModels();
            }
        });
        // PanelUtils.setupComponent(cloneModelButton);
        cloneModelButton
                .setToolTipText("<html>Clones all the settings of source model<br>"
                        + "to the selected models in the table.</html>");
    }

    private void cloneModels() {
        PartitionSubstitutionModel sourceModel = (PartitionSubstitutionModel)sourceModelCombo.getSelectedItem();
        for (PartitionSubstitutionModel model : models) {
            if (!model.equals(sourceModel)) {
                model.copyFrom(sourceModel);
            }
        }
    }

    public void setOptions(List<PartitionSubstitutionModel> models, List<PartitionSubstitutionModel> sourceModels) {
        this.models = models;

        removeAll();

        sourceModelCombo = new JComboBox();
        for (PartitionSubstitutionModel model : sourceModels) {
            sourceModelCombo.addItem(model);
        }

        PanelUtils.setupComponent(sourceModelCombo);
        sourceModelCombo
                .setToolTipText("<html>Select the substitution model to act as a source<br>to copy to the other selected models.</html>");

        setupPanel();
    }

    /**
     * Lays out the appropriate components in the panel for this partition
     * model.
     */
    private void setupPanel() {

        addSpanningComponent(new JLabel("<html>Select the substitution model to act as a source<br>to copy to the other selected models.</html>"));
        addComponentWithLabel("Source Model:", sourceModelCombo);
        addSeparator();
        addSpanningComponent(cloneModelButton);
    }

}
