/*
 * SelectTraitDialog.java
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

package dr.app.beauti.datapanel;

import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.options.TraitData;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.util.Collection;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class SelectTraitDialog {

    private JFrame frame;

    JComboBox traitCombo;
    JCheckBox copyCheck;
    JTextField nameField;

    OptionsPanel optionPanel;

    public SelectTraitDialog(JFrame frame) {
        this.frame = frame;

        traitCombo = new JComboBox();

        copyCheck = new JCheckBox("Name trait partition:");
        nameField = new JTextField();
        nameField.setColumns(20);

        optionPanel = new OptionsPanel(12, 12);

        copyCheck.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        nameField.setEnabled(copyCheck.isSelected());
                    }
                }
        );

    }

    public int showDialog(Collection<TraitData> traits, String defaultName) {
        optionPanel.removeAll();
        if (traits == null) {
            optionPanel.addSpanningComponent(new JLabel("Create a new data partition using the selected trait(s)."));
            optionPanel.addComponentWithLabel("Name trait partition:", nameField);
            nameField.setText(defaultName != null ? defaultName : "untitled_traits");
            nameField.setEnabled(true);
            nameField.selectAll();
        } else {
            traitCombo.removeAllItems();
            for (Object model : traits) {
                traitCombo.addItem(model);
            }
            optionPanel.addSpanningComponent(new JLabel("Create a new data partition using the following trait."));
            optionPanel.addComponentWithLabel("Trait:", traitCombo);
            optionPanel.addComponents(copyCheck, nameField);
            nameField.setEnabled(copyCheck.isSelected());
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Partition");
        dialog.pack();

        int result;
        boolean isValid;

        do {
            dialog.setVisible(true);

            isValid = true;
            result = JOptionPane.CANCEL_OPTION;

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }
            if (result != JOptionPane.CANCEL_OPTION) {
                String name = getName().trim();
                if (name.isEmpty()) {
                    isValid = false;
                }
            }
        } while (!isValid);

        return result;
    }

    public TraitData getTrait() {
        return (TraitData) traitCombo.getSelectedItem();
    }

    public boolean getMakeCopy() {
        return copyCheck.isSelected();
    }

    public String getName() {
        return nameField.getText();
    }

}