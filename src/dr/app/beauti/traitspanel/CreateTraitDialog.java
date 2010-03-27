/*
 * PriorDialog.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.beauti.traitspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.TraitData;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

/**
 * @author Andrew Rambaut
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class CreateTraitDialog {

    private final BeautiFrame frame;

//    JTextField nameField;
    JComboBox nameCombo;
    JComboBox typeCombo;

    OptionsPanel optionPanel;

    public CreateTraitDialog(BeautiFrame frame) {
        this.frame = frame;

//        nameField = new JTextField(TraitGuesser.Traits.TRAIT_SPECIES.toString());
//        nameField.setColumns(20);

        nameCombo = new JComboBox(TraitData.Traits.values());
        typeCombo = new JComboBox(TraitData.TraitType.values());

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Name:", nameCombo);
        optionPanel.addComponentWithLabel("Type:", typeCombo);


    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        boolean done;
        int result = JOptionPane.CANCEL_OPTION;

        do {
            final JDialog dialog = optionPane.createDialog(frame, "Create New Trait");
            dialog.pack();

            dialog.setVisible(true);

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }
            done = true;
            if (result != JOptionPane.CANCEL_OPTION) {
                done = frame.validateTraitName(getName());
            }
        } while (!done);

        return result;
    }

    public String getName() {
//        return nameField.getText();
        return nameCombo.getSelectedItem().toString();
    }

    public TraitData.TraitType getType() {
        return (TraitData.TraitType) typeCombo.getSelectedItem();
    }
}