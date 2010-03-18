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

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.options.TraitsOptions;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class CreateTraitDialog {

    private JFrame frame;


//    JTextField nameField;
    JComboBox nameCombo;
    JComboBox typeCombo;

    OptionsPanel optionPanel;

    public CreateTraitDialog(JFrame frame) {
        this.frame = frame;

//        nameField = new JTextField(TraitGuesser.Traits.TRAIT_SPECIES.toString());
//        nameField.setColumns(20);

        nameCombo = new JComboBox(TraitGuesser.Traits.values());
        typeCombo = new JComboBox(TraitGuesser.TraitType.values());

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Name:", nameCombo);
        optionPanel.addComponentWithLabel("Type:", typeCombo);


    }

    public int showDialog(BeautiOptions options) {

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
                done = validate(options);
            }
        } while (!done);

        return result;
    }

    private boolean validate(BeautiOptions options) {
        // check that the name is valid
        if (getName().trim().length() == 0) {
            Toolkit.getDefaultToolkit().beep();
            return false;
        }

        // disallow a trait called 'date'
        if (getName().equalsIgnoreCase("date")) {
            JOptionPane.showMessageDialog(frame,
                    "This trait name has a special meaning. Use the 'Tip Date' panel\n" +
                            " to set dates for taxa.",
                    "Reserved trait name",
                    JOptionPane.WARNING_MESSAGE);

            return false;
        }

        // check that the trait name doesn't exist
        if (TraitsOptions.containTrait(getName())) {
            int option = JOptionPane.showConfirmDialog(frame,
                    "A trait of this name already exists. Do you wish to replace\n" +
                            "it with this new trait? This may result in the loss or change\n" +
                            "in trait values for the taxa.",
                    "Overwrite trait?",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE);

            if (option == JOptionPane.NO_OPTION) {
                return false;
            }

        }

        return true;
    }

    public String getName() {
//        return nameField.getText();
        return nameCombo.getSelectedItem().toString();
    }

    public TraitGuesser.TraitType getType() {
        return (TraitGuesser.TraitType) typeCombo.getSelectedItem();
    }
}