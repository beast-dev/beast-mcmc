/*
 * GuessTraitDialog.java
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

package dr.app.beauti.traitspanel;

import dr.app.beauti.options.TraitGuesser;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: GuessTraitDialog.java,v 1.4 2009/05/25 13:29:34 rambaut Exp $
 */
public class GuessTraitDialog {

    private JFrame frame;
    private final OptionsPanel optionPanel;

    JLabel descriptionText = new JLabel();

    private final JRadioButton orderRadio = new JRadioButton("Defined by its order", true);
    private final JComboBox orderCombo = new JComboBox(new String[]{"first", "second", "third",
            "fourth", "fourth from last", "third from last", "second from last", "last"});
    private final JTextField delimiterText = new JTextField(6);

    private final JRadioButton regexRadio = new JRadioButton("Defined by regular expression (REGEX)", false);
    private final JTextField regexText = new JTextField(16);

    public GuessTraitDialog(JFrame frame) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addSpanningComponent(descriptionText);

        optionPanel.addLabel("The trait value is given by a part of string in the taxon label that is:");

        optionPanel.addComponents(orderRadio, orderCombo);
        optionPanel.addComponentWithLabel("with delimiter ", delimiterText);
        delimiterText.setEnabled(true);
        optionPanel.addSeparator();

        regexText.setEnabled(false);
        optionPanel.addComponents(regexRadio, regexText);

        ButtonGroup group = new ButtonGroup();
        group.add(orderRadio);
        group.add(regexRadio);
        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
            	delimiterText.setEditable(orderRadio.isSelected());
            	orderCombo.setEnabled(orderRadio.isSelected());
                regexText.setEnabled(regexRadio.isSelected());
            }
        };
        orderRadio.addItemListener(listener);
        regexRadio.addItemListener(listener);
    }

    public void setDescription(String description) {
        descriptionText.setText(description);
    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Guess Trait Value for Taxa");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public void setupGuesserFromDialog(final TraitGuesser guesser) {
        if (orderRadio.isSelected()) {
            int order = orderCombo.getSelectedIndex();
            if (order > 3) {
                order = order - 8; // http://code.google.com/p/beast-mcmc/issues/detail?id=394
            }
            guesser.setGuessType(TraitGuesser.GuessType.DELIMITER);
            guesser.setOrder(order);
            guesser.setDelimiter(delimiterText.getText());
        } else if (regexRadio.isSelected()) {
            guesser.setGuessType(TraitGuesser.GuessType.REGEX);
            guesser.setRegex(regexText.getText());
        } else {
            throw new IllegalArgumentException("unknown radio button selected");
        }
    }

//    private void setTrait(String selectedTrait, TraitGuesser.TraitType selectedTraitType) {
//        guesser.setTraitName(selectedTrait);
//        if (selectedTrait.equalsIgnoreCase(TraitGuesser.Traits.TRAIT_SPECIES.toString())) {
//            guesser.setTraitType(TraitGuesser.TraitType.DISCRETE);
//        } else {
//            guesser.setTraitType(selectedTraitType);
//        }
//    }
}