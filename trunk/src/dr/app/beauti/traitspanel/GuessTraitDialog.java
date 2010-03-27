/*
 * GuessTraitDialog.java
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

import dr.app.beauti.options.TraitData;
import dr.app.beauti.options.TraitGuesser;
import org.virion.jam.panels.OptionsPanel;

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

    private TraitGuesser guesser;
    private JFrame frame;
    private final OptionsPanel optionPanel;

    //    private JComboBox traitTypeComb = new JComboBox(TraitGuesser.TraitType.values());
    
    private final JRadioButton suffixRadio = new JRadioButton("Defined by a suffix, after the", true);
    private final JComboBox suffixOrderCombo = new JComboBox(new String[]{"last", "second from last", "third from last", "fourth from last"});
    private final JTextField suffixText = new JTextField(6);

    private final JRadioButton prefixRadio = new JRadioButton("Defined by a prefix, before the", false);
    private final JComboBox prefixOrderCombo = new JComboBox(new String[]{"first", "second", "third", "fourth"});
    private final JTextField prefixText = new JTextField(6);

    private final JRadioButton regexRadio = new JRadioButton("Defined by regular expression (REGEX)", false);
    private final JTextField regexText = new JTextField(16);

    public GuessTraitDialog(JFrame frame, TraitGuesser guesser) {
        this.frame = frame;
        this.guesser = guesser;

        optionPanel = new OptionsPanel(12, 12);

        JTextField traitNameField = new JTextField(18);
        traitNameField.setText(guesser.getTraitData().getName());
        traitNameField.setEnabled(false);
        optionPanel.addComponentWithLabel("The trait name: ", traitNameField);
//        traitNameField.addKeyListener(new java.awt.event.KeyListener() {
//            public void keyTyped(KeyEvent e) {}
//            public void keyPressed(KeyEvent e) {}
//
//            public void keyReleased(KeyEvent e) {
//               guesser.setTraitName(traitNameField.getText());
//            }
//        } );

//        optionPanel.addComponentWithLabel("The trait type: ", traitTypeComb);
//        traitTypeComb.setEnabled(false);
//        traitTypeComb.addItemListener(
//                new ItemListener() {
//                    public void itemStateChanged(ItemEvent ev) {
//                    	setTrait(traitNameField.getText(), (TraitGuesser.TraitType) traitTypeComb.getSelectedItem());
//                    }
//                }
//        );
        JComboBox traitTypeComb = new JComboBox(TraitData.TraitType.values());
        traitTypeComb.setSelectedItem(guesser.getTraitData().getTraitType());
        traitTypeComb.setEnabled(false);
        optionPanel.addComponentWithLabel("The trait type: ", traitTypeComb);

        optionPanel.addSeparator();
        
        optionPanel.addLabel("The trait value is given by a part of string in the taxon label that is:");
        
        optionPanel.addComponents(suffixRadio, suffixOrderCombo);
        optionPanel.addComponentWithLabel("separator ", suffixText);
        suffixText.setEnabled(true);
        optionPanel.addSeparator();

        optionPanel.addComponents(prefixRadio, prefixOrderCombo);
        optionPanel.addComponentWithLabel("separator", prefixText);
        prefixText.setEnabled(false);
        optionPanel.addSeparator();

        regexText.setEnabled(false);
        optionPanel.addComponents(regexRadio, regexText);
 
        ButtonGroup group = new ButtonGroup();
        group.add(suffixRadio);
        group.add(prefixRadio);
        group.add(regexRadio);
        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
            	suffixText.setEditable(suffixRadio.isSelected());
            	suffixOrderCombo.setEnabled(suffixRadio.isSelected());
            	prefixOrderCombo.setEnabled(prefixRadio.isSelected());
                prefixText.setEnabled(prefixRadio.isSelected());
                regexText.setEnabled(regexRadio.isSelected());
            }
        };
        suffixRadio.addItemListener(listener);
        prefixRadio.addItemListener(listener);
        regexRadio.addItemListener(listener);

//        setupGuesser();
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

    public void setupGuesser() {
        if (suffixRadio.isSelected()) {
            guesser.setGuessType(TraitGuesser.GuessType.SUFFIX);
            guesser.setIndex(suffixOrderCombo.getSelectedIndex());
            guesser.setSeparator(suffixText.getText());
        } else if (prefixRadio.isSelected()) {
            guesser.setGuessType(TraitGuesser.GuessType.PREFIX);
            guesser.setIndex(prefixOrderCombo.getSelectedIndex());
            guesser.setSeparator(prefixText.getText());
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