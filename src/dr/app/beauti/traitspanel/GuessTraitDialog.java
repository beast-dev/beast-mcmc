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

import dr.app.beauti.PanelUtils;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.PartitionModel;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.options.TreePrior;
import dr.gui.chart.Axis;
import dr.gui.chart.JChart;
import dr.gui.chart.LinearAxis;
import dr.gui.chart.PDFPlot;
import dr.math.*;
import dr.util.NumberFormatter;
import dr.evolution.datatype.*;
import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.ExpansionModel;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.LogisticGrowthModel;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.YuleModel;

import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.Calendar;
import java.util.GregorianCalendar;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: GuessTraitDialog.java,v 1.4 2009/05/25 13:29:34 rambaut Exp $
 */
public class GuessTraitDialog {

    private JFrame frame;
    private TraitGuesser guesser;
    
    private final OptionsPanel optionPanel;

    private final JComboBox traitSelectionCombo;
    
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
        
        traitSelectionCombo = new JComboBox(TraitGuesser.Traits.values());
        optionPanel.addComponentWithLabel("The selected trait is: ", traitSelectionCombo);
        traitSelectionCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	setTrait((TraitGuesser.Traits) traitSelectionCombo.getSelectedItem());
                    }
                }
        );
        optionPanel.addSeparator();
        
        optionPanel.addLabel("The trait value is given by a part of string in the taxon label that is:");
        
        optionPanel.addComponents(suffixRadio, suffixOrderCombo);
        optionPanel.addComponentWithLabel("seperator ", suffixText);
        suffixText.setEnabled(true);
        optionPanel.addSeparator();

        optionPanel.addComponents(prefixRadio, prefixOrderCombo);
        optionPanel.addComponentWithLabel("seperator", prefixText);
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
    	
    	setTrait((TraitGuesser.Traits) traitSelectionCombo.getSelectedItem());   	
    	
        if (suffixRadio.isSelected()) {
            guesser.guessType = TraitGuesser.GuessType.SUFFIX;
            guesser.index = suffixOrderCombo.getSelectedIndex();
            guesser.separator = suffixText.getText();
        } else if (prefixRadio.isSelected()) {
            guesser.guessType = TraitGuesser.GuessType.PREFIX;
            guesser.index = prefixOrderCombo.getSelectedIndex();
            guesser.separator = prefixText.getText();
        } else if (regexRadio.isSelected()) {
            guesser.guessType = TraitGuesser.GuessType.REGEX;
            guesser.regex = regexText.getText();
        } else {
            throw new IllegalArgumentException("unknown radio button selected");
        }
    }
    
    private void setTrait(TraitGuesser.Traits selectedTrait) {
    	
		switch (selectedTrait) {
			case TRAIT_SPECIES:
				guesser.traitAnalysisType = TraitGuesser.Traits.TRAIT_SPECIES;
	    		guesser.traitType = TraitGuesser.TraitType.DISCRETE;
				break;
			default:
				throw new IllegalArgumentException("unknown trait selected");
		}
    	
    }
}