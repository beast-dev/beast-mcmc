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
    
    private final JRadioButton orderRadio = new JRadioButton("Defined by its order", true);
    private final JComboBox orderCombo = new JComboBox(new String[]{"first", "second", "third",
            "fourth", "fourth from last", "third from last", "second from last", "last"});

    private final JRadioButton prefixRadio = new JRadioButton("Defined by a prefix", false);
    private final JTextField prefixText = new JTextField(16);

    private final JRadioButton regexRadio = new JRadioButton("Defined by regular expression (REGEX)", false);
    private final JTextField regexText = new JTextField(16);

    public GuessTraitDialog(JFrame frame, TraitGuesser guesser) {
        this.frame = frame;
        this.guesser = guesser;

        optionPanel = new OptionsPanel(12, 12);

        traitSelectionCombo = new JComboBox(TraitGuesser.TraitAnalysisType.values());
        optionPanel.addComponentWithLabel("The selected trait is: ", traitSelectionCombo);
        traitSelectionCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	setTrait((TraitGuesser.TraitAnalysisType) traitSelectionCombo.getSelectedItem());
                    }
                }
        );
        optionPanel.addSeparator();
        
        optionPanel.addLabel("The trait is given by a part of string in the taxon label that is:");

        optionPanel.addComponents(orderRadio, orderCombo);
        optionPanel.addSeparator();

        prefixText.setEnabled(false);
        optionPanel.addComponents(prefixRadio, prefixText);
        optionPanel.addSeparator();

        regexText.setEnabled(false);
        optionPanel.addComponents(regexRadio, regexText);
 
        ButtonGroup group = new ButtonGroup();
        group.add(orderRadio);
        group.add(prefixRadio);
        group.add(regexRadio);
        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                orderCombo.setEnabled(orderRadio.isSelected());
                prefixText.setEnabled(prefixRadio.isSelected());
                regexText.setEnabled(regexRadio.isSelected());
            }
        };
        orderRadio.addItemListener(listener);
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

        final JDialog dialog = optionPane.createDialog(frame, "Guess Trait for Taxa");
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
    	
    	setTrait((TraitGuesser.TraitAnalysisType) traitSelectionCombo.getSelectedItem());   	
    	
        if (orderRadio.isSelected()) {
            guesser.guessType = TraitGuesser.GuessType.ORDER;
            guesser.order = orderCombo.getSelectedIndex();
            guesser.fromLast = false;
            if (guesser.order > 3) {
                guesser.fromLast = true;
                guesser.order = 8 - guesser.order - 1;
            }

        } else if (prefixRadio.isSelected()) {
            guesser.guessType = TraitGuesser.GuessType.PREFIX;
            guesser.prefix = prefixText.getText();
        } else if (regexRadio.isSelected()) {
            guesser.guessType = TraitGuesser.GuessType.REGEX;
            guesser.regex = regexText.getText();
        } else {
            throw new IllegalArgumentException("unknown radio button selected");
        }
    }
    
    private void setTrait(TraitGuesser.TraitAnalysisType selectedTrait) {
    	
		switch (selectedTrait) {
			case SPECIES_ANALYSIS:
				guesser.traitAnalysisType = TraitGuesser.TraitAnalysisType.SPECIES_ANALYSIS;
	    		guesser.traitType = TraitGuesser.TraitType.DISCRETE;
				break;
			default:
				throw new IllegalArgumentException("unknown trait selected");
		}
    	
    }
}