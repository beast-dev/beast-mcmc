/*
 * DiscretePriorDialog.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.priorsPanel;

import dr.app.beauti.enumTypes.PriorType;
import dr.app.beauti.options.Parameter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class DiscretePriorDialog {

    private JFrame frame;

    public static String[] priors = {
            "Uniform",
            "Poisson"};

    private String[] argumentNames = new String[]{
            "Poisson Mean", "Zero Offset"
    };

    private JComboBox priorCombo;
    private int[][] argumentIndices = {{}, {0, 1}};
    private RealNumberField initialField = new RealNumberField();
    private RealNumberField[] argumentFields = new RealNumberField[argumentNames.length];
    private OptionsPanel optionPanel;

    private Parameter parameter;

    public DiscretePriorDialog(JFrame frame) {
        this.frame = frame;

        priorCombo = new JComboBox(priors);

        initialField.setColumns(8);
        for (int i = 0; i < argumentNames.length; i++) {
            argumentFields[i] = new RealNumberField();
            argumentFields[i].setColumns(8);
        }

        optionPanel = new OptionsPanel(12, 12);
    }

    public int showDialog(final Parameter parameter) {

        this.parameter = parameter;

        priorCombo.setSelectedIndex(parameter.priorType == PriorType.POISSON_PRIOR ? 1 : 0);

        if (!parameter.isStatistic) {
            initialField.setRange(parameter.lower, parameter.upper);
            initialField.setValue(parameter.initial);
        }

        setArguments();
        setupComponents();

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Prior for Parameter");
        dialog.pack();

        priorCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                setupComponents();
                dialog.pack();
                dialog.repaint();
            }
        });

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        if (result == JOptionPane.OK_OPTION) {
            getArguments();
        }

        return result;
    }

    private void setArguments() {
        argumentFields[0].setRange(0.0, Double.MAX_VALUE);
        argumentFields[0].setValue(parameter.poissonMean);
        argumentFields[1].setValue(parameter.poissonOffset);

    }

    private void getArguments() {
        parameter.priorType = priorCombo.getSelectedIndex() == 0 ? PriorType.UNIFORM_PRIOR : PriorType.POISSON_PRIOR;

        if (initialField.getValue() != null) parameter.initial = initialField.getValue();

        switch (parameter.priorType) {
            case UNIFORM_PRIOR:
                if (argumentFields[0].getValue() != null) parameter.uniformLower = argumentFields[0].getValue();
                if (argumentFields[1].getValue() != null) parameter.uniformUpper = argumentFields[1].getValue();
                break;
            case POISSON_PRIOR:
                if (argumentFields[0].getValue() != null) parameter.poissonMean = argumentFields[0].getValue();
                if (argumentFields[1].getValue() != null) parameter.poissonOffset = argumentFields[1].getValue();
                break;
            default:
                throw new IllegalArgumentException("Unknown prior index");
        }
    }

    private void setupComponents() {
        optionPanel.removeAll();

        optionPanel.addSpanningComponent(new JLabel("Select prior distribution for " + parameter.getName()));

        optionPanel.addComponents(new JLabel("Prior Distribution:"), priorCombo);
        int priorType = priorCombo.getSelectedIndex();

        optionPanel.addSeparator();

        for (int i = 0; i < argumentIndices[priorType].length; i++) {
            int k = argumentIndices[priorType][i];
            optionPanel.addComponentWithLabel(argumentNames[k] + ":", argumentFields[k]);
        }


        if (!parameter.isStatistic) {
            optionPanel.addSeparator();
            optionPanel.addComponents(new JLabel("Initial Value:"), initialField);
        }
    }

}
