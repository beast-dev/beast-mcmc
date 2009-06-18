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

package dr.app.beauti.tipdatepanel;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.PartitionModel;
import dr.app.beauti.options.DateGuesser;
import dr.gui.chart.Axis;
import dr.gui.chart.JChart;
import dr.gui.chart.LinearAxis;
import dr.gui.chart.PDFPlot;
import dr.math.*;
import dr.util.NumberFormatter;
import dr.evolution.datatype.*;
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
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class GuessDatesDialog {

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private final JRadioButton orderRadio = new JRadioButton("Defined by its order", true);
    private final JComboBox orderCombo = new JComboBox(new String[]{"first", "second", "third",
            "fourth", "fourth from last",
            "third from last", "second from last", "last"});

    private final JRadioButton prefixRadio = new JRadioButton("Defined by a prefix", false);
    private final JTextField prefixText = new JTextField(16);

    private final JRadioButton regexRadio = new JRadioButton("Defined by regular expression (REGEX)", false);
    private final JTextField regexText = new JTextField(16);

    private final JCheckBox offsetCheck = new JCheckBox("Add the following value to each: ", false);
    private final RealNumberField offsetText = new RealNumberField();

    private final JCheckBox unlessCheck = new JCheckBox("...unless less than:", false);
    private final RealNumberField unlessText = new RealNumberField();

    private final RealNumberField offset2Text = new RealNumberField();

    public GuessDatesDialog(JFrame frame) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);

        optionPanel.addLabel("The date is given by a numerical field in the taxon label that is:");


        optionPanel.addComponents(orderRadio, orderCombo);
        optionPanel.addSeparator();

        prefixText.setEnabled(false);
        optionPanel.addComponents(prefixRadio, prefixText);
        optionPanel.addSeparator();

        regexText.setEnabled(false);
        optionPanel.addComponents(regexRadio, regexText);
        optionPanel.addSeparator();

        offsetText.setValue(1900);
        offsetText.setColumns(16);
        offsetText.setEnabled(false);
        offsetCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                offsetText.setEnabled(offsetCheck.isSelected());
            }
        });
        optionPanel.addComponents(offsetCheck, offsetText);

        Calendar calendar = GregorianCalendar.getInstance();

        int year = calendar.get(Calendar.YEAR) - 1999;
        unlessText.setValue(year);
        unlessText.setColumns(16);
        unlessText.setEnabled(false);
        optionPanel.addComponents(unlessCheck, unlessText);

        offset2Text.setValue(2000);
        offset2Text.setColumns(16);
        offset2Text.setEnabled(false);
        optionPanel.addComponentWithLabel("...in which case add:", offset2Text);

        unlessCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                unlessText.setEnabled(unlessCheck.isSelected());
                offset2Text.setEnabled(unlessCheck.isSelected());
            }
        });

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

        final JDialog dialog = optionPane.createDialog(frame, "Guess Dates for Taxa");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public void setupGuesser(DateGuesser guesser) {
        if (orderRadio.isSelected()) {
            guesser.guessType = DateGuesser.GuessType.ORDER;
            guesser.order = orderCombo.getSelectedIndex();
            guesser.fromLast = false;
            if (guesser.order > 3) {
                guesser.fromLast = true;
                guesser.order = 8 - guesser.order - 1;
            }

        } else if (prefixRadio.isSelected()) {
            guesser.guessType = DateGuesser.GuessType.PREFIX;
            guesser.prefix = prefixText.getText();
        } else if (regexRadio.isSelected()) {
            guesser.guessType = DateGuesser.GuessType.REGEX;
            guesser.regex = regexText.getText();
        } else {
            throw new IllegalArgumentException("unknown radio button selected");
        }

        guesser.offset = 0.0;
        guesser.unlessLessThan = 0.0;
        if (offsetCheck.isSelected()) {
            guesser.offset = offsetText.getValue();
            if (unlessCheck.isSelected()) {
                guesser.unlessLessThan = unlessText.getValue();
                guesser.offset2 = offset2Text.getValue();

            }
        }

    }
}