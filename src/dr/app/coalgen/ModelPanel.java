/*
 * ModelPanel.java
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

package dr.app.coalgen;

import dr.app.gui.components.RealNumberField;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class ModelPanel extends JPanel implements Exportable {

    private final CoalGenFrame frame;
    private final CoalGenData data;

    private JComboBox demographicCombo;

    private OptionsPanel optionPanel;

    private RealNumberField[] argumentFields = new RealNumberField[CoalGenData.argumentNames.length];
    private JCheckBox[] argumentCheckBoxes = new JCheckBox[CoalGenData.argumentNames.length];
    private JComboBox[] argumentCombos = new JComboBox[CoalGenData.argumentNames.length];

    public ModelPanel(final CoalGenFrame frame, final CoalGenData data) {

        super();

        this.frame = frame;
        this.data = data;

        setOpaque(false);
        setLayout(new BorderLayout());

        optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
        add(optionPanel, BorderLayout.NORTH);

        demographicCombo = new JComboBox();
        demographicCombo.setOpaque(false);

        for (String demographicModel : CoalGenData.demographicModels) {
            demographicCombo.addItem(demographicModel);
        }

        demographicCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                setDemographicArguments();
                frame.fireModelChanged();
            }
        });

        for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
            argumentFields[i] = new RealNumberField();
            argumentFields[i].setColumns(8);
            argumentFields[i].setValue(data.argumentValues[i]);
            argumentCheckBoxes[i] = new JCheckBox("From Trace:");
            argumentCheckBoxes[i].setOpaque(false);
            argumentCombos[i] = new JComboBox();
            argumentCombos[i].setEnabled(false);
            argumentCombos[i].setOpaque(false);
            argumentCheckBoxes[i].addActionListener(
                    new ArgumentActionListener(argumentCheckBoxes[i], argumentFields[i], argumentCombos[i]));
        }

        setDemographicArguments();
    }

    class ArgumentActionListener implements ActionListener {
        JCheckBox checkBox;
        RealNumberField field;
        JComboBox combo;

        ArgumentActionListener(JCheckBox checkBox, RealNumberField field, JComboBox combo) {
            this.checkBox = checkBox;
            this.field = field;
            this.combo = combo;
        }

        public void actionPerformed(ActionEvent ae) {
            field.setEnabled(!checkBox.isSelected());
            combo.setEnabled(checkBox.isSelected());
        }
    }

    private int findArgument(JComboBox comboBox, String argument) {
        for (int i = 0; i < comboBox.getItemCount(); i++) {
            String item = ((String) comboBox.getItemAt(i)).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    private void setDemographicArguments() {
        optionPanel.removeAll();

        optionPanel.addComponents(new JLabel("Demographic Model:"), demographicCombo);

        optionPanel.addSeparator();

        optionPanel.addLabel("Select the parameter values (or obtain from a trace file):");

        int demo = demographicCombo.getSelectedIndex();

        for (int i = 0; i < data.argumentIndices[demo].length; i++) {
            int k = data.argumentIndices[demo][i];

            JPanel panel = new JPanel(new BorderLayout(6, 6));
            panel.add(argumentFields[k], BorderLayout.WEST);
            panel.add(argumentCheckBoxes[k], BorderLayout.CENTER);
            panel.add(argumentCombos[k], BorderLayout.EAST);
            panel.setOpaque(false);
            optionPanel.addComponentWithLabel(CoalGenData.argumentNames[k] + ":",
                    panel);
        }
        validate();
        repaint();
    }

    public final void tracesChanged() {
        demographicCombo.setSelectedIndex(data.demographicModel);

        if (data.traces != null) {
            for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
                argumentCombos[i].removeAllItems();
                for (int j = 0; j < data.traces.getTraceCount(); j++) {
                    String statistic = data.traces.getTraceName(j);
                    argumentCombos[i].addItem(statistic);
                }

                int index = data.argumentTraces[i];

                for (int j = 0; j < CoalGenData.argumentGuesses[i].length; j++) {
                    if (index != -1) break;

                    index = findArgument(argumentCombos[i], CoalGenData.argumentGuesses[i][j]);
                }
                if (index == -1) {
                    argumentCheckBoxes[i].setEnabled(false);
                    index = 0;
                } else {
                    argumentCheckBoxes[i].setEnabled(true);
                }

                argumentCombos[i].setSelectedIndex(index);
            }
        }
        for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
            argumentFields[i].setValue(data.argumentValues[i]);
        }
        setDemographicArguments();
    }

    public final void collectSettings() {
        data.demographicModel = demographicCombo.getSelectedIndex();
        for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
            data.argumentValues[i] = argumentFields[i].getValue();
            if (argumentCheckBoxes[i].isSelected()) {
                data.argumentTraces[i] = argumentCombos[i].getSelectedIndex();
            } else {
                data.argumentTraces[i] = -1;
            }
        }
    }

    public JComponent getExportableComponent() {
        return this;
    }
}
