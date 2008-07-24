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

package dr.app.beauti;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.PartitionModel;
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

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class SelectModelDialog {

    private JFrame frame;

    JComboBox modelCombo;
    JCheckBox copyCheck;
    JTextField nameField;

    OptionsPanel optionPanel;

    public SelectModelDialog(JFrame frame) {
        this.frame = frame;

        modelCombo = new JComboBox();

        copyCheck = new JCheckBox("Make copy with name:");
        nameField = new JTextField();
        nameField.setColumns(32);
        nameField.setEnabled(false);

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Partition Model:", modelCombo);
        optionPanel.addComponents(copyCheck, nameField);

        copyCheck.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        nameField.setEnabled(copyCheck.isSelected());
                    }
                }
        );

    }

    public int showDialog(Object[] models) {

        modelCombo.removeAllItems();
        for (Object model : models) {
            modelCombo.addItem(model);
        }
        
        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Model");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public PartitionModel getModel() {
        return (PartitionModel)modelCombo.getSelectedItem();
    }

    public boolean getMakeCopy() {
        return copyCheck.isSelected();
    }

    public String getName() {
        return nameField.getText();
    }

}