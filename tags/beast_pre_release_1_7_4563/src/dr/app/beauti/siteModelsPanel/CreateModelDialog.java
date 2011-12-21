/*
 * CreateModelDialog.java
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

package dr.app.beauti.siteModelsPanel;

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.TwoStates;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 * @deprecated
 */
public class CreateModelDialog {

    private JFrame frame;

    public static DataType[] dataTypes = {
            Nucleotides.INSTANCE,
            AminoAcids.INSTANCE,
            TwoStates.INSTANCE
    };


    JTextField nameField;
    JComboBox dataTypeCombo;

    OptionsPanel optionPanel;

    public CreateModelDialog(JFrame frame) {
        this.frame = frame;

        nameField = new JTextField();
        nameField.setColumns(20);

        dataTypeCombo = new JComboBox(dataTypes);

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Name:", nameField);
        optionPanel.addComponentWithLabel("Data Type:", dataTypeCombo);


    }

    public int showDialog() {

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

    public String getName() {
        return nameField.getText();
    }

    public DataType getDataType() {
        return (DataType) dataTypeCombo.getSelectedItem();
    }
}