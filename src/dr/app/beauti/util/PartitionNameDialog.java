/*
 * PartitionNameDialog.java
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

package dr.app.beauti.util;

import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionNameDialog {

    private JFrame frame;
    private final OptionsPanel optionPanel;

    JLabel descriptionText = new JLabel();

    private final JTextField valueText = new JTextField(16);

    public PartitionNameDialog(JFrame frame) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);

        optionPanel.addSpanningComponent(descriptionText);
        optionPanel.addComponentWithLabel("New Name: ", valueText);
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

        final JDialog dialog = optionPane.createDialog(frame, "Rename partition");
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
        return valueText.getText().trim();
    }
}