/*
 * BadTraitFormatDialog.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.datapanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.traitspanel.CreateTraitDialog;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class BadTraitFormatDialog {

    private final BeautiFrame frame;
    private final JButton exampleButton = new JButton("Show example of mapping file format");


    OptionsPanel optionPanel;

    public BadTraitFormatDialog(final BeautiFrame frame) {
        this.frame = frame;
        optionPanel = new OptionsPanel(12, 12);

        optionPanel.addSpanningComponent(new JLabel("An error occurred when importing the traits file."));
        optionPanel.addSpanningComponent(new JLabel("It is likely the file is not properly formatted."));

        optionPanel.addComponent(exampleButton);

        exampleButton.setEnabled(true);
        exampleButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                CreateTraitDialog.showExampleTraitFormat(frame);
            }
        });
    }

    public int showDialog() {
        String[] options = {"OK"};

        int choice = JOptionPane.showOptionDialog(frame,
                optionPanel,
                "Could not import traits.",
                JOptionPane.DEFAULT_OPTION,
                JOptionPane.ERROR_MESSAGE,
                null,
                options,
                options[0]);
        return choice;
    }

}
