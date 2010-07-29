/*
 * SpecialNumberPanel.java
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


package dr.app.beauti.priorsPanel;

import javax.swing.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.text.NumberFormat;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: SpecialNumberPanel.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class SpecialNumberPanel extends JPanel implements ActionListener {
    protected JButton b1, b2, b3, b4, b5;

    private PriorDialog priorDialog;

    SpecialNumberPanel (PriorDialog pd) {
        this.priorDialog = pd;

        b1 = new JButton(NumberFormat.getNumberInstance().format(Double.POSITIVE_INFINITY));
        b1.setActionCommand(RealNumberField.POSITIVE_INFINITY);
        b1.addActionListener(this);

        b2 = new JButton(NumberFormat.getNumberInstance().format(Double.NEGATIVE_INFINITY));
        b2.setActionCommand(RealNumberField.NEGATIVE_INFINITY);
        b2.addActionListener(this);

        b3 = new JButton("MAX");
        b3.setActionCommand(RealNumberField.MAX_VALUE);
        b3.addActionListener(this);

        b4 = new JButton("MIN");
        b4.setActionCommand(RealNumberField.MIN_VALUE);
        b4.addActionListener(this);

        b5 = new JButton(RealNumberField.NaN);
        b5.setActionCommand(RealNumberField.NaN);
        b5.addActionListener(this);

        b1.setToolTipText("Click to set Double.POSITIVE_INFINITY in the selected text field.");
        b2.setToolTipText("Click to set Double.POSITIVE_INFINITY in the selected text field.");
        b3.setToolTipText("Click to set Double.MAX_VALUE in the selected text field.");
        b3.setToolTipText("Click to set Double.MIN_VALUE in the selected text field.");
        b3.setToolTipText("Click to set Double.NaN in the selected text field.");

        //Add Components to this container, using the default FlowLayout.
        add(b1);
        add(b2);
        add(b3);
        add(b4);
        add(b5);
    }

    public void actionPerformed(ActionEvent e) {
        if (priorDialog.getSelectedField() == null) {
            JOptionPane.showMessageDialog(this, "Please select a text field above !",
                    "No text field selected",
                    JOptionPane.INFORMATION_MESSAGE);
        } else {
            priorDialog.getSelectedField().setText(e.getActionCommand());
        }

        priorDialog.setSelectedField(null); // to force user to select the correct field
    }
}
