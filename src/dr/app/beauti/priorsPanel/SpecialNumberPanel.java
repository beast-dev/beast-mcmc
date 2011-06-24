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

import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.text.NumberFormat;

/**
 * @author Walter Xie
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SpecialNumberPanel.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class SpecialNumberPanel extends JPanel implements ActionListener {
    private final JButton b1, b2, b3, b4, b5;
    private final JLabel label;

    private RealNumberField selectedField = null;

    SpecialNumberPanel() {
        super(new BorderLayout());

        label = new JLabel("Set a special value in the text fields above:");
        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        panel.setOpaque(false);

        b1 = new JButton(NumberFormat.getNumberInstance().format(Double.POSITIVE_INFINITY));
        PanelUtils.setupComponent(b1);
        b1.setFocusable(false);
        b1.setActionCommand(RealNumberField.POSITIVE_INFINITY);
        b1.addActionListener(this);

        b2 = new JButton(NumberFormat.getNumberInstance().format(Double.NEGATIVE_INFINITY));
        PanelUtils.setupComponent(b2);
        b2.setFocusable(false);
        b2.setActionCommand(RealNumberField.NEGATIVE_INFINITY);
        b2.addActionListener(this);

        b3 = new JButton("MAX");
        PanelUtils.setupComponent(b3);
        b3.setFocusable(false);
        b3.setActionCommand(RealNumberField.MAX_VALUE);
        b3.addActionListener(this);

        b4 = new JButton("MIN");
        PanelUtils.setupComponent(b4);
        b4.setFocusable(false);
        b4.setActionCommand(RealNumberField.MIN_VALUE);
        b4.addActionListener(this);

        b5 = new JButton(RealNumberField.NaN);
        PanelUtils.setupComponent(b5);
        b5.setFocusable(false);
        b5.setActionCommand(RealNumberField.NaN);
        b5.addActionListener(this);

        b1.setToolTipText("Click to set 'Positive Infinity' in the selected text field.");
        b2.setToolTipText("Click to set 'Negative Infinity' in the selected text field.");
        b3.setToolTipText("Click to set the 'Maximum numerical value' in the selected text field.");
        b3.setToolTipText("Click to set the 'Minimum numerical value' in the selected text field.");
        b3.setToolTipText("Click to set 'Not a Number' in the selected text field.");

        //Add Components to this container, using the default FlowLayout.
        panel.add(b1);
        panel.add(b2);
        panel.add(b3);
        panel.add(b4);
        panel.add(b5);

        add(label, BorderLayout.NORTH);
        add(panel, BorderLayout.CENTER);
    }

    public void setSelectedField(RealNumberField selectedField) {
        this.selectedField = selectedField;
    }

    @Override
    public void setEnabled(boolean enabled) {
        super.setEnabled(enabled);
        b1.setEnabled(enabled);
        b2.setEnabled(enabled);
        b3.setEnabled(enabled);
        b4.setEnabled(enabled);
        b5.setEnabled(enabled);
        label.setEnabled(enabled);
    }

    public void actionPerformed(ActionEvent e) {
        if (selectedField != null) {
            selectedField.setText(e.getActionCommand());
        }
    }
}
