/*
 * NewTemporalAnalysisDialog.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.app.tracer.analysis;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import jam.framework.DocumentFrame;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

public class NewTemporalAnalysisDialog {

    private JFrame frame;

    private JTextField titleField;

    private WholeNumberField binCountField;

    private RealNumberField minTimeField;
    private RealNumberField maxTimeField;

    private OptionsPanel optionPanel;

    public NewTemporalAnalysisDialog(JFrame frame) {
        this.frame = frame;

        titleField = new JTextField();
        titleField.setColumns(32);

        binCountField = new WholeNumberField(2, 2000);
        binCountField.setValue(100);
        binCountField.setColumns(4);

        maxTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
        maxTimeField.setColumns(12);

        minTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
        minTimeField.setColumns(12);

        optionPanel = new OptionsPanel(12, 12);
    }

    public int showDialog() {

        setArguments();

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Demographic Analysis");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        if (result == JOptionPane.OK_OPTION) {
        }

        return result;
    }

    private void setArguments() {
        optionPanel.removeAll();

        optionPanel.addComponentWithLabel("Title:", titleField);

        optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Number of bins:", binCountField);

        optionPanel.addComponentWithLabel("Minimum time:", minTimeField);
        optionPanel.addComponentWithLabel("Maximum time:", maxTimeField);

//		JLabel label3 = new JLabel(
//				"<html>You can set the age of sampling of the most recent tip in<br>" +
//						"the tree. If this is set to zero then the plot is shown going<br>" +
//						"backwards in time, otherwise forwards in time.</html>");
//		label3.setFont(label3.getFont().deriveFont(((float)label3.getFont().getSize() - 2)));
//		optionPanel.addSpanningComponent(label3);
    }

    public TemporalAnalysisFrame createTemporalAnalysisFrame(DocumentFrame parent) {

        TemporalAnalysisFrame frame = new TemporalAnalysisFrame(parent, titleField.getText(),
                binCountField.getValue(), minTimeField.getValue(), maxTimeField.getValue());
        frame.initialize();
        return frame;
    }

}
