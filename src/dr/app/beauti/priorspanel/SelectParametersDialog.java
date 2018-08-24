/*
 * SelectParametersDialog.java
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

package dr.app.beauti.priorspanel;

import dr.app.beauti.options.Parameter;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.util.List;


/**
 * @author Andrew Rambaut
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class SelectParametersDialog {

    private final JFrame frame;

    private JComboBox parameterCombo = new JComboBox();

    OptionsPanel optionPanel = new OptionsPanel(12, 12);


    public SelectParametersDialog(final JFrame frame) {
        this.frame = frame;
    }

    public int showDialog(String message, List<Parameter> parameterList) {

        optionPanel.removeAll();

        if (message != null && !message.isEmpty()) {
            optionPanel.addSpanningComponent(new JLabel(message));
        }

        parameterCombo.removeAllItems();
        for (Parameter parameter : parameterList) {
            parameterCombo.addItem(parameter);
        }
        optionPanel.addComponentWithLabel("Parameter:", parameterCombo);

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        int result = JOptionPane.CANCEL_OPTION;

        final JDialog dialog = optionPane.createDialog(frame, "Add Parameter");
        dialog.pack();

        dialog.setVisible(true);

        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public Parameter getSelectedParameter() {
        return (Parameter)parameterCombo.getSelectedItem();
    }

}