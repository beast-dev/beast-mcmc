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

package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.BeautiFrame;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class DefaultPriorDialog {

    private BeautiFrame frame;
    private PriorsPanel priorsPanel;

    public DefaultPriorDialog(BeautiFrame frame) {
        this.frame = frame;
        priorsPanel = new PriorsPanel(frame, true);
    }

    public JButton findButton(String label, Container container) {
        for (Component comp : container.getComponents()) {
            if (comp instanceof JButton) {
                JButton button = (JButton)comp;
                if (button.getText().equals(label)) {
                    return button;
                }
            } else if (comp instanceof Container) {
                JButton button = findButton(label, (Container) comp);
                if (button != null) {
                    return button;
                }
            }

        }
        return null;
    }

    public boolean showDialog(BeautiOptions options) {
        priorsPanel.setParametersList(options);

        Object[] buttons;
        JOptionPane optionPane;

        String title;

        if (priorsPanel.hasUndefinedPrior) {
            title = "Undefined Priors";
        } else {
            title = "Unchanged Default Priors";
        }
        buttons = new String[] {"Continue", "Cancel"};
        optionPane = new JOptionPane(priorsPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                buttons,
                buttons[0]);

        JButton button = findButton("Continue", optionPane);
        if (button != null) {
            priorsPanel.setContinueButton(button);
        }

//       if (priorsPanel.hasUndefinedPrior) {
//            buttons = new String[] {"OK"};
//            title = "Undefined Priors";
//            optionPane = new JOptionPane(priorsPanel,
//                    JOptionPane.ERROR_MESSAGE,
//                    JOptionPane.OK_OPTION,
//                    null,
//                    buttons,
//                    buttons[0]);
//        } else {
//            buttons = new String[] {"Continue", "Cancel"};
//            title = "Unchanged Default Priors";
//            optionPane = new JOptionPane(priorsPanel,
//                    JOptionPane.PLAIN_MESSAGE,
//                    JOptionPane.OK_CANCEL_OPTION,
//                    null,
//                    buttons,
//                    buttons[0]);
//        }

        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        optionPane.setPreferredSize(new java.awt.Dimension(800, 600));

        final JDialog dialog = optionPane.createDialog(frame, title);
        dialog.pack();
        dialog.setResizable(true);
        dialog.setVisible(true);

        return optionPane.getValue() != null && optionPane.getValue().equals("Continue");
    }
}