/*
 * PriorDialog.java
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
import dr.app.util.OSType;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorDialog implements AbstractPriorDialog {

    private final JFrame frame;

    private final PriorSettingsPanel priorSettingsPanel;
    private Parameter parameter;

    public PriorDialog(JFrame frame) {
        this.frame = frame;

         priorSettingsPanel = new PriorSettingsPanel(frame);

    }

    /**
     * Set the parameter to be controlled
     *
     * @param parameter
     */
    public void setParameter(final Parameter parameter) {
        this.parameter = parameter;
        priorSettingsPanel.setParameter(parameter);
    }

    public int showDialog() {

        JPanel panel = new JPanel(new BorderLayout(0, 6));
        panel.add(new JLabel("Select prior distribution for " + parameter.getName()), BorderLayout.NORTH);
        panel.add(priorSettingsPanel, BorderLayout.CENTER);

        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        JOptionPane optionPane = new JOptionPane(scrollPane,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Prior for Parameter " + parameter.getName());

        priorSettingsPanel.setDialog(dialog);

        if (OSType.isMac()) {
            dialog.setMinimumSize(new Dimension(dialog.getBounds().width, 300));
        } else {
            Toolkit tk = Toolkit.getDefaultToolkit();
            Dimension d = tk.getScreenSize();
            if (d.height < 700 && priorSettingsPanel.getHeight() > 450) {
                dialog.setSize(new java.awt.Dimension(priorSettingsPanel.getWidth() + 100, 550));
            } else {
                // setSize because optionsPanel is shrunk in dialog
                dialog.setSize(new java.awt.Dimension(priorSettingsPanel.getWidth() + 100, priorSettingsPanel.getHeight() + 100));
            }

//            System.out.println("panel width = " + panel.getWidth());
//            System.out.println("panel height = " + panel.getHeight());
        }

        dialog.pack();
        dialog.setResizable(true);
        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public void getArguments(Parameter parameter) {
            priorSettingsPanel.getArguments(parameter);
    }

    public boolean hasInvalidInput(boolean showError) {
        return priorSettingsPanel.hasInvalidInput(showError);
    }
}