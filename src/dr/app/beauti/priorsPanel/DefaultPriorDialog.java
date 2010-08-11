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

	public boolean showDialog(BeautiOptions options) {
        priorsPanel.setParametersList(options);

        Object[] buttons = {"Continue", "Back to BEAUTi"};
        JOptionPane optionPane = new JOptionPane(priorsPanel,
				JOptionPane.PLAIN_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				buttons,
				buttons[0]);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));
//        optionPane.setPreferredSize(new java.awt.Dimension(800, 600));

		final JDialog dialog = optionPane.createDialog(frame, "Check the default priors");
		dialog.pack();
		dialog.setVisible(true);
        dialog.setResizable(true);
        
        return optionPane.getValue() != null && optionPane.getValue().equals(buttons[0]);

    }
}