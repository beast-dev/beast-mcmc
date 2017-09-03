/*
 * PanelUtils.java
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

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class PanelUtils {
	public static JPanel createAddRemoveButtonPanel(Action addAction,
			Icon addIcon, String addToolTip, Action removeAction,
			Icon removeIcon, String removeToolTip, int axis) {

		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel, axis));
		buttonPanel.setOpaque(false);
		JButton addButton = new JButton(addAction);
		if (addIcon != null) {
			addButton.setIcon(addIcon);
			addButton.setText(null);
		}
		addButton.setToolTipText(addToolTip);
		addButton.putClientProperty("JButton.buttonType", "toolbar");
		addButton.setOpaque(false);
		addAction.setEnabled(false);

		JButton removeButton = new JButton(removeAction);
		if (removeIcon != null) {
			removeButton.setIcon(removeIcon);
			removeButton.setText(null);
		}
		removeButton.setToolTipText(removeToolTip);
		removeButton.putClientProperty("JButton.buttonType", "toolbar");
		removeButton.setOpaque(false);
		removeAction.setEnabled(false);

		buttonPanel.add(addButton);
		buttonPanel.add(new JToolBar.Separator(new Dimension(6, 6)));
		buttonPanel.add(removeButton);

		return buttonPanel;
	}

	public static void setupComponent(JComponent comp) {
		comp.setOpaque(false);

		// comp.setFont(UIManager.getFont("SmallSystemFont"));
		// comp.putClientProperty("JComponent.sizeVariant", "small");
		if (comp instanceof JButton) {
			comp.putClientProperty("JButton.buttonType", "textured");
//			comp.putClientProperty("JButton.buttonType", "bevel");
//			comp.putClientProperty("JButton.buttonType", "roundRect");
		}
		if (comp instanceof JComboBox) {
			comp.putClientProperty("JComboBox.isPopDown", Boolean.TRUE);
			// comp.putClientProperty("JComboBox.isSquare", Boolean.TRUE);
		}
        if (comp instanceof JTextArea) {
            ((JTextArea) comp).setEditable(false);
            ((JTextArea) comp).setOpaque(false);
            ((JTextArea) comp).setLineWrap(true);
            ((JTextArea) comp).setWrapStyleWord(true);
        }
	}

	public static Frame getActiveFrame() {
		Frame result = null;
		Frame[] frames = Frame.getFrames();
		for (int i = 0; i < frames.length; i++) {
			Frame frame = frames[i];
			if (frame.isVisible()) {
				result = frame;
				break;
			}
		}
		return result;
	}// END: getActiveFrame

}
