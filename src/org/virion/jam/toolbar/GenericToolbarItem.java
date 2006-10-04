/*
 * GenericToolbarItem.java
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

package org.virion.jam.toolbar;

import javax.swing.*;
import java.awt.*;

/**
 * @author rambaut
 *         Date: Oct 18, 2005
 *         Time: 10:09:21 PM
 */
public class GenericToolbarItem extends JPanel implements ToolbarItem {

    public GenericToolbarItem(String title, String toolTipText, JComponent component) {
        setLayout(new BorderLayout());
        add(component, BorderLayout.NORTH);

        label = new JLabel(title);
        label.setHorizontalAlignment(SwingConstants.CENTER);
        add(label, BorderLayout.SOUTH);
        setToolTipText(toolTipText);
    }

    public void setToolbarOptions(ToolbarOptions options) {
        switch (options.getDisplay()) {
            case ToolbarOptions.ICON_AND_TEXT:
            case ToolbarOptions.TEXT_ONLY:
                label.setVisible(true);
                break;
            case ToolbarOptions.ICON_ONLY:
                label.setVisible(false);
                break;
        }
    }

    public void setAction(Action action) {
        throw new UnsupportedOperationException("Method setAction() not supported in GenericToolBarItem");
    }

    private JLabel label;
}
