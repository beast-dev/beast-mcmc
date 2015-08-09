/*
 * MultiLineTableCellRenderer.java
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

package dr.app.gui.table;

import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.*;
import java.awt.*;

import jam.mac.Utils;

public class MultiLineTableCellRenderer extends DefaultTableCellRenderer {

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        JLabel label = (JLabel) super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
        if ((row % 2 == 0) && !isSelected) {
            Color newColour = alternateRowColor(table.getSelectionBackground());
            label.setBackground(newColour);
        } else if (!isSelected) {
            label.setBackground(table.getBackground());
        }
        return label;
    }

    public Color alternateRowColor(Color selectionColor) {
        Color useColor = selectionColor;
        if (Utils.isMacOSX()) {
            useColor = UIManager.getColor("Focus.color");
        }
        return new Color(255 - (255 - useColor.getRed()) / 10, 255 - (255 - useColor.getGreen()) / 10, 255 - (255 - useColor.getBlue()) / 10);

    }

    protected void setValue(Object value) {
        MultiLineTableCellContent content = (MultiLineTableCellContent) value;
        setIcon(content.getTableCellIcon());
        setText(content.getTableCellContent());
        setToolTipText(content.getToolTipContent());
    }
}
