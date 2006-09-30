/*
 * TableRenderer.java
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

package org.virion.jam.table;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.DefaultTableCellRenderer;
import java.awt.*;


public class TableRenderer extends DefaultTableCellRenderer {
    protected Color bg1 = new Color(0xED, 0xF3, 0xFE);
    protected Color bg2 = Color.white;
    protected boolean striped;

    public TableRenderer(int alignment, Insets insets) {

        this(true, alignment, insets);
    }

    public TableRenderer(boolean striped, int alignment, Insets insets) {
        super();
        this.striped = striped;
        setOpaque(true);
        setHorizontalAlignment(alignment);
        if (insets != null) {
            setBorder(new BorderUIResource.EmptyBorderUIResource(insets));
        }
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                   boolean hasFocus, int row, int column) {

        if (value != null) {
            setText(value.toString());
        }
        setEnabled(table.isEnabled());
        setFont(table.getFont());

        // if cell is selected, set background color to default cell selection background color
        if (isSelected) {
            setBackground(table.getSelectionBackground());
            setForeground(table.getSelectionForeground());
        } else {
            if (striped) {
                if (row % 2 == 0) {
                    setBackground(bg1);
                } else {
                    setBackground(bg2);
                }
            } else {
                setBackground(table.getBackground());
            }
            setForeground(table.getForeground());
        }

        return this;
    }
}
