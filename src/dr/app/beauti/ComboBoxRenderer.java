/*
 * ComboBoxRenderer.java
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

package dr.app.beauti;

import javax.swing.*;
import javax.swing.table.TableCellRenderer;
import java.awt.*;

/**
 * @author Andrew Rambaut
* @version $Id$
*/
public class ComboBoxRenderer extends JComboBox implements TableCellRenderer {
    private final boolean isListAll;

    public ComboBoxRenderer() {
        super();
        setOpaque(true);
        isListAll = false;
        putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
    }

    public <T> ComboBoxRenderer(T[] allValues) {
        super(allValues);
        setOpaque(true);
        isListAll = true;
        putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
    }

    public Component getTableCellRendererComponent(JTable table, Object value,
                                                   boolean isSelected, boolean hasFocus,
                                                   int row, int column) {

        if (isSelected) {
            this.setForeground(table.getSelectionForeground());
            this.setBackground(table.getSelectionBackground());
        } else {
            this.setForeground(table.getForeground());
            this.setBackground(table.getBackground());
        }

        if (isListAll) {
            setSelectedItem(value);
        } else {
            if (value != null) {
                removeAllItems();
                addItem(value);
            }
        }
        return this;
    }

}
