/*
 * DateCellEditor.java
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

import dr.app.gui.components.RealNumberField;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


public class DateCellEditor extends DefaultCellEditor {

    /**
     *
     */
    private static final long serialVersionUID = 5067833373685886590L;
    private RealNumberField editor;

    public DateCellEditor() {
        this(false);
    }

    public DateCellEditor(boolean allowEmpty) {
        super(new RealNumberField(0.0, Double.MAX_VALUE));

        editor = (RealNumberField) getComponent();
        editor.setAllowEmpty(allowEmpty);

        setClickCountToStart(2); //This is usually 1 or 2.

        // Must do this so that editing stops when appropriate.
        editor.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                fireEditingStopped();
            }
        });
    }

    protected void fireEditingStopped() {
        super.fireEditingStopped();
    }

    public Object getCellEditorValue() {
        return editor.getValue();
    }

    public Component getTableCellEditorComponent(JTable table,
                                                 Object value,
                                                 boolean isSelected,
                                                 int row,
                                                 int column) {
        FontMetrics metrics = table.getFontMetrics(table.getFont());
        int fontHeight = metrics.getHeight();
        table.setRowHeight(row, fontHeight + fontHeight / 2);
//      System.out.println(editor.getPreferredSize() + "\t" + table.getRowHeight(row) + "\t" + table.getHeight());
        editor.setFont(table.getFont());
        if (value != null) {
            editor.setValue(((Double) value).doubleValue());
        }
        return editor;
    }
}