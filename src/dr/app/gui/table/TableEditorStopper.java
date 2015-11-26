/*
 * TableEditorStopper.java
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

import java.awt.Component;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import javax.swing.JTable;
import javax.swing.table.TableCellEditor;

/**
 * This class is a workaround for a bug in Java (4503845)
 *      see http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4503845
 * Although the state appears as fixed, it does not yet work on 1.5
 *      (see the comments at the end on the previous page)
 *
 * The bug refers to tables with editing cells. When the user interacts
 *   with other UI element -or other cell in the table-, the editor component
 *   should receive a stopCellEditing call.
 *   This does not happen if the user selects a UI element outside the table.
 *
 * For example:
 *   - A dialog contains a table and a Cancel button that, before closing,
 *      verifies whether there are changes to save.
 *   -The user is modifying one cell and then presses the Cancel button
 *   -Then the cell editor does not receive an stopCellEditing call, so
 *      there seems to be no changes to save!
 *
 * A basic solution for this last example would be to call to
 *   table.getCellEditor().stopCellEditing() when reacting to events on the
 *   Cancel button, but this call should be made on every event handler for
 *   each UI element!
 *
 * This class works by listening to focus events when a cell editor is selected.
 * If the focus is transfered to any place outside the table, it calls to
 *  stopCellEditing()
 * Why is it not enough to check the focus on just the editor component?
 *   => Because the TableCellEditor could be composed of multiple Components
 *
 */
public class TableEditorStopper extends FocusAdapter implements PropertyChangeListener
{
    public static void ensureEditingStopWhenTableLosesFocus(JTable table)
    {
        new TableEditorStopper(table);
    }
    private TableEditorStopper(JTable table)
    {
        this.table=table;
        table.addPropertyChangeListener("tableCellEditor", this);
    }
    public void propertyChange(PropertyChangeEvent evt)
    {
        if (focused!=null)
        {
            focused.removeFocusListener(this);
        }
        focused = table.getEditorComponent();
        if (focused!=null)
        {
            focused.addFocusListener(this);
        }
    }
    public void focusLost(FocusEvent e)
    {
        if (focused!=null)
        {
            focused.removeFocusListener(this);
            focused = e.getOppositeComponent();
            if (table==focused || table.isAncestorOf(focused))
            {
                focused.addFocusListener(this);
            }
            else
            {
                focused=null;
                TableCellEditor editor = table.getCellEditor();
                if (editor!=null)
                {
                    editor.stopCellEditing();
                }
            }
        }
    }
    private Component focused;
    private JTable table;
}
