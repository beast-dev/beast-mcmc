package dr.app.gui.table;

import dr.app.gui.components.RealNumberField;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


public class RealNumberCellEditor extends DefaultCellEditor {

    private RealNumberField editor;

    public RealNumberCellEditor(double minValue, double maxValue) {
        super(new RealNumberField(minValue, maxValue));

        editor = (RealNumberField) getComponent();

        editor.setBorder(BorderFactory.createLineBorder(Color.black));

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
        editor.setValue(((Double) value).doubleValue());
        return editor;
    }
}