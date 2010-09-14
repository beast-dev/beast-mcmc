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
