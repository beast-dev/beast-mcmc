package dr.app.tracer.traces;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * @author Walter Xie
 */
public class FilterDiscretePanel extends FilterAbstractPanel {
    JList allValues;
    JList selectedValues;
//        JButton selectButton;

    FilterDiscretePanel(String[] allValuesArray, String[] selectedValuesArray) {
        setLayout(new FlowLayout(FlowLayout.CENTER, 20, 50));

        allValues = new JList(allValuesArray);
        allValues.setVisibleRowCount(6);
        allValues.setFixedCellWidth(100);
        allValues.setFixedCellHeight(15);
        allValues.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        add(new JScrollPane(allValues));

        if (selectedValuesArray != null) {
            int[] indices = new int[selectedValuesArray.length];
            for (int i = 0; i < indices.length; i++) {
                for (int j = 0; j < allValuesArray.length; j++) {
                    if (selectedValuesArray[i].equals(allValuesArray[j])) {
                        indices[i] = j;
                        break;
                    }
                }
            }

            allValues.setSelectedIndices(indices);
        }

        JButton selectButton = new JButton("Select >>>");
        selectButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                selectedValues.setListData(allValues.getSelectedValues());
            }
        });
        add(selectButton);

        if (selectedValuesArray == null) {
            selectedValues = new JList();
        } else {
            selectedValues = new JList(selectedValuesArray);
        }
        selectedValues.setVisibleRowCount(6);
        selectedValues.setFixedCellWidth(100);
        selectedValues.setFixedCellHeight(15);
        selectedValues.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
        add(new JScrollPane(selectedValues));

        add(new JLabel("Hold Shift or Ctrl key for multi-selection"));
    }

    public Object[] getSelectedValues() {
        int size = selectedValues.getModel().getSize();
        if (size < 0) return null;
        Object[] sel = new Object[size];
        for (int i=0; i < size; i++) {
           sel[i] = selectedValues.getModel().getElementAt(i); 
        }
        return sel;
    }

}

