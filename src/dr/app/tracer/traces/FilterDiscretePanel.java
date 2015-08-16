/*
 * FilterDiscretePanel.java
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

package dr.app.tracer.traces;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.TreeSet;

/**
 * @author Walter Xie
 */
public class FilterDiscretePanel extends FilterAbstractPanel {
    JList allValues;
    JList selectedValues;
//        JButton selectButton;

    FilterDiscretePanel(TreeSet<String> allValuesSet, String[] selectedValuesArray) {
        JPanel panel = new JPanel(new FlowLayout(FlowLayout.CENTER));

        allValues = new JList(allValuesSet.toArray());
        allValues.setVisibleRowCount(8);
        allValues.setFixedCellWidth(80);
        allValues.setFixedCellHeight(20);
        allValues.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        panel.add(new JScrollPane(allValues));

//        if (selectedValuesArray != null) {
//            int[] indices = new int[selectedValuesArray.length];
//            for (int i = 0; i < indices.length; i++) {
//                for (int j = 0; j < allValuesSet.length; j++) {
//                    if (selectedValuesArray[i].equals(allValuesSet[j])) {
//                        indices[i] = j;
//                        break;
//                    }
//                }
//            }
//
//            allValues.setSelectedIndices(indices);
//        }

        JButton selectButton = new JButton("Select >>>");
        selectButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                selectedValues.setListData(allValues.getSelectedValues());
            }
        });
        panel.add(selectButton);

        if (selectedValuesArray == null) {
            selectedValues = new JList();
        } else {
            selectedValues = new JList(selectedValuesArray);
        }
        selectedValues.setVisibleRowCount(8);
        selectedValues.setFixedCellWidth(80);
        selectedValues.setFixedCellHeight(20);
        selectedValues.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
        panel.add(new JScrollPane(selectedValues));

        JPanel panel1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel1.add(new JLabel("<html>Hold Shift or Ctrl key for multi-selection,<br> and re-select values from left to correct" +
                "<br> the selection list on the left.</html>"));

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(panel, BorderLayout.CENTER);
        add(panel1, BorderLayout.SOUTH);
    }

    public String[] getSelectedValues() {
        int size = selectedValues.getModel().getSize();
        if (size < 0) return null;
        String[] sel = new String[size];
        for (int i=0; i < size; i++) {
            sel[i] = selectedValues.getModel().getElementAt(i).toString();
        }
        return sel;
    }

}

