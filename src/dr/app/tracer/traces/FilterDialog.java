/*
 * FilterDialog.java
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

import dr.app.gui.util.LongTask;
import dr.app.tracer.application.TracerFrame;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class FilterDialog extends JDialog {

    private TracerFrame frame;
    Object[] options = {"Apply Filter Changes", "Remove All Filters", "Cancel"};

    public FilterDialog(TracerFrame frame) {
        this.frame = frame;
    }

    public String showDialog(FilterListPanel filterListPanel, String previousMessage) {
        String message = "";

        JOptionPane optionPane = new JOptionPane(filterListPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.YES_NO_CANCEL_OPTION,
                null,
                options,
                options[0]);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Filter Editor");
//        dialog.setModal(true);
//        dialog.setResizable(true);
        dialog.pack();
        dialog.setVisible(true);

        Object result = optionPane.getValue();
        if (result == null) return previousMessage;

        if (result.equals(options[0])) {
            message = "";

            filterListPanel.applyFilterChanges();
            frame.statisticTableSelectionChanged();

//            if (filterPanel.containsNullValue()) {
//                JOptionPane.showMessageDialog(frame, "The selected value for filter is invalid \ror no value is selected !",
//                        "Invalid Filter Input",
//                        JOptionPane.ERROR_MESSAGE);
//                return previousMessage;
//            }


//            for (int i = 0; i < filteredTraceListGroup.size(); i++) {
//                FilteredTraceList fTL = (FilteredTraceList) filteredTraceListGroup.get(i);
//                f = fTL.getFilter(traceIndex);
//
//                if (f == null) {
//                    f = new Filter(filterPanel.getSelectedValues());
//                } else {
//                    f.setIn(filterPanel.getSelectedValues());
//                }
//
//                fTL.setFilter(traceIndex, f);
//            }
//            message += f.getStatusMessage(); // todo
//
//            for (int i = 0; i < filteredTraceListGroup.size(); i++) {
//                if (i == 0) message += " in file(s) ";
//                if (i > 0) message += " and ";
//                message += "\'" + filteredTraceListGroup.get(i).getName() + "\'";
//            }
        } else if (result.equals(options[1])) {
            filterListPanel.removeAllFilters();
            frame.statisticTableSelectionChanged();
            message = "";

        } else if (result.equals(options[2])) {
            return previousMessage;
        }

        return message;
    }

    class AnalyseFilterTask extends LongTask {

        public Object doWork() throws Exception {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getLengthOfTask() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getCurrent() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public String getMessage() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}