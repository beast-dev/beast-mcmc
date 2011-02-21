package dr.app.tracer.traces;

import dr.inference.trace.FilteredTraceList;
import jam.framework.DocumentFrame;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class FilterDialog extends JDialog {

    private DocumentFrame frame;
    Object[] options = {"Apply Filter Changes", "Remove All Filters", "Cancel"};
            
    public FilterDialog(DocumentFrame frame) {
        this.frame = frame;


    }

    public String showDialog(FilteredTraceList selectedTraceList, String previousMessage) {
        String message = "";

        FilterListPanel filterListPanel = new FilterListPanel(selectedTraceList, this);

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
            message = "";

        } else if (result.equals(options[2])) {
            return previousMessage;
        }

        return message;
    }

}