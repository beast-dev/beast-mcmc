package dr.app.tracer.traces;

import dr.inference.trace.TraceList;
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
    Object[] options = {"Apply Filter", "Update Filter Changes", "Remove All Filters", "Cancel"};
            
    public FilterDialog(DocumentFrame frame) {
        this.frame = frame;


    }

    public String showDialog(TraceList selectedTraceList, String previousMessage) {

//        this.traceName = traceName;
        String message = "";
//        initComponents(filteredTraceListGroup.get(0));
//        FilteredTraceList filteredTraceList = (FilteredTraceList) filteredTraceListGroup.get(fileList.getSelectedIndex()); // only pick up the 1st one
//        int traceIndex = filteredTraceList.getTraceIndex(traceName);
//        TraceDistribution td = filteredTraceList.getDistributionStatistics(traceIndex);
//
//        typeField.setText(td.getTraceType().toString());
//        nameField.setText(traceName);

//        Filter f = filteredTraceList.getFilter(traceIndex);
//
//        String[] sel;
//        if (f == null) {
//            sel = null;
//        } else {
//            sel = f.getIn();
//        }
//        if (td.getTraceType() == TraceFactory.TraceType.DOUBLE) {
//            String[] minMax = new String[]{Double.toString(td.getMinimum()), Double.toString(td.getMaximum())};
//            filterPanel = filterPanel.new FilterContinuousPanel(minMax, sel);
//        } else {// integer and string
//            List<String> allNames = td.getRange();
//            String[] all = allNames.toArray(new String[allNames.size()]);
//            filterPanel = filterPanel.new FilterDiscretePanel(all, sel);
//        }

        FilterListPanel filterListPanel = new FilterListPanel(selectedTraceList, this);

        JOptionPane optionPane = new JOptionPane(filterListPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.YES_NO_CANCEL_OPTION,
                null,
                options,
                options[0]);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Filter Summary");
//        dialog.setModal(true);
//        dialog.setResizable(true);
        dialog.pack();
        dialog.setVisible(true);



        Object result = optionPane.getValue();
        if (result == null) return previousMessage;
        
//        FilteredTraceList filteredTraceList = filteredTraceListGroup.get(treeFileCombo.getSelectedIndex());
//        TraceDistribution td = filteredTraceList.getDistributionStatistics(traceIndex);

        if (result.equals(options[0])) {
             message = "";

        } else if (result.equals(options[1])) {

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
        } else if (result.equals(options[2])) {
//            for (int i = 0; i < filteredTraceListGroup.size(); i++) {
//                ((FilteredTraceList) filteredTraceListGroup.get(i)).removeFilter(traceIndex);
//            }
            message = "";

        } else if (result.equals(options[3])) {
            return previousMessage;
        }

        return message;
    }


//    private void initComponents(FilteredTraceList filteredTraceList) {
//
//        TraceDistribution td = filteredTraceList.getDistributionStatistics(filteredTraceList.getTraceIndex(traceName));
//
//        typeField.setText(td.getTraceType().toString());
//        nameField.setText(traceName);
//
//        Filter f = filteredTraceList.getFilter(traceName);
//
//        if (td.getTraceType() == TraceFactory.TraceType.DOUBLE) {
//
//        } else {// integer and string
//            String[] all = td.getRangeAll();
//            String[] sel;
//
//            if (f == null) {
//                sel = null;
//            } else {
//                sel = f.getIn();
//            }
//
//            filterPanel = new FilterDiscretePanel(all, sel);
//        }
//
//
//    }
//
//    public String getName() {
//        return nameField.getText();
//    }


}