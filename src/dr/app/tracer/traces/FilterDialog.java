package dr.app.tracer.traces;

import dr.inference.trace.*;
import jam.framework.DocumentFrame;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class FilterDialog {
    private DocumentFrame frame;
    private static final int MINIMUM_TABLE_WIDTH = 160;

    TraceList selectedTraceList;
    String currentTraceName = null;
    Map<String, FilterPanel> filterPanels = new HashMap<String, FilterPanel>();
    JPanel panelParent;
    TitledBorder panelBorder;

    JTable traceFilterTable = null;
    TraceFilterTableModel traceFilterTableModel = null;

    public FilterDialog(DocumentFrame frame) {
        this.frame = frame;

        traceFilterTableModel = new TraceFilterTableModel();
        traceFilterTable = new JTable(traceFilterTableModel);

        TableRenderer renderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
        traceFilterTable.getColumnModel().getColumn(0).setCellRenderer(renderer);

        traceFilterTable.getColumnModel().getColumn(1).setPreferredWidth(20);

        traceFilterTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
//        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(traceFilterTable);

        traceFilterTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        setCurrentFilter(null);
    }

    public String showDialog(TraceList selectedTraceList, String previousMessage) {
        this.selectedTraceList = selectedTraceList;
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

        JScrollPane scrollPane = new JScrollPane(traceFilterTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setOpaque(false);

        JTextField fileFiled = new JTextField(selectedTraceList.getName());
        fileFiled.setColumns(30);
        fileFiled.setEditable(false);
        JPanel tittlePanel = new JPanel(new BorderLayout(12, 12));
        tittlePanel.add(new JLabel("Selected Log File : "), BorderLayout.WEST);
        tittlePanel.add(fileFiled, BorderLayout.CENTER);

        panelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panelParent.setOpaque(false);
        panelBorder = new TitledBorder("");
        panelParent.setBorder(panelBorder);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, scrollPane, panelParent);
        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        JPanel basePanel = new JPanel(new BorderLayout(20, 20));
        basePanel.add(tittlePanel, BorderLayout.NORTH);
        basePanel.add(splitPane, BorderLayout.CENTER);

        Object[] options = {"Apply Filter", "Update Filter Changes", "Remove All Filters", "Cancel"};
        JOptionPane optionPane = new JOptionPane(basePanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.YES_NO_CANCEL_OPTION,
                null,
                options,
                options[0]);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Filter Summary");
        dialog.setModal(true);
        dialog.setResizable(true);
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

    private void selectionChanged() {
        int selRow = traceFilterTable.getSelectedRow();
        if (selRow >= 0) {
            setCurrentFilter(traceFilterTable.getValueAt(selRow, 0).toString());
            traceFilterTableModel.fireTableDataChanged();
        }
    }

    private void setCurrentFilter(String traceName) {
        if (traceName != null) {
            if (currentTraceName != null) panelParent.removeAll();

            FilterPanel panel = filterPanels.get(traceName);
            if (panel == null) {

                int traceIndex = selectedTraceList.getTraceIndex(traceName);
                TraceDistribution td = selectedTraceList.getDistributionStatistics(traceIndex);
                Filter f = ((FilteredTraceList) selectedTraceList).getFilter(traceIndex);

                String[] sel;
                if (f == null) {
                    sel = null;
                } else {
                    sel = f.getIn();
                }

                if (td.getTraceType() == TraceFactory.TraceType.DOUBLE) {
                    String[] minMax = new String[]{Double.toString(td.getMinimum()), Double.toString(td.getMaximum())};
                    panel = new FilterContinuousPanel(minMax, sel);
                } else {// integer and string
                    List<String> allNames = td.getRange();
                    String[] all = allNames.toArray(new String[allNames.size()]);
                    panel = new FilterDiscretePanel(all, sel);
                }

                filterPanels.put(traceName, panel);
            }

            currentTraceName = traceName;
            panelParent.add(panel);

            updateBorder();
        }
    }

    private void updateBorder() {
        String title;
        panelBorder.setTitle("");
        panelParent.repaint();
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

    class TraceFilterTableModel extends AbstractTableModel {

        String[] columnNames = {"Trace Name", "Filter"};

        public TraceFilterTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (selectedTraceList == null) return 0;
            return selectedTraceList.getTraceCount();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return selectedTraceList.getTraceName(row);
                case 1:
                    return !(selectedTraceList.getTrace(row).getFilter() == null);
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 1) return true;
            return false;
        }

        public void setValueAt(Object value, int row, int col) {
            // fire add filter function
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            if (getRowCount() == 0) {
                return Object.class;
            }
            return getValueAt(0, c).getClass();
        }

        public String toString() {
            StringBuffer buffer = new StringBuffer();

            buffer.append(getColumnName(0));
            for (int j = 1; j < getColumnCount(); j++) {
                buffer.append("\t");
                buffer.append(getColumnName(j));
            }
            buffer.append("\n");

            for (int i = 0; i < getRowCount(); i++) {
                buffer.append(getValueAt(i, 0));
                for (int j = 1; j < getColumnCount(); j++) {
                    buffer.append("\t");
                    buffer.append(getValueAt(i, j));
                }
                buffer.append("\n");
            }

            return buffer.toString();
        }
    }
}