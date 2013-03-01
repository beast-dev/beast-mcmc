package dr.app.tracer.traces;

import dr.app.gui.util.LongTask;
import dr.inference.trace.Filter;
import dr.inference.trace.FilteredTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceFactory;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * @author Walter Xie
 */
public class FilterListPanel extends JPanel {
    private static final int MINIMUM_TABLE_WIDTH = 200;
    private static final int PREFERRED_HIGHT = 400;

    FilteredTraceList selectedTraceList;
    String currentTraceName = null;
    Map<String, FilterAbstractPanel> filterPanels = new HashMap<String, FilterAbstractPanel>();
    TitledBorder panelBorder;
    JPanel modelPanelParent;

    JTable traceFilterTable = null;

    public FilterListPanel(FilteredTraceList selectedTraceList) {
        this.selectedTraceList = selectedTraceList;

        TraceFilterTableModel traceFilterTableModel = new TraceFilterTableModel();
        traceFilterTable = new JTable(traceFilterTableModel);

        traceFilterTable.getTableHeader().setReorderingAllowed(false);
        traceFilterTable.setColumnSelectionAllowed(false);

        TableRenderer renderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
        traceFilterTable.getColumnModel().getColumn(0).setCellRenderer(renderer);

        traceFilterTable.getColumnModel().getColumn(1).setPreferredWidth(20);

        traceFilterTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
//        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(traceFilterTable);

        traceFilterTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });


        JTextField fileFiled = new JTextField(selectedTraceList.getName());
        fileFiled.setColumns(30);
        fileFiled.setEditable(false);
        JPanel tittlePanel = new JPanel(new BorderLayout(12, 12));
        tittlePanel.add(new JLabel("Selected Log File : "), BorderLayout.WEST);
        tittlePanel.add(fileFiled, BorderLayout.CENTER);

        modelPanelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        modelPanelParent.setOpaque(false);
        panelBorder = new TitledBorder("Select a trace");
        modelPanelParent.setBorder(panelBorder);

        JScrollPane scrollPane = new JScrollPane(traceFilterTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setOpaque(false);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);
        scrollPane.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, PREFERRED_HIGHT));

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, scrollPane, modelPanelParent);
        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setPreferredSize(new Dimension(600, PREFERRED_HIGHT));
        setLayout(new BorderLayout(20, 20));
        add(tittlePanel, BorderLayout.NORTH);
        add(splitPane, BorderLayout.CENTER);

        final LoadFiltersTask analyseTask = new LoadFiltersTask(selectedTraceList);

        final ProgressMonitor progressMonitor = new ProgressMonitor(this,
                "Loading filters ...",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

        analyseTask.go();

//        Runnable r = new Runnable() {
//            public void run() {
//                initFilterPanels();
//            }
//        };
//        if(SwingUtilities.isEventDispatchThread()) {
//            r.run();
//        }
//        else {
//            SwingUtilities.invokeLater(r);
//        }
    }

    class LoadFiltersTask extends LongTask {
        private int lengthOfTask = 0;
        private int current = 0;

        public LoadFiltersTask(FilteredTraceList selectedTraceList) {
            lengthOfTask = selectedTraceList.getTraceCount();
        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Loading filters ...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {
            for (int i = 0; i < selectedTraceList.getTraceCount(); i++) {
                String traceName = selectedTraceList.getTraceName(i);
                Trace trace = selectedTraceList.getTrace(i);
                //            TraceDistribution td = selectedTraceList.getDistributionStatistics(i);
                Filter f = selectedTraceList.getFilter(i);

                String[] sel;
                if (f == null) {
                    sel = null;
                } else {
                    sel = f.getIn();
                }

                FilterAbstractPanel panel;
                if (trace.getTraceType() == TraceFactory.TraceType.STRING) {
                    panel = new FilterDiscretePanel(trace.getRange(), sel);
                } else {// integer and double
                    panel = new FilterContinuousPanel(trace.getRange(), sel);
                }
                //            System.out.println("traceName = " + traceName + ";  i = " + i);
                filterPanels.put(traceName, panel);

                current += 1;
            }

            return null;
        }

    }

    private void selectionChanged() {
        int selRow = traceFilterTable.getSelectedRow();
//        System.out.println("selRow = " + selRow);
        if (selRow >= 0) {
            String traceName = selectedTraceList.getTraceName(selRow);
            setCurrentFilter(traceName);
        }
    }

    private void setCurrentFilter(String traceName) {
        if (traceName != null) {
            modelPanelParent.removeAll();

            FilterAbstractPanel panel = filterPanels.get(traceName);

            if (panel == null) {
                throw new RuntimeException("null filter panel, " + traceName);
            }

            panelBorder.setTitle(traceName);
            modelPanelParent.add(panel);

            validate();
            repaint();
        }
    }

    public void applyFilterChanges() {
        Iterator iterator = filterPanels.keySet().iterator();

        while (iterator.hasNext()) {
            String traceName = iterator.next().toString();
            FilterAbstractPanel fp = filterPanels.get(traceName);

            int traceIndex = selectedTraceList.getTraceIndex(traceName);
//            System.out.println("traceName = " + traceName + ";  traceIndex = " + traceIndex);
            if (fp.containsNullValue()) {
                if (selectedTraceList.hasFilter(traceIndex))
                    selectedTraceList.removeFilter(traceIndex);
            } else {
//                System.out.println("traceIndex = " + traceIndex + "; fp.getSelectedValues() " + fp.getSelectedValues().length);
//                for (Object o : fp.getSelectedValues()) {
//                    System.out.print("; " + o);
//                }
//                System.out.println();

                Filter f;
                String[] in = fp.getSelectedValues();
//                if (selectedTraceList.getTrace(traceIndex).getTraceType() == TraceFactory.TraceType.INTEGER) {
//                    // as Integer is stored as Double in Trace
//                    Object[] inInt = new Object[in.length];
//                    for (int i = 0; i < in.length; i++) {
//                        inInt[i] = Double.valueOf(in[i]);
//                    }
//                    f = new Filter(inInt, selectedTraceList.getTrace(traceIndex).getValuesSize());
//                } else {
                    f = new Filter(in, selectedTraceList.getTrace(traceIndex).getTraceType());
//                }
                selectedTraceList.setFilter(traceIndex, f);
            }
        }
    }

    public void removeAllFilters() {
        selectedTraceList.removeAllFilters();
    }

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
                    return selectedTraceList.hasFilter(row);
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 1) return true;
            return false;
        }

        public void setValueAt(Object value, int row, int col) {
//            if (col == 1 && (Boolean) value) {
//                System.out.println("row = " + row);
//            }
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
