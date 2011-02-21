package dr.app.tracer.traces;

import dr.app.gui.table.TableEditorStopper;
import dr.inference.trace.*;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Walter Xie
 */
public class FilterListPanel extends JPanel {
    private static final int MINIMUM_TABLE_WIDTH = 160;

    TraceList selectedTraceList;
    String currentTraceName = null;
    Map<String, FilterAbstractPanel> filterPanels = new HashMap<String, FilterAbstractPanel>();
    JPanel panelParent;
    TitledBorder panelBorder;

    JTable traceFilterTable = null;
    TraceFilterTableModel traceFilterTableModel = null;

    final JDialog dialog; // use for repanit

    public FilterListPanel(TraceList selectedTraceList, JDialog dialog) {
        this.selectedTraceList = selectedTraceList;
        this.dialog = dialog;

        traceFilterTableModel = new TraceFilterTableModel();
        traceFilterTable = new JTable(traceFilterTableModel);

        TableRenderer renderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
        traceFilterTable.getColumnModel().getColumn(0).setCellRenderer(renderer);

        traceFilterTable.getColumnModel().getColumn(1).setPreferredWidth(20);

        traceFilterTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(traceFilterTable);

        traceFilterTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

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
        panelBorder = new TitledBorder("Select a trace");
        panelParent.setBorder(panelBorder);
        
        setCurrentFilter(null);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, scrollPane, panelParent);
        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setLayout(new BorderLayout(20, 20));
        add(tittlePanel, BorderLayout.NORTH);
        add(splitPane, BorderLayout.CENTER);
    }

    private void selectionChanged() {
        int selRow = traceFilterTable.getSelectedRow();
        System.out.println("selRow = " + selRow);
        if (selRow >= 0) {
            String traceName = selectedTraceList.getTraceName(selRow);
            setCurrentFilter(traceName);
        }
    }

    private void setCurrentFilter(String traceName) {
        if (traceName != null) {
            panelParent.removeAll();

            FilterAbstractPanel panel = filterPanels.get(traceName);
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
                    java.util.List<String> allNames = td.getRange();
                    String[] all = allNames.toArray(new String[allNames.size()]);
                    panel = new FilterDiscretePanel(all, sel);
                }

                filterPanels.put(traceName, panel);
            }

            currentTraceName = traceName;
            panelParent.add(panel);

            updateBorder(traceName);
//            panelParent.repaint();
            repaint(); //todo why not working?
            dialog.repaint();
        }
    }

    private void updateBorder(String traceName) {
        String title;
        if (traceName != null)  panelBorder.setTitle(traceName);
//        panelParent.repaint();
//        traceFilterTable.repaint();
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
            if (col == 1 && (Boolean) value) {
                System.out.println("row = " + row);
            }
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
