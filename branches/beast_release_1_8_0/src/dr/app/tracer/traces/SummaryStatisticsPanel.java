package dr.app.tracer.traces;

import dr.inference.trace.TraceAnalysis;
import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceFactory;
import dr.inference.trace.TraceList;
import jam.framework.Exportable;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;


public class SummaryStatisticsPanel extends JPanel implements Exportable {

    static final String NAME_ROW = "name";
    static final String MEAN_ROW = "mean";
    static final String MODE_ROW = "mode";
    static final String STDEV_ROW = "stderr of mean";
    static final String STDEV = "stdev";
    static final String FREQ_MODE_ROW = "frequency of mode";
    static final String VARIANCE_ROW = "variance";
    //    static final String STDEV_VAR_ROW = "stderr of variance";
    static final String GEOMETRIC_MEAN_ROW = "geometric mean";
    static final String MEDIAN_ROW = "median";
    static final String LOWER_UPPER_ROW = "95% HPD Interval";
    static final String CRED_SET_ROW = "95% Credible Set";
    static final String INCRED_SET_ROW = "5% Incredible Set";
    static final String ACT_ROW = "auto-correlation time (ACT)";
    static final String ESS_ROW = "effective sample size (ESS)";
    static final String SUM_ESS_ROW = "effective sample size (sum of ESS)";

    TraceList[] traceLists = null;
    java.util.List<String> traceNames = null;

    StatisticsModel statisticsModel;
    JTable statisticsTable = null;
    JScrollPane scrollPane1 = null;
    JPanel topPanel = null;
    JSplitPane splitPane1 = null;

    int dividerLocation = -1;

    FrequencyPanel frequencyPanel = null;
    IntervalsPanel intervalsPanel = null;
    JComponent currentPanel = null;

    public SummaryStatisticsPanel(JFrame frame) {

        setOpaque(false);

        statisticsModel = new StatisticsModel();
        statisticsTable = new JTable(statisticsModel);

        statisticsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.RIGHT, new Insets(0, 4, 0, 4)));
        statisticsTable.getColumnModel().getColumn(1).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        statisticsTable.addMouseMotionListener(new MouseMotionAdapter(){
            public void mouseMoved(MouseEvent e){
                Point p = e.getPoint();
                int row = statisticsTable.rowAtPoint(p);
                if (row == 6) {
                    Object t = statisticsModel.getValueAt(9,1);
                    if (!t.equals("-")) {
                        statisticsTable.setToolTipText("Incredible set : " + t);
                    }
                }
            }//end MouseMoved
        }); // end MouseMotionAdapter

        scrollPane1 = new JScrollPane(statisticsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

        topPanel = new JPanel(new BorderLayout(0, 0));
        topPanel.setOpaque(false);
        topPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
                new java.awt.Insets(0, 0, 6, 0)));
        topPanel.add(scrollPane1, BorderLayout.CENTER);

        frequencyPanel = new FrequencyPanel(frame);
        frequencyPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
                new java.awt.Insets(6, 0, 0, 0)));

        intervalsPanel = new IntervalsPanel();
        intervalsPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
                new java.awt.Insets(6, 0, 0, 0)));

        splitPane1 = new JSplitPane(JSplitPane.VERTICAL_SPLIT, true, topPanel, frequencyPanel);
        splitPane1.setOpaque(false);
        splitPane1.setBorder(null);
//		splitPane1.setBorder(new BorderUIResource.EmptyBorderUIResource(
//								new java.awt.Insets(12, 12, 12, 12)));

        setLayout(new BorderLayout(0, 0));
        add(splitPane1, BorderLayout.CENTER);

        splitPane1.setDividerLocation(2000);

    }

    private void setupDividerLocation() {

        if (dividerLocation == -1 || dividerLocation == splitPane1.getDividerLocation()) {
            int h0 = topPanel.getHeight();
            int h1 = scrollPane1.getViewport().getHeight();
            int h2 = statisticsTable.getPreferredSize().height;
            dividerLocation = h2 + h0 - h1;

            splitPane1.setDividerLocation(dividerLocation);
        }
    }

    public void setTraces(TraceList[] traceLists, java.util.List<String> traceNames) {

        this.traceLists = traceLists;
        this.traceNames = traceNames;

        int divider = splitPane1.getDividerLocation();

        statisticsModel.fireTableStructureChanged();
        if (traceLists != null && traceNames != null) {
            if (traceLists.length == 1 && traceNames.size() == 1) {
                statisticsTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

                currentPanel = frequencyPanel;
                frequencyPanel.setTrace(traceLists[0], traceNames.get(0));
                intervalsPanel.setTraces(null, null);
                splitPane1.setBottomComponent(frequencyPanel);
            } else {
                statisticsTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
                statisticsTable.getColumnModel().getColumn(0).setPreferredWidth(200);
                for (int i = 1; i < statisticsTable.getColumnCount(); i++) {
                    statisticsTable.getColumnModel().getColumn(i).setPreferredWidth(100);
                }

                currentPanel = intervalsPanel;
                frequencyPanel.setTrace(null, null);
                intervalsPanel.setTraces(traceLists, traceNames);
                splitPane1.setBottomComponent(intervalsPanel);
            }
        } else {
            currentPanel = statisticsTable;
            frequencyPanel.setTrace(null, null);
            splitPane1.setBottomComponent(frequencyPanel);
        }

        splitPane1.setDividerLocation(divider);

        statisticsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.RIGHT, new Insets(0, 4, 0, 4)));
        for (int i = 1; i < statisticsTable.getColumnCount(); i++) {
            statisticsTable.getColumnModel().getColumn(i).setCellRenderer(
                    new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        }

        setupDividerLocation();

        validate();
        repaint();
    }

    public JComponent getExportableComponent() {
        if (currentPanel instanceof Exportable) {
            return ((Exportable) currentPanel).getExportableComponent();
        }
        return currentPanel;
    }

    public String toString() {
        return statisticsModel.toString();
    }

    class StatisticsModel extends AbstractTableModel {

        String[] rowNames = {MEAN_ROW, STDEV_ROW, STDEV, VARIANCE_ROW, MEDIAN_ROW, MODE_ROW, GEOMETRIC_MEAN_ROW,
                LOWER_UPPER_ROW, ACT_ROW, ESS_ROW};

        public StatisticsModel() {
        }

        public int getColumnCount() {
            if (traceLists != null && traceNames != null) {
                return (traceLists.length * traceNames.size()) + 1;
            } else {
                return 2;
            }
        }

        public int getRowCount() {
            return rowNames.length;
        }

        public Object getValueAt(int row, int col) {

            TraceCorrelation tc = null;
            if (traceLists != null && traceNames != null && traceNames.size() > 0) {
                int colNew = col;
                if (col == 0) {
                    colNew++;
                }
                int n1 = (colNew - 1) / traceNames.size();
                int n2 = (colNew - 1) % traceNames.size();

                TraceList tl = traceLists[n1];
                int index = tl.getTraceIndex(traceNames.get(n2));
                tc = tl.getCorrelationStatistics(index);
            } else {
                return "-";
            }

            if (col == 0) {
                if (tc != null && tc.getTraceType() != TraceFactory.TraceType.DOUBLE && row == 6) {
                    return CRED_SET_ROW;
//                } else if (tc != null && tc.getTraceType() != TraceFactory.TraceType.DOUBLE && row == 6) {
//                    return INCRED_SET_ROW;
                } else {
                    return rowNames[row];
                }
            }

            double value = 0.0;

            if (tc != null) {
                if (row != 0 && !tc.isValid()) return "n/a";

                switch (row) {
                    case 0:
                        if (tc.getTraceType() == TraceFactory.TraceType.STRING) {
                            return "n/a";
                        } else {
                            value = tc.getMean();
                        }
                        break;
                    case 1:
                        if (tc.getTraceType() == TraceFactory.TraceType.STRING) {
                            return "n/a";
                        } else {
                            value = tc.getStdErrorOfMean();
                        }
                        break;
                    case 2:
                        if (tc.getTraceType() == TraceFactory.TraceType.STRING) {
                            return "n/a";
                        } else {
                            value = tc.getStdError();
                        }
                        break;
                    case 3:
                        if (tc.getTraceType() == TraceFactory.TraceType.STRING) {
                            return "n/a";
                        } else {
                            value = tc.getVariance();
                        }
                        break;
                    case 4:
                        if (tc.getTraceType() == TraceFactory.TraceType.STRING) {
                            return "n/a";
                        } else {
                            value = tc.getMedian();
                        }
                        break;
                    case 5:
                        if (tc.getTraceType() == TraceFactory.TraceType.DOUBLE) {
                            return "n/a";
                        } else {
                            return tc.getMode();
                        }
                    case 6:
                        if (!tc.hasGeometricMean()) return "n/a";
                        value = tc.getGeometricMean();
                        break;
                    case 7:
                        if (tc.getTraceType() == TraceFactory.TraceType.DOUBLE) {
                            return "[" + TraceAnalysis.formattedNumber(tc.getLowerHPD()) + ", " + TraceAnalysis.formattedNumber(tc.getUpperHPD()) + "]";
                        } else {
                            return tc.printCredibleSet();
                        }
                    case 8:
                        value = tc.getACT();
                        break;
                    case 9:
                        value = tc.getESS();
                        break;
                    case 10:
                        if (tc.getTraceType() == TraceFactory.TraceType.DOUBLE) {
                            return "-";
                        } else {
                            return tc.printInCredibleSet();
                        }

                }
            } else {
                return "-";
            }

            return TraceAnalysis.formattedNumber(value);
        }

        public String getColumnName(int column) {
            if (column == 0) return "Summary Statistic";
            if (traceLists != null && traceNames != null) {
                int n1 = (column - 1) / traceNames.size();
                int n2 = (column - 1) % traceNames.size();
                String columnName = "";
                if (traceLists.length > 1) {
                    columnName += traceLists[n1].getName();
                    if (traceNames.size() > 1) {
                        columnName += ": ";
                    }
                }
                if (traceNames.size() > 1) {
                    columnName += traceNames.get(n2);
                }
                return columnName;
            }
            return "-";
        }

        public Class getColumnClass(int c) {
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
