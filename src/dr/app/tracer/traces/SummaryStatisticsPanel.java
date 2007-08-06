package dr.app.tracer.traces;

import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;
import org.virion.jam.framework.Exportable;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.text.DecimalFormat;


public class SummaryStatisticsPanel extends JPanel implements Exportable {

    static final String NAME_ROW = "name";
    static final String MEAN_ROW = "mean";
    static final String STDEV_ROW = "stdev of mean";
    static final String MEDIAN_ROW = "median";
    static final String LOWER_ROW = "95% HPD lower";
    static final String UPPER_ROW = "95% HPD upper";
    static final String ACT_ROW = "auto-correlation time (ACT)";
    static final String ESS_ROW = "effective sample size (ESS)";
    static final String SUM_ESS_ROW = "effective sample size (sum of ESS)";

    TraceList[] traceLists = null;
    int[] traceIndices = null;

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

        currentPanel = statisticsTable;
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

    public void setTraces(TraceList[] traceLists, int[] traces) {

        this.traceLists = traceLists;
        this.traceIndices = traces;

        int divider = splitPane1.getDividerLocation();

        statisticsModel.fireTableStructureChanged();
        if (traceLists != null && traces != null) {
            if (traceLists.length == 1 && traces.length == 1) {
                statisticsTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

                currentPanel = frequencyPanel;
                frequencyPanel.setTrace(traceLists[0], traces[0]);
                intervalsPanel.setTraces(null, null);
                splitPane1.setBottomComponent(frequencyPanel);
            } else {
                statisticsTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
                statisticsTable.getColumnModel().getColumn(0).setPreferredWidth(200);
                for (int i = 1; i < statisticsTable.getColumnCount(); i++) {
                    statisticsTable.getColumnModel().getColumn(i).setPreferredWidth(100);
                }

                currentPanel = intervalsPanel;
                frequencyPanel.setTrace(null, -1);
                intervalsPanel.setTraces(traceLists, traces);
                splitPane1.setBottomComponent(intervalsPanel);
            }
        } else {
            currentPanel = statisticsTable;
            frequencyPanel.setTrace(null, -1);
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

    public void copyToClipboard() {
        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(statisticsModel.toString());

        clipboard.setContents(selection, selection);
    }

    public JComponent getExportableComponent() {
        return currentPanel;
    }

    public String toString() {
        return statisticsModel.toString();
    }

    class StatisticsModel extends AbstractTableModel {

        String[] rowNames = {MEAN_ROW, STDEV_ROW, MEDIAN_ROW, LOWER_ROW, UPPER_ROW, ACT_ROW, ESS_ROW};

        private DecimalFormat formatter = new DecimalFormat("0.###E0");
        private DecimalFormat formatter2 = new DecimalFormat("####0.###");

        public StatisticsModel() {
        }

        public int getColumnCount() {
            if (traceLists != null && traceIndices != null) {
                return (traceLists.length * traceIndices.length) + 1;
            } else {
                return 2;
            }
        }

        public int getRowCount() {
            return rowNames.length;
        }

        public Object getValueAt(int row, int col) {

            if (col == 0) {
                return rowNames[row];
            }

            TraceDistribution td;
            TraceCorrelation tc = null;

            if (traceLists != null && traceIndices != null) {
                int n1 = (col - 1) / traceIndices.length;
                int n2 = (col - 1) % traceIndices.length;

                TraceList tl = traceLists[n1];
                if (tl instanceof CombinedTraces) {
                    td = tl.getDistributionStatistics(traceIndices[n2]);
                } else {
                    tc = tl.getCorrelationStatistics(traceIndices[n2]);
                    td = tc;
                }
            } else {
                return "-";
            }

            double value = 0.0;

            if (tc != null) {
                if (row != 0 && !tc.isValid()) return "n/a";

                switch (row) {
                    case 0:
                        value = tc.getMean();
                        break;
                    case 1:
                        value = tc.getStdErrorOfMean();
                        break;
                    case 2:
                        value = tc.getMedian();
                        break;
                    case 3:
                        value = tc.getLowerHPD();
                        break;
                    case 4:
                        value = tc.getUpperHPD();
                        break;
                    case 5:
                        value = tc.getACT();
                        break;
                    case 6:
                        value = tc.getESS();
                        break;
                }
            } else if (td != null) {
                if (row != 0 && !td.isValid()) return "n/a";

                switch (row) {
                    case 0:
                        value = td.getMean();
                        break;
                    case 1:
                        return "n/a";
                    case 2:
                        value = td.getMedian();
                        break;
                    case 3:
                        value = td.getLowerHPD();
                        break;
                    case 4:
                        value = td.getUpperHPD();
                        break;
                    case 5:
                        return "n/a";
                    case 6:
                        value = td.getESS();
                        break;
                }
            } else {
                return "-";
            }

            if (Math.abs(value) < 0.1 || Math.abs(value) >= 10000.0) {
                return formatter.format(value);
            } else return formatter2.format(value);
        }

        public String getColumnName(int column) {
            if (column == 0) return "Summary Statistic";
            if (traceLists != null && traceIndices != null) {
                int n1 = (column - 1) / traceIndices.length;
                int n2 = (column - 1) % traceIndices.length;
                return traceLists[n1].getTraceName(traceIndices[n2]);
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
