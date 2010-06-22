package dr.app.tracer.traces;

import dr.inference.trace.TraceList;
import org.virion.jam.framework.Exportable;
import org.virion.jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that displays information about traces
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TracePanel.java,v 1.2 2006/05/15 23:32:30 rambaut Exp $
 */
public class TracePanel extends javax.swing.JPanel implements Exportable {

    private final JTabbedPane tabbedPane = new JTabbedPane();

    private final SummaryStatisticsPanel summaryPanel;
    private final DensityPanel densityPanel;
    private final BetterDensityPanel newDensityPanel;
    private final CorrelationPanel correlationPanel;
    private final RawTracePanel tracePanel;

    private static final boolean USE_KDE = false;


    /**
     * Creates new form TracePanel
     */
    public TracePanel(JFrame parent) {

        //JFrame parent1 = parent;

        Icon traceIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/trace-small-icon.gif"));
        Icon frequencyIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/frequency-small-icon.gif"));
        Icon densityIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/density-small-icon.gif"));
        Icon summaryIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/summary-small-icon.png"));
        Icon correlationIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/correlation-small-icon.gif"));

        summaryPanel = new SummaryStatisticsPanel(parent);
        densityPanel = new DensityPanel(parent);
        if (USE_KDE) {
            newDensityPanel = new BetterDensityPanel(parent);
        } else {
            newDensityPanel = null;
        }
        correlationPanel = new CorrelationPanel(parent);
        tracePanel = new RawTracePanel(parent);

        tabbedPane.addTab("Estimates", summaryIcon, summaryPanel);
        tabbedPane.addTab("Marginal Prob Distribution", densityIcon, densityPanel);
        if (USE_KDE) {
            tabbedPane.addTab("Better Prob Distribution", densityIcon, newDensityPanel);
        }
        tabbedPane.addTab("Joint-Marginal", correlationIcon, correlationPanel);
        tabbedPane.addTab("Trace", traceIcon, tracePanel);

        setLayout(new BorderLayout());
        add(tabbedPane, BorderLayout.CENTER);
    }

    /**
     * This function takes a multiple statistics in a single log files
     */
    public void setTraces(TraceList[] traceLists, java.util.List<String> traces) {

        summaryPanel.setTraces(traceLists, traces);
        densityPanel.setTraces(traceLists, traces);
        if (USE_KDE) {
            newDensityPanel.setTraces(traceLists, traces);
        }
        correlationPanel.setTraces(traceLists, traces);
        tracePanel.setTraces(traceLists, traces);
    }

    public void doCopy() {

        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(getExportText());

        clipboard.setContents(selection, selection);

    }

    public String getExportText() {
        switch (tabbedPane.getSelectedIndex()) {
            case 0:
                return summaryPanel.toString();
            case 1:
                return densityPanel.toString();
            case 2:
                return correlationPanel.toString();
            case 3:
                return tracePanel.toString();
        }
        return "";
    }


    public JComponent getExportableComponent() {

        JComponent exportable = null;
        Component comp = tabbedPane.getSelectedComponent();

        if (comp instanceof Exportable) {
            exportable = ((Exportable) comp).getExportableComponent();
        } else if (comp instanceof JComponent) {
            exportable = (JComponent) comp;
        }

        return exportable;
    }
}
