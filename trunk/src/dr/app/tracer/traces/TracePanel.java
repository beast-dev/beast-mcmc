package dr.app.tracer.traces;

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

	private JTabbedPane tabbedPane = new JTabbedPane();
	private JFrame parent = null;

	private SummaryStatisticsPanel summaryPanel;
	private DensityPanel densityPanel;
	private CorrelationPanel correlationPanel;
	private RawTracePanel tracePanel;

	Icon traceIcon = null;
	Icon frequencyIcon = null;
	Icon densityIcon = null;
	Icon summaryIcon = null;
	Icon correlationIcon = null;

	/** Creates new form TracePanel */
	public TracePanel(JFrame parent) {

		this.parent = parent;

		traceIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/trace-small-icon.gif"));
		frequencyIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/frequency-small-icon.gif"));
		densityIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/density-small-icon.gif"));
		summaryIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/summary-small-icon.png"));
		correlationIcon = new ImageIcon(IconUtils.getImage(TracePanel.class, "images/correlation-small-icon.gif"));

		summaryPanel = new SummaryStatisticsPanel();
		densityPanel = new DensityPanel();
		correlationPanel = new CorrelationPanel();
		tracePanel = new RawTracePanel();

		tabbedPane.addTab("Estimates", summaryIcon, summaryPanel);
		tabbedPane.addTab("Marginal Density", densityIcon, densityPanel);
		tabbedPane.addTab("Joint-Marginal", correlationIcon, correlationPanel);
		tabbedPane.addTab("Trace", traceIcon, tracePanel);

		setLayout(new BorderLayout());
		add(tabbedPane, BorderLayout.CENTER);
	}

    /** This function takes a multiple statistics in a single log files */
    public void setTraces(TraceList[] traceLists, int[] traces) {

        summaryPanel.setTraces(traceLists, traces);
        densityPanel.setTraces(traceLists, traces);
        correlationPanel.setTraces(traceLists, traces);
        tracePanel.setTraces(traceLists, traces);
    }

	public String getExportText() {
		switch (tabbedPane.getSelectedIndex()) {
			case 0: return summaryPanel.toString();
	        case 1: return densityPanel.toString();
			case 2: return correlationPanel.toString();
			case 3: return tracePanel.toString();
		}
	    return "";
	}

    public void doCopy() {
    	summaryPanel.copyToClipboard();
    	switch (tabbedPane.getSelectedIndex()) {
			case 0: summaryPanel.copyToClipboard(); break;
			case 1: densityPanel.copyToClipboard(); break;
			case 2: correlationPanel.copyToClipboard(); break;
			case 3: tracePanel.copyToClipboard(); break;
		}
	}

    public JComponent getExportableComponent() {

		JComponent exportable = null;
		Component comp = tabbedPane.getSelectedComponent();

		if (comp instanceof Exportable) {
			exportable = ((Exportable)comp).getExportableComponent();
		} else if (comp instanceof JComponent) {
			exportable = (JComponent)comp;
		}

		return exportable;
	}

	//************************************************************************
	// private methods
	//************************************************************************

}
