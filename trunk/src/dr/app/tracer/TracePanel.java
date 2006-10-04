/*
 * TracePanel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.tracer;

import org.virion.jam.framework.Exportable;
import org.virion.jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that displays information about traces
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TracePanel.java,v 1.38 2005/10/03 16:27:36 rambaut Exp $
 */
public class TracePanel extends javax.swing.JPanel implements Exportable {

	private JTabbedPane tabbedPane = new JTabbedPane(); 
	private JFrame parent = null;
	
	private SummaryStatisticsPanel summaryPanel;
	private DensityPanel densityPanel;
	private CorrelationPanel correlationPanel;
	private RawTracePanel tracePanel;
	
	/** Creates new form TracePanel */
	public TracePanel(JFrame parent) {

		this.parent = parent;
	
        Icon traceIcon = IconUtils.getIcon(TracerApp.class, "images/trace-small-icon.gif");
        Icon frequencyIcon = IconUtils.getIcon(TracerApp.class, "images/frequency-small-icon.gif");
        Icon densityIcon = IconUtils.getIcon(TracerApp.class, "images/density-small-icon.gif");
        Icon summaryIcon = IconUtils.getIcon(TracerApp.class, "images/summary-small-icon.png");
        Icon correlationIcon = IconUtils.getIcon(TracerApp.class, "images/correlation-small-icon.gif");
		
		summaryPanel = new SummaryStatisticsPanel();
		densityPanel = new DensityPanel();
		correlationPanel = new CorrelationPanel();
		tracePanel = new RawTracePanel();

		tabbedPane.addTab("Estimates", summaryIcon, summaryPanel); 
		tabbedPane.addTab("Density", densityIcon, densityPanel); 
		tabbedPane.addTab("Joint-Marginal", correlationIcon, correlationPanel);
		tabbedPane.addTab("Trace", traceIcon, tracePanel); 

		setLayout(new BorderLayout());
		add(tabbedPane, BorderLayout.CENTER);
	}
	
	/** This function takes a single statistic across a number of log files */
	public void setCombinedTraces(CombinedTraces combinedTraces, int[] traces) {

		summaryPanel.setCombinedTraces(combinedTraces, traces);
		densityPanel.setCombinedTraces(combinedTraces, traces);
		correlationPanel.setCombinedTraces();
		tracePanel.setCombinedTraces(combinedTraces, traces);
	}

	/** This function takes a multiple statistics in a single log files */
	public void setTraces(TraceList traceList, int[] traces) {
	
		summaryPanel.setTraces(traceList, traces);
		densityPanel.setTraces(traceList, traces);
		
		if (traceList != null && traces != null && traces.length == 2) {
			correlationPanel.setTraces(traceList, traces[0], traces[1]);
		} else {
			correlationPanel.setTraces(null, -1, -1);
		}
		
		tracePanel.setTraces(traceList, traces);
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
        JComponent component = null;

    	switch (tabbedPane.getSelectedIndex()) {
			case 0: component = summaryPanel; break;
			case 1: component = densityPanel; break;
			case 2: component = correlationPanel; break;
			case 3: component = tracePanel; break;
		}

        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(component.toString());

        clipboard.setContents(selection, selection);
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
