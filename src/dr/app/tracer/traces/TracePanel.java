/*
 * TracePanel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.tracer.traces;

import dr.inference.trace.TraceList;
import jam.framework.Exportable;
import jam.util.IconUtils;

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
    private final JointDensityPanel jointDensityPanel;
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
        jointDensityPanel = new JointDensityPanel(parent);
        tracePanel = new RawTracePanel(parent);

        tabbedPane.addTab("Estimates", summaryIcon, summaryPanel);
        tabbedPane.addTab("Marginal Prob Distribution", densityIcon, densityPanel);
        tabbedPane.addTab("Joint-Marginal", correlationIcon, jointDensityPanel);
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
        jointDensityPanel.setTraces(traceLists, traces);
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
                return jointDensityPanel.toString();
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
