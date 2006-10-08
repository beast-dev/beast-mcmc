/*
 * DemographicFrame.java
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

import dr.util.Variate;
import org.virion.jam.framework.AuxilaryFrame;

import javax.swing.*;
import java.awt.*;


public class DemographicFrame extends AuxilaryFrame {
	
	private Variate xData;
	private Variate yDataMean;
	private Variate yDataMedian;
	private Variate yDataUpper;
	private Variate yDataLower;

	DemographicPlotPanel demographicPlotPanel = null;
	
	public DemographicFrame(TracerFrame tracerFrame) {
	
		super(tracerFrame);
		
		demographicPlotPanel = new DemographicPlotPanel();

		setContentsPanel(demographicPlotPanel);
	
		getSaveAction().setEnabled(false);
		getSaveAsAction().setEnabled(false);

		getCutAction().setEnabled(false);
		getCopyAction().setEnabled(true);
		getPasteAction().setEnabled(false);
		getDeleteAction().setEnabled(false);
		getSelectAllAction().setEnabled(false);
		getFindAction().setEnabled(false);

		getZoomWindowAction().setEnabled(false);
	}
	
	public void initializeComponents() {

		setSize(new java.awt.Dimension(640, 480));
	}

    public void setupDemographic(String title, Variate xData,
                                 Variate yDataMean, Variate yDataMedian,
                                 Variate yDataUpper, Variate yDataLower,
                                 double timeMean, double timeMedian,
                                 double timeUpper, double timeLower) {

		this.xData = xData;
		this.yDataMean = yDataMean;
		this.yDataMedian = yDataMedian;
		this.yDataUpper = yDataUpper;
		this.yDataLower = yDataLower;
		
		demographicPlotPanel.setupPlot(title, xData, yDataMean, yDataMedian, yDataUpper, yDataLower,
                timeMean, timeMedian, timeUpper, timeLower);
		show();
	}
	
	public boolean useExportAction() { return true; }

    public JComponent getExportableComponent() {
		return demographicPlotPanel.getExportableComponent();
	} 	
      	
	public void doCopy() {
		java.awt.datatransfer.Clipboard clipboard = 
			Toolkit.getDefaultToolkit().getSystemClipboard();

		java.awt.datatransfer.StringSelection selection = 
			new java.awt.datatransfer.StringSelection(this.toString());

		clipboard.setContents(selection, selection);
	}

	public String toString() {
		StringBuffer buffer = new StringBuffer();

		buffer.append("Time\tMean\tMedian\tUpper\tLower\n");
		
		for (int i = 0; i < xData.getCount(); i++) {
			buffer.append(String.valueOf(xData.get(i)));
			buffer.append("\t");
			buffer.append(String.valueOf(yDataMean.get(i)));
			buffer.append("\t");
			buffer.append(String.valueOf(yDataMedian.get(i)));
			buffer.append("\t");
			buffer.append(String.valueOf(yDataUpper.get(i)));
			buffer.append("\t");
			buffer.append(String.valueOf(yDataLower.get(i)));
			buffer.append("\n");
		}
		
		return buffer.toString();
	}
}
