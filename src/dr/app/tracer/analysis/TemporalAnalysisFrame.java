package dr.app.tracer.analysis;

import org.virion.jam.framework.AuxilaryFrame;
import org.virion.jam.framework.DocumentFrame;
import dr.util.Variate;
import dr.app.tracer.application.TracerFileMenuHandler;

import javax.swing.*;
import java.awt.*;

public class TemporalAnalysisFrame extends AuxilaryFrame {
	private int binCount;
	private double minTime;
	private double maxTime;

	private boolean rangeSet;

	TemporalAnalysisPlotPanel temporalAnalysisPlotPanel = null;

	public TemporalAnalysisFrame(DocumentFrame frame, String title, int binCount) {
		this(frame, title, binCount, 0.0, 0.0);
		rangeSet = false;
	}

	public TemporalAnalysisFrame(DocumentFrame frame, String title, int binCount, double minTime, double maxTime) {

		super(frame);

		setTitle(title);

		this.binCount = binCount;
		this.minTime = minTime;
		this.maxTime = maxTime;

		rangeSet = true;

		temporalAnalysisPlotPanel = new TemporalAnalysisPlotPanel(this);

		setContentsPanel(temporalAnalysisPlotPanel);

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

    public void addDemographic(String title, Variate xData,
                                 Variate yDataMean, Variate yDataMedian,
                                 Variate yDataUpper, Variate yDataLower,
                                 double timeMean, double timeMedian,
                                 double timeUpper, double timeLower) {

	    if (!rangeSet) {
		    throw new RuntimeException("Range not set");
	    }

	    if (getTitle().length() == 0) {
		    setTitle(title);
	    }

		temporalAnalysisPlotPanel.addDemographicPlot(title, xData, yDataMean, yDataMedian, yDataUpper, yDataLower,
                timeMean, timeMedian, timeUpper, timeLower);
		setVisible(true);
	}

	public void addDensity(String title, Variate xData,  Variate yData) {

		if (!rangeSet) {
			throw new RuntimeException("Range not set");
		}

		temporalAnalysisPlotPanel.addDensityPlot(title, xData, yData);
		setVisible(true);
	}

	public boolean useExportAction() { return true; }

    public JComponent getExportableComponent() {
		return temporalAnalysisPlotPanel.getExportableComponent();
	}

	public void doCopy() {
		java.awt.datatransfer.Clipboard clipboard =
			Toolkit.getDefaultToolkit().getSystemClipboard();

		java.awt.datatransfer.StringSelection selection =
			new java.awt.datatransfer.StringSelection(this.toString());

		clipboard.setContents(selection, selection);
	}

	public int getBinCount() {
		return binCount;
	}

	public double getMinTime() {
		if (!rangeSet) {
			throw new RuntimeException("Range not set");
		}

		return minTime;
	}

	public double getMaxTime() {
		if (!rangeSet) {
			throw new RuntimeException("Range not set");
		}

		return maxTime;
	}

	public void setRange(double minTime, double maxTime) {
		if (rangeSet) {
			throw new RuntimeException("Range already set");
		}

		this.minTime = minTime;
		this.maxTime = maxTime;
		rangeSet = true;
	}

	public boolean isRangeSet() {
		return rangeSet;
	}

	public String toString() {
		StringBuffer buffer = new StringBuffer();

        java.util.List<TemporalAnalysisPlotPanel.AnalysisData> analyses = temporalAnalysisPlotPanel.getAnalysisData();
        buffer.append("Time");
        for (TemporalAnalysisPlotPanel.AnalysisData analysis : analyses) {
            if (analysis.isDemographic) {
                buffer.append("\t").append(analysis.title).append("\tMedian\tUpper\tLower");
            } else {
                buffer.append("\t").append(analysis.title);
            }
        }
        buffer.append("\n");

        Variate timeScale = temporalAnalysisPlotPanel.getTimeScale();
        for (int i = 0; i < timeScale.getCount(); i++) {
            buffer.append(String.valueOf(timeScale.get(i)));

            for (TemporalAnalysisPlotPanel.AnalysisData analysis : analyses) {
                if (analysis.isDemographic) {
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataMean.get(i)));
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataMedian.get(i)));
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataUpper.get(i)));
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataLower.get(i)));
                } else {
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataMean.get(i)));
                }
            }
			buffer.append("\n");
		}

		return buffer.toString();
	}
}
