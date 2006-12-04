package dr.app.tracer.analysis;

import org.virion.jam.framework.AuxilaryFrame;
import org.virion.jam.framework.DocumentFrame;
import dr.util.Variate;

import javax.swing.*;
import java.awt.*;


public class DemographicFrame extends AuxilaryFrame {

	private Variate xData;
	private Variate yDataMean;
	private Variate yDataMedian;
	private Variate yDataUpper;
	private Variate yDataLower;

	DemographicPlotPanel demographicPlotPanel = null;

	public DemographicFrame(DocumentFrame frame) {

		super(frame);

		demographicPlotPanel = new DemographicPlotPanel(this);

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
		setVisible(true);
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
