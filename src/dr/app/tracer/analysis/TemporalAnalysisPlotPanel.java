package dr.app.tracer.analysis;

import dr.gui.chart.*;
import dr.util.Variate;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;

/**
 * A panel that displays demographic plot
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TemporalAnalysisPlotPanel.java,v 1.2 2006/12/01 09:36:17 rambaut Exp $
 */
public class TemporalAnalysisPlotPanel extends JPanel {


	private JChart demoChart = new JChart(new LinearAxis(Axis.AT_DATA, Axis.AT_DATA), new LogAxis());
	private JChartPanel chartPanel = new JChartPanel(demoChart, null, "", "");

	private JComboBox meanMedianComboBox = new JComboBox(new String[] { "Median", "Mean" });
	private JCheckBox solidIntervalCheckBox = new JCheckBox("Solid interval");

	private ChartSetupDialog chartSetupDialog = null;


	/** Creates new FrequencyPanel */
	public TemporalAnalysisPlotPanel(final JFrame frame) {
		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());

		JToolBar toolBar = new JToolBar();
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
		toolBar.setFloatable(false);

		JLabel label = new JLabel("Show:");
		label.setLabelFor(meanMedianComboBox);
		toolBar.add(label);

		toolBar.add(meanMedianComboBox);

		JButton chartSetupButton = new JButton("Setup Axes");
		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		toolBar.add(chartSetupButton);

		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		toolBar.add(solidIntervalCheckBox);

		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));

		meanMedianComboBox.setFont(UIManager.getFont("SmallSystemFont"));

		add(chartPanel, BorderLayout.CENTER);
		add(toolBar, BorderLayout.SOUTH);

		meanMedianComboBox.addItemListener(
				new java.awt.event.ItemListener() {
					public void itemStateChanged(java.awt.event.ItemEvent ev) {
						updatePlots();
					}
				}
		);

		chartSetupButton.addActionListener(
				new java.awt.event.ActionListener() {
					public void actionPerformed(ActionEvent actionEvent) {
						if (chartSetupDialog == null) {
							chartSetupDialog = new ChartSetupDialog(frame, false, true,
									Axis.AT_DATA, Axis.AT_DATA, Axis.AT_DATA, Axis.AT_DATA);
						}

						chartSetupDialog.showDialog(demoChart);
						validate();
						repaint();
					}
				}
		);

		solidIntervalCheckBox.addItemListener(
				new java.awt.event.ItemListener() {
					public void itemStateChanged(java.awt.event.ItemEvent ev) {
						updatePlots();
					}
				}
		);

		chartPanel.setTitle("");
		chartPanel.setXAxisTitle("Time");
		chartPanel.setYAxisTitle("");
	}

	public void addDemographicPlot(String title, Variate xData,
	                               Variate yDataMean, Variate yDataMedian,
	                               Variate yDataUpper, Variate yDataLower,
	                               double timeMean, double timeMedian,
	                               double timeUpper, double timeLower) {

		analysisData.add( new AnalysisData(title, xData, yDataMean, yDataMedian, yDataUpper, yDataLower,
				timeMean, timeMedian, timeUpper, timeLower));
		updatePlots();

		setVisible(true);
	}

	public void addDensityPlot(String title, Variate xData,  Variate yData) {

		analysisData.add( new AnalysisData(title, xData, yData));

		updatePlots();

		setVisible(true);
	}

	public void updatePlots() {

		demoChart.removeAllPlots();
		for (AnalysisData analysis : analysisData) {
			updatePlot(analysis);
		}

		validate();
		repaint();
	}

	public void updatePlot(AnalysisData analysis) {
		if (solidIntervalCheckBox.isSelected()) {

			AreaPlot areaPlot = new AreaPlot(analysis.xData, analysis.yDataUpper, analysis.xData, analysis.yDataLower);
			areaPlot.setLineColor(new Color(0x9999FF));
			demoChart.addPlot(areaPlot);
		} else {
			LinePlot plot = new LinePlot(analysis.xData, analysis.yDataLower);
			plot.setLineStyle(new BasicStroke(1.0F), new Color(0x9999FF));
			demoChart.addPlot(plot);

			plot = new LinePlot(analysis.xData, analysis.yDataUpper);
			plot.setLineStyle(new BasicStroke(1.0F), new Color(0x9999FF));
			demoChart.addPlot(plot);

		}

		LinePlot linePlot = null;
		if (meanMedianComboBox.getSelectedItem().equals("Median")) {
			linePlot = new LinePlot(analysis.xData, analysis.yDataMedian);
		} else {
			linePlot = new LinePlot(analysis.xData, analysis.yDataMean);
		}
		linePlot.setLineStyle(new BasicStroke(2.0F), Color.black);
		demoChart.addPlot(linePlot);

		Variate y1 = new Variate.Double();
		y1.add(demoChart.getYAxis().getMinAxis());
		y1.add(demoChart.getYAxis().getMaxAxis());

		if (analysis.timeMean > 0.0 && analysis.timeMedian > 0.0) {
			Variate x1 = new Variate.Double();
			if (meanMedianComboBox.getSelectedItem().equals("Median")) {
				x1.add(analysis.timeMedian);
				x1.add(analysis.timeMedian);

			} else {
				x1.add(analysis.timeMean);
				x1.add(analysis.timeMean);
			}
			LinePlot linePlot2 = new LinePlot(x1, y1);
			linePlot2.setLineStyle(new BasicStroke(2F, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL, 0.0F,
					new float[] {0.5F, 3.0F}, 0.0F), Color.black);
			demoChart.addPlot(linePlot2);
		}

		if (analysis.timeLower > 0.0) {
			Variate x2 = new Variate.Double();
			x2.add(analysis.timeLower);
			x2.add(analysis.timeLower);

			LinePlot linePlot3 = new LinePlot(x2, y1);
			linePlot3.setLineStyle(new BasicStroke(1.0F, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL, 0.0F,
					new float[] {0.5F, 2.0F}, 0.0F), Color.black);
			demoChart.addPlot(linePlot3);
		}

		if (analysis.timeUpper > 0.0) {
			Variate x3 = new Variate.Double();
			x3.add(analysis.timeUpper);
			x3.add(analysis.timeUpper);

			LinePlot linePlot4 = new LinePlot(x3, y1);
			linePlot4.setLineStyle(new BasicStroke(1.0F, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL, 0.0F,
					new float[] {0.5F, 2.0F}, 0.0F), Color.black);
			demoChart.addPlot(linePlot4);
		}
	}

	public JComponent getExportableComponent() {
		return chartPanel;
	}

	java.util.List<AnalysisData> analysisData = new ArrayList<AnalysisData>();
	class AnalysisData {
		public AnalysisData(String title, Variate xData, Variate yDataMean, Variate yDataMedian, Variate yDataUpper, Variate yDataLower,
		                    double timeMedian, double timeMean, double timeUpper, double timeLower) {

			this.title = title;
			this.xData = xData;
			this.yDataMean = yDataMean;
			this.yDataMedian = yDataMedian;
			this.yDataUpper = yDataUpper;
			this.yDataLower = yDataLower;
			this.timeMedian = timeMedian;
			this.timeMean = timeMean;
			this.timeUpper = timeUpper;
			this.timeLower = timeLower;
		}

		public AnalysisData(String title, Variate xData, Variate yDataMean) {
			this.title = title;
			this.xData = xData;
			this.yDataMean = yDataMean;
		}

		String title;

		Variate xData = null;
		Variate yDataMean = null;
		Variate yDataMedian = null;
		Variate yDataUpper = null;
		Variate yDataLower = null;

		double timeMedian = -1;
		double timeMean = -1;
		double timeUpper = -1;
		double timeLower = -1;
	}
}
