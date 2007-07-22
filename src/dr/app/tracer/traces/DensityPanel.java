package dr.app.tracer.traces;

import dr.gui.chart.*;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

/**
 * A panel that displays density plots of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DensityPanel.java,v 1.3 2006/11/29 09:54:30 rambaut Exp $
 */
public class DensityPanel extends JPanel implements Exportable {
	public static int COLOUR_BY_TRACE = 0;
	public static int COLOUR_BY_FILE = 1;
	public static int COLOUR_BY_ALL = 2;

	private static final Paint[] paints = new Paint[] {
			Color.BLACK,
			new Color(64,35,225),
			new Color(229,35,60),
			new Color(255,174,34),
			new Color(86,255,34),
			new Color(35,141,148),
			new Color(146,35,142),
			new Color(255,90,34),
			new Color(239,255,34),
			Color.DARK_GRAY
	};

    private int minimumBins = 100;

    private ChartSetupDialog chartSetupDialog = null;

	private JChart traceChart = new JChart(new LinearAxis(Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK), new LinearAxis());
	private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "");

    private JComboBox binsCombo = new JComboBox(
        new Integer[] { 10, 20, 50, 100, 200, 500, 1000 });

	private JCheckBox relativeDensityCheckBox = new JCheckBox("Relative density");
	private JCheckBox solidCheckBox = new JCheckBox("Fill plot");
	private JComboBox legendCombo = new JComboBox(
		new String[] { "None", "Top-Left", "Top", "Top-Right", "Left",
						"Right", "Bottom-Left", "Bottom", "Bottom-Right" }
	);
	private JComboBox colourByCombo = new JComboBox(
		new String[] { "Trace", "Trace File", "All" }
	);
	private JLabel messageLabel = new JLabel("No data loaded");

	private int colourBy = COLOUR_BY_TRACE;

	/** Creates new FrequencyPanel */
	public DensityPanel(final JFrame frame) {

		setOpaque(false);

		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());

		JToolBar toolBar = new JToolBar();
		toolBar.setOpaque(false);
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
		toolBar.setFloatable(false);

        JButton chartSetupButton = new JButton("Axes...");
		chartSetupButton.putClientProperty(
			"Quaqua.Button.style", "placard"
		);
		chartSetupButton.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(chartSetupButton);

		binsCombo.setFont(UIManager.getFont("SmallSystemFont"));
        binsCombo.setOpaque(false);
        binsCombo.setSelectedItem(100);
        JLabel label = new JLabel("Bins:");
		label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setLabelFor(binsCombo);
        toolBar.add(label);
        toolBar.add(binsCombo);

        relativeDensityCheckBox.setOpaque(false);
		relativeDensityCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
		toolBar.add(relativeDensityCheckBox);

		// Probably don't need this as an option - takes up space and
		// solid (translucent) plots look cool...
//		solidCheckBox.setOpaque(false);
//		solidCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
//		solidCheckBox.setSelected(true);
//		toolBar.add(solidCheckBox);

		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		label = new JLabel("Legend:");
		label.setFont(UIManager.getFont("SmallSystemFont"));
		label.setLabelFor(legendCombo);
		toolBar.add(label);
		legendCombo.setFont(UIManager.getFont("SmallSystemFont"));
		legendCombo.setOpaque(false);
		toolBar.add(legendCombo);

		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		label = new JLabel("Colour by:");
		label.setFont(UIManager.getFont("SmallSystemFont"));
		label.setLabelFor(colourByCombo);
		toolBar.add(label);
		colourByCombo.setFont(UIManager.getFont("SmallSystemFont"));
		colourByCombo.setOpaque(false);
		toolBar.add(colourByCombo);

		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));

		add(messageLabel, BorderLayout.NORTH);
		add(toolBar, BorderLayout.SOUTH);
		add(chartPanel, BorderLayout.CENTER);

        chartSetupButton.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (chartSetupDialog == null) {
                            chartSetupDialog = new ChartSetupDialog(frame, true, false,
                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_ZERO, Axis.AT_MAJOR_TICK);
                        }

                        chartSetupDialog.showDialog(traceChart);
                        validate();
                        repaint();
                    }
                }
        );

        binsCombo.addItemListener(
            new java.awt.event.ItemListener() {
                public void itemStateChanged(java.awt.event.ItemEvent ev) {
                    minimumBins = (Integer)binsCombo.getSelectedItem();
                    setupTraces(traceLists, traceIndices);
                }
            }
        );

		relativeDensityCheckBox.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					for (int i = 0; i < traceChart.getPlotCount(); i++) {
						((DensityPlot)traceChart.getPlot(i)).setRelativeDensity(relativeDensityCheckBox.isSelected());
					}
					traceChart.recalibrate();
					validate();
					repaint();
				}
			}
		);

		solidCheckBox.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					for (int i = 0; i < traceChart.getPlotCount(); i++) {
						((DensityPlot)traceChart.getPlot(i)).setSolid(solidCheckBox.isSelected());
					}
					traceChart.recalibrate();
					validate();
					repaint();
				}
			}
		);

		legendCombo.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					switch (legendCombo.getSelectedIndex()) {
						case 0: break;
						case 1: traceChart.setLegendAlignment(SwingConstants.NORTH_WEST); break;
						case 2: traceChart.setLegendAlignment(SwingConstants.NORTH); break;
						case 3: traceChart.setLegendAlignment(SwingConstants.NORTH_EAST); break;
						case 4: traceChart.setLegendAlignment(SwingConstants.WEST); break;
						case 5: traceChart.setLegendAlignment(SwingConstants.EAST); break;
						case 6: traceChart.setLegendAlignment(SwingConstants.SOUTH_WEST); break;
						case 7: traceChart.setLegendAlignment(SwingConstants.SOUTH); break;
						case 8: traceChart.setLegendAlignment(SwingConstants.SOUTH_EAST); break;
					}
					traceChart.setShowLegend(legendCombo.getSelectedIndex() != 0);
					validate();
					repaint();
				}
			}
		);

		colourByCombo.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					colourBy = colourByCombo.getSelectedIndex();
					setupTraces(traceLists, traceIndices);
				}
			}
		);

	}

	private TraceList[] traceLists = null;
	private int[] traceIndices = null;

    public void setTraces(TraceList[] traceLists, int[] traceIndices) {
	    this.traceLists = traceLists;
	    this.traceIndices = traceIndices;
	    setupTraces(traceLists, traceIndices);
    }


	private void setupTraces(TraceList[] traceLists, int[] traceIndices) {
        traceChart.removeAllPlots();

        if (traceLists == null || traceIndices == null || traceIndices.length == 0) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("No traces selected");
            add(messageLabel, BorderLayout.NORTH);
            return;
        }

        remove(messageLabel);

        int i = 0;
        for (TraceList tl : traceLists) {
            int n = tl.getStateCount();

            for (int j = 0; j < traceIndices.length; j++) {
	            double values[] = new double[n];
                tl.getValues(traceIndices[j], values);
                String name = tl.getTraceName(traceIndices[j]);
                if (traceLists.length > 1) {
                    name = tl.getName() + " - " + name;
                }
                DensityPlot plot = new DensityPlot(values, minimumBins);
                plot.setName(name);
                if (tl instanceof CombinedTraces) {
                    plot.setLineStyle(new BasicStroke(2.0f), paints[i]);
                } else {
                    plot.setLineStyle(new BasicStroke(1.0f), paints[i]);
                }

                traceChart.addPlot(plot);

	            if (colourBy == COLOUR_BY_TRACE || colourBy == COLOUR_BY_ALL) {
		            i++;
	            }
	            if (i == paints.length) i = 0;
            }
	        if (colourBy == COLOUR_BY_FILE) {
		        i++;
	        } else if (colourBy == COLOUR_BY_TRACE) {
		        i = 0;
	        }
            if (i == paints.length) i = 0;
        }

        if (traceLists.length == 1) {
            chartPanel.setXAxisTitle(traceLists[0].getName());
        } else if (traceIndices.length == 1) {
            chartPanel.setXAxisTitle(traceLists[0].getTraceName(traceIndices[0]));
        } else {
            chartPanel.setXAxisTitle("Multiple Traces");
        }

		validate();
		repaint();
	}


    public JComponent getExportableComponent() {
		return chartPanel;
	}

	public void copyToClipboard() {
		java.awt.datatransfer.Clipboard clipboard =
			Toolkit.getDefaultToolkit().getSystemClipboard();

		java.awt.datatransfer.StringSelection selection =
			new java.awt.datatransfer.StringSelection(this.toString());

		clipboard.setContents(selection, selection);
	}

	public String toString() {
		return "";
	}

}
