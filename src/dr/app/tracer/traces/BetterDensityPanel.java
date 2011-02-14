package dr.app.tracer.traces;

import dr.app.gui.chart.Axis;
import dr.app.gui.chart.KDENumericalDensityPlot;
import dr.app.gui.chart.KDESetupDialog;
import dr.app.gui.chart.Plot;
import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceList;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.HashMap;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class BetterDensityPanel extends DensityPanel {

    private JButton kdeSetupButton;
    private KDESetupDialog kdeSetupDialog = null;

    public BetterDensityPanel(final JFrame frame) {
        super(frame);
    }

    protected JToolBar setupToolBar(final JFrame frame) {
        JToolBar toolBar = super.setupToolBar(frame);

        kdeSetupButton = new JButton("KDE...");

                kdeSetupButton.putClientProperty(
                "Quaqua.Button.style", "placard"
        );
        kdeSetupButton.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(kdeSetupButton);

        kdeSetupButton.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (kdeSetupDialog == null) {
                            kdeSetupDialog = new KDESetupDialog(frame, true, false,
                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_ZERO, Axis.AT_MAJOR_TICK);
                        }

                        kdeSetupDialog.showDialog(densityChart);
                        validate();
                        repaint();
                    }
                }
        );

        return toolBar;
    }

    protected Plot setupDensityPlot(TraceList tl, int traceIndex, TraceCorrelation td) {
        Double values[] = tl.getValues(traceIndex, tl.getStateCount());

        Plot plot = new KDENumericalDensityPlot(values, minimumBins, td);

        densityChart.setXAxis(false, new HashMap<Integer, String>());// make HashMap empty
        chartPanel.setYAxisTitle("Density");

        relativeDensityCheckBox.setVisible(true);
        labelBins.setVisible(true);
        binsCombo.setVisible(true);
        return plot;
    }
}



