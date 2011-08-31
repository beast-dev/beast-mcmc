/*
 * PriorsPanel.java
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

package dr.app.pathogen;

import dr.app.gui.util.LongTask;
import dr.evolution.tree.*;
import dr.app.gui.chart.*;
import dr.stats.DiscreteStatistics;
import dr.stats.Regression;
import dr.stats.Variate;
import dr.util.NumberFormatter;
import jam.framework.Exportable;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.PrintWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.*;

import figtree.panel.FigTreePanel;
import figtree.treeviewer.TreePaneSelector;
import figtree.treeviewer.TreeSelectionListener;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.graphs.Node;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreesPanel extends JPanel implements Exportable {

    StatisticsModel statisticsModel;
    JTable statisticsTable = null;

    private Tree tree = null;
    private Tree currentTree = null;
    private Tree bestFittingRootTree = null;

    PathogenFrame frame = null;
    JTabbedPane tabbedPane = new JTabbedPane();
    JTextArea textArea = new JTextArea();

    //    JTreeDisplay treePanel;
    FigTreePanel treePanel;
    JChartPanel rootToTipPanel;
    JChart rootToTipChart;
    ScatterPlot rootToTipPlot;
    JChartPanel residualPanel;
    JChart residualChart;
    ScatterPlot residualPlot;

    Map<Node, Integer> pointMap = new HashMap<Node, Integer>();

    private boolean bestFittingRoot;
    private TemporalRooting.RootingFunction rootingFunction;
    private TemporalRooting temporalRooting = null;

    public TreesPanel(PathogenFrame parent, Tree tree) {
        statisticsModel = new StatisticsModel();
        statisticsTable = new JTable(statisticsModel);

        statisticsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.RIGHT, new Insets(0, 4, 0, 4)));
        statisticsTable.getColumnModel().getColumn(1).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        JScrollPane scrollPane = new JScrollPane(statisticsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

//        JToolBar toolBar1 = new JToolBar();
//        toolBar1.setFloatable(false);
//        toolBar1.setOpaque(false);
//        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        Box controlPanel1 = new Box(BoxLayout.PAGE_AXIS);
        controlPanel1.setOpaque(false);

        JPanel panel3 = new JPanel(new BorderLayout(0,0));
        panel3.setOpaque(false);
        rootingCheck = new JCheckBox("Best-fitting root");
        panel3.add(rootingCheck, BorderLayout.CENTER);

        controlPanel1.add(panel3);

        final JComboBox rootingFunctionCombo = new JComboBox(TemporalRooting.RootingFunction.values());
        /*
        JPanel panel4 = new JPanel(new BorderLayout(0,0));
        panel4.setOpaque(false);
        panel4.add(new JLabel("Function: "), BorderLayout.WEST);
        panel4.add(rootingFunctionCombo, BorderLayout.CENTER);
        controlPanel1.add(panel4);
        */

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));

        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.NORTH);


        treePanel = new FigTreePanel(FigTreePanel.Style.SIMPLE);
        tabbedPane.add("Tree", treePanel);

        treePanel.getTreeViewer().setSelectionMode(TreePaneSelector.SelectionMode.TAXA);
        treePanel.getTreeViewer().addTreeSelectionListener(new TreeSelectionListener() {
            public void selectionChanged() {
                treeSelectionChanged();
            }
        });

        rootToTipChart = new JChart(new LinearAxis(), new LinearAxis(Axis.AT_ZERO, Axis.AT_MINOR_TICK));

        ChartSelector selector1 = new ChartSelector(rootToTipChart);

        rootToTipPanel = new JChartPanel(rootToTipChart, "", "time", "divergence");
        rootToTipPanel.setOpaque(false);

        tabbedPane.add("Root-to-tip", rootToTipPanel);

        residualChart = new JChart(new LinearAxis(), new LinearAxis(Axis.AT_ZERO, Axis.AT_MINOR_TICK));

        ChartSelector selector2 = new ChartSelector(residualChart);

        residualPanel = new JChartPanel(residualChart, "", "time", "residual");
        residualPanel.setOpaque(false);

        tabbedPane.add("Residuals", residualPanel);

//        textArea.setEditable(false);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(tabbedPane, BorderLayout.CENTER);
//        panel2.add(textArea, BorderLayout.SOUTH);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, panel2);
        splitPane.setDividerLocation(220);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        add(splitPane, BorderLayout.CENTER);

        rootingCheck.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setBestFittingRoot(rootingCheck.isSelected(), (TemporalRooting.RootingFunction)rootingFunctionCombo.getSelectedItem());
            }
        });

        rootingFunctionCombo.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setBestFittingRoot(rootingCheck.isSelected(), (TemporalRooting.RootingFunction)rootingFunctionCombo.getSelectedItem());
            }
        });

        setTree(tree);
    }

    private void treeSelectionChanged() {
        if (rootToTipPlot != null) {
            Set<Node> selectedTips = treePanel.getTreeViewer().getSelectedTips();
            Set<Integer> selectedPoints = new HashSet<Integer>();
            for (Node node : selectedTips) {
                selectedPoints.add(pointMap.get(node));
            }
            rootToTipPlot.setSelectedPoints(selectedPoints);
            residualPlot.setSelectedPoints(selectedPoints);
        }
    }

    private void plotSelectionChanged(final Set<Integer> selectedPoints) {
        Set<String> selectedTaxa = new HashSet<String>();
        for (Integer i : selectedPoints) {
            selectedTaxa.add(tree.getTaxon(i).toString());
        }

        treePanel.getTreeViewer().selectTaxa(selectedTaxa);
    }

    public void timeScaleChanged() {
        bestFittingRootTree = null;
        if (rootingCheck.isSelected()) {
            rootingCheck.setSelected(false);
        } else {
            setupPanel();
        }
    }

    public JComponent getExportableComponent() {
        return (JComponent)tabbedPane.getSelectedComponent();
    }

    public void setTree(Tree tree) {
        this.tree = tree;
        setupPanel();
    }

    public void setBestFittingRoot(boolean bestFittingRoot, final TemporalRooting.RootingFunction rootingFunction) {
        this.bestFittingRoot = bestFittingRoot;
        if (this.rootingFunction != rootingFunction) {
            bestFittingRootTree = null;
            this.rootingFunction = rootingFunction;
        }
        if (this.bestFittingRoot && bestFittingRootTree == null) {
            findRoot();
        }

        setupPanel();
    }

    public Tree getTree() {
        return tree;
    }

    public Tree getTreeAsViewed() {
        return currentTree;
    }

    public void writeDataFile(Writer writer) {
        PrintWriter pw = new PrintWriter(writer);
        String labels[] = temporalRooting.getTipLabels(currentTree);
        double yValues[] = temporalRooting.getRootToTipDistances(currentTree);

        if (temporalRooting.isContemporaneous()) {
            double meanY = DiscreteStatistics.mean(yValues);
            pw.println("tip\tdistance\tdeviation");
            for (int i = 0; i < yValues.length; i++) {
                pw.println(labels[i] + "\t" + "\t" + yValues[i] + "\t" + (yValues[i] - meanY));
            }
        } else {
            double xValues[] = temporalRooting.getTipDates(currentTree);
            Regression r = temporalRooting.getRootToTipRegression(currentTree);
            double[] residuals = temporalRooting.getRootToTipResiduals(currentTree, r);
            pw.println("tip\tdate\tdistance\tresidual");
            for (int i = 0; i < xValues.length; i++) {
                pw.println(labels[i] + "\t" + xValues[i] + "\t" + yValues[i] + "\t" + residuals[i]);
            }
        }
    }

    public void setupPanel() {
        StringBuilder sb = new StringBuilder();
        NumberFormatter nf = new NumberFormatter(6);

        if (tree != null) {
            temporalRooting = new TemporalRooting(tree);
            currentTree = this.tree;

            if (bestFittingRoot && bestFittingRootTree != null) {
                currentTree = bestFittingRootTree;
                sb.append("Best-fitting root");
            } else {
                sb.append("User root");
            }

            if (temporalRooting.isContemporaneous()) {
                if (tabbedPane.getSelectedIndex() == 2) {
                    tabbedPane.setSelectedIndex(1);
                }
                tabbedPane.setEnabledAt(2, false);
            } else {
                tabbedPane.setEnabledAt(2, true);
            }

            RootedTree jtree = Tree.Utils.asJeblTree(currentTree);

            if (temporalRooting.isContemporaneous()) {
                double values[] = temporalRooting.getRootToTipDistances(currentTree);

                rootToTipChart.removeAllPlots();
                NumericalDensityPlot dp = new NumericalDensityPlot(values, 20, null);
                dp.setLineColor(new Color(9,70,15));

                double yOffset = dp.getYData().getMax() / 2;
                double[] dummyValues = new double[values.length];
                for (int i = 0; i < dummyValues.length; i++) {
                    dummyValues[i] = yOffset;
                }

                rootToTipPlot = new ScatterPlot(values, dummyValues);
                rootToTipPlot.setMarkStyle(Plot.CIRCLE_MARK, 5, new BasicStroke(0.5F), new Color(44,44,44), new Color(249,202,105));
                rootToTipPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44,44,44), UIManager.getColor("List.selectionBackground"));
                rootToTipPlot.addListener(new Plot.Adaptor() {
                    public void selectionChanged(final Set<Integer> selectedPoints) {
                        plotSelectionChanged(selectedPoints);
                    }
                });

                rootToTipChart.addPlot(rootToTipPlot);
                rootToTipChart.addPlot(dp);
                rootToTipPanel.setXAxisTitle("root-to-tip divergence");
                rootToTipPanel.setYAxisTitle("proportion");

                residualChart.removeAllPlots();

                sb.append(", contemporaneous tips");
                sb.append(", mean root-tip distance: " + nf.format(DiscreteStatistics.mean(values)));
                sb.append(", coefficient of variation: " + nf.format(DiscreteStatistics.stdev(values) / DiscreteStatistics.mean(values)));
                sb.append(", stdev: " + nf.format(DiscreteStatistics.stdev(values)));
                sb.append(", variance: " + nf.format(DiscreteStatistics.variance(values)));
            } else {
                Regression r = temporalRooting.getRootToTipRegression(currentTree);
                double[] residuals = temporalRooting.getRootToTipResiduals(currentTree, r);
                pointMap.clear();
                for (int i = 0; i < currentTree.getExternalNodeCount(); i++) {
                    NodeRef tip = currentTree.getExternalNode(i);
                    Node node = jtree.getNode(Taxon.getTaxon(currentTree.getNodeTaxon(tip).getId()));
                    node.setAttribute("residual", residuals[i]);

                    pointMap.put(node, i);
                }

                rootToTipChart.removeAllPlots();
                rootToTipPlot = new ScatterPlot(r.getXData(), r.getYData());
                rootToTipPlot.addListener(new Plot.Adaptor() {
                    public void selectionChanged(final Set<Integer> selectedPoints) {
                        plotSelectionChanged(selectedPoints);
                    }
                });
                rootToTipPlot.setMarkStyle(Plot.CIRCLE_MARK, 5, new BasicStroke(0.5F), new Color(44,44,44), new Color(249,202,105));
                rootToTipPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44,44,44), UIManager.getColor("List.selectionBackground"));
                rootToTipChart.addPlot(rootToTipPlot);
                rootToTipChart.addPlot(new RegressionPlot(r));
                rootToTipChart.getXAxis().addRange(r.getXIntercept(), r.getXData().getMax());
                rootToTipPanel.setXAxisTitle("time");
                rootToTipPanel.setYAxisTitle("root-to-tip divergence");

                residualChart.removeAllPlots();
                Variate values = r.getYResidualData();
                NumericalDensityPlot dp = new NumericalDensityPlot(values, 20);
                dp.setLineColor(new Color(103,128,144));

                double yOffset = dp.getYData().getMax() / 2;
                double[] dummyValues = new double[values.getCount()];
                for (int i = 0; i < dummyValues.length; i++) {
                    dummyValues[i] = yOffset;
                }
                Variate yOffsetValues = new Variate.Double(dummyValues);
                residualPlot = new ScatterPlot(values, yOffsetValues);
                residualPlot.addListener(new Plot.Adaptor() {
                    public void selectionChanged(final Set<Integer> selectedPoints) {
                        plotSelectionChanged(selectedPoints);
                    }
                });
                residualPlot.setMarkStyle(Plot.CIRCLE_MARK, 5, new BasicStroke(0.5F), new Color(44,44,44), new Color(249,202,105));
                residualPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44,44,44), UIManager.getColor("List.selectionBackground"));

                residualChart.addPlot(residualPlot);
                residualChart.addPlot(dp);
                residualPanel.setXAxisTitle("residual");
                residualPanel.setYAxisTitle("proportion");

//                residualChart.removeAllPlots();
//                residualPlot = new ScatterPlot(r.getXData(), r.getYResidualData());
//                residualPlot.addListener(new Plot.Adaptor() {
//                    public void selectionChanged(final Set<Integer> selectedPoints) {
//                        plotSelectionChanged(selectedPoints);
//                    }
//                });
//                residualChart.addPlot(residualPlot);
//                residualPanel.setXAxisTitle("residual");
//                residualPanel.setYAxisTitle("proportion");

                sb.append(", dated tips");
                sb.append(", date range: " + nf.format(temporalRooting.getDateRange()));
                sb.append(", slope (rate): " + nf.format(r.getGradient()));
                sb.append(", x-intercept (TMRCA): " + nf.format(r.getXIntercept()));
                sb.append(", corr. coeff: " + nf.format(r.getCorrelationCoefficient()));
                sb.append(", R^2: " + nf.format(r.getRSquared()));
            }

            treePanel.setTree(jtree);
            treePanel.setColourBy("residual");

        } else {
            treePanel.setTree(null);
            rootToTipChart.removeAllPlots();
            sb.append("No trees loaded");
        }

        textArea.setText(sb.toString());

        statisticsModel.fireTableStructureChanged();
        repaint();
    }

    private javax.swing.Timer timer = null;


    private void findRoot() {

//        bestFittingRootTree = temporalRooting.findRoot(tree);
        final FindRootTask analyseTask = new FindRootTask();

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Finding best-fit root",
                "", 0, tree.getNodeCount());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

        timer = new javax.swing.Timer(10, new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                progressMonitor.setProgress(analyseTask.getCurrent());
                if (progressMonitor.isCanceled() || analyseTask.done()) {
                    progressMonitor.close();
                    analyseTask.stop();
                    timer.stop();
                }
            }
        });

        analyseTask.go();
        timer.start();

    }

    class FindRootTask extends LongTask {

        public FindRootTask() {
        }

        public int getCurrent() {
            return temporalRooting.getCurrentRootBranch();
        }

        public int getLengthOfTask() {
            return temporalRooting.getTotalRootBranches();
        }

        public String getDescription() {
            return "Calculating demographic reconstruction...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {
            bestFittingRootTree = temporalRooting.findRoot(tree, rootingFunction);
            EventQueue.invokeLater(
                    new Runnable() {
                        public void run() {
                            setupPanel();
                        }
                    });

            return null;
        }

    }


    public TemporalRooting getTemporalRooting() {
        return temporalRooting;
    }

    class StatisticsModel extends AbstractTableModel {

        String[] rowNamesDatedTips = {"Date range", "Slope (rate)", "X-Intercept (TMRCA)", "Correlation Coefficient", "R squared", "Residual Mean Squared"};
        String[] rowNamesContemporaneousTips = {"Mean root-tip", "Coefficient of variation", "Stdev", "Variance"};

        private DecimalFormat formatter = new DecimalFormat("0.####E0");
        private DecimalFormat formatter2 = new DecimalFormat("####0.####");

        public StatisticsModel() {
        }

        public int getColumnCount() {
            return 2;
        }

        public int getRowCount() {
            if (temporalRooting == null) {
                return 0;
            } else if (temporalRooting.isContemporaneous()) {
                return rowNamesContemporaneousTips.length;
            } else {
                return rowNamesDatedTips.length;
            }
        }

        public Object getValueAt(int row, int col) {

            double value = 0;
            if (temporalRooting.isContemporaneous()) {
                if (col == 0) {
                    return rowNamesContemporaneousTips[row];
                }
                double values[] = temporalRooting.getRootToTipDistances(currentTree);

                switch (row) {
                    case 0:
                        value =DiscreteStatistics.mean(values);
                        break;
                    case 1:
                        value = DiscreteStatistics.stdev(values) / DiscreteStatistics.mean(values);
                        break;
                    case 2:
                        value = DiscreteStatistics.stdev(values);
                        break;
                    case 3:
                        value = DiscreteStatistics.variance(values);
                        break;
                }
            } else {
                Regression r = temporalRooting.getRootToTipRegression(currentTree);
                if (col == 0) {
                    return rowNamesDatedTips[row];
                }
                switch (row) {
                    case 0:
                        value = temporalRooting.getDateRange();
                        break;
                    case 1:
                        value = r.getGradient();
                        break;
                    case 2:
                        value = r.getXIntercept();
                        break;
                    case 3:
                        value = r.getCorrelationCoefficient();
                        break;
                    case 4:
                        value = r.getRSquared();
                        break;
                    case 5:
                        value = r.getResidualMeanSquared();
                        break;
                }
            }

            if (value > 0 && (Math.abs(value) < 0.1 || Math.abs(value) >= 100000.0)) {
                return formatter.format(value);
            } else return formatter2.format(value);
        }

        public String getColumnName(int column) {
            if (column > 0) {
                return "";
            }
            if (temporalRooting == null) {
                return "No tree loaded";
            } else if (temporalRooting.isContemporaneous()) {
                return "Contemporaneous Tips";
            } else {
                return "Dated Tips";
            }
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

        public String toString() {
            StringBuffer buffer = new StringBuffer();

            buffer.append(getColumnName(0));
            for (int j = 1; j < getColumnCount(); j++) {
                buffer.append("\t");
                buffer.append(getColumnName(j));
            }
            buffer.append("\n");

            for (int i = 0; i < getRowCount(); i++) {
                buffer.append(getValueAt(i, 0));
                for (int j = 1; j < getColumnCount(); j++) {
                    buffer.append("\t");
                    buffer.append(getValueAt(i, j));
                }
                buffer.append("\n");
            }

            return buffer.toString();
        }
    }

    private JCheckBox rootingCheck;
}
