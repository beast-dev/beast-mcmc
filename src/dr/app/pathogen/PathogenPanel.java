/*
 * PathogenPanel.java
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

package dr.app.pathogen;

import dr.app.gui.chart.*;
import dr.app.gui.util.LongTask;
import dr.evolution.tree.*;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;
import dr.stats.Regression;
import dr.stats.Variate;
import dr.util.NumberFormatter;
import figtree.panel.FigTreePanel;
import figtree.treeviewer.TreePaneSelector;
import figtree.treeviewer.TreeSelectionListener;
import figtree.treeviewer.TreeViewer;
import jam.framework.Exportable;
import jam.panels.SearchPanel;
import jam.panels.SearchPanelListener;
import jam.table.TableRenderer;
import jam.toolbar.Toolbar;
import jam.toolbar.ToolbarAction;
import jam.toolbar.ToolbarButton;
import jam.util.IconUtils;
import jebl.evolution.graphs.Node;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;

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
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PathogenPanel extends JPanel implements Exportable {

    public static final String COLOUR = "Colour...";
    public static final String CLEAR_COLOURING = "Clear Colouring...";

    StatisticsModel statisticsModel;
    JTable statisticsTable = null;

    private Tree tree = null;
    private Tree currentTree = null;
    private Tree bestFittingRootTree = null;

    private final PathogenFrame frame;
    private final JTabbedPane tabbedPane = new JTabbedPane();
    private final JTextArea textArea = new JTextArea();
    private final JCheckBox showMRCACheck = new JCheckBox("Show ancestor traces");

    //    JTreeDisplay treePanel;
    private final SamplesPanel samplesPanel;
    private final FigTreePanel treePanel;

    private SearchPanel filterPanel;
    private JPopupMenu filterPopup;

    JChartPanel rootToTipPanel;
    JChart rootToTipChart;
    ScatterPlot rootToTipPlot;

    private static final boolean SHOW_NODE_DENSITY = true;
    JChartPanel nodeDensityPanel;
    JChart nodeDensityChart;
    ScatterPlot nodeDensityPlot;

    JChartPanel residualPanel;
    JChart residualChart;
    ScatterPlot residualPlot;

    ErrorBarPlot errorBarPlot;
    ParentPlot mrcaPlot;

    Map<Node, Integer> pointMap = new HashMap<Node, Integer>();

    Set<Integer> selectedPoints = new HashSet<Integer>();

    private boolean bestFittingRoot;
    private TemporalRooting.RootingFunction rootingFunction;
    private TemporalRooting temporalRooting = null;

    public PathogenPanel(PathogenFrame parent, TaxonList taxa, Tree tree) {
        frame = parent;

        samplesPanel = new SamplesPanel(parent, taxa);

        tabbedPane.addTab("Sample Dates", samplesPanel);

        statisticsModel = new StatisticsModel();
        statisticsTable = new JTable(statisticsModel);

        statisticsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.RIGHT, new Insets(0, 4, 0, 4)));
        statisticsTable.getColumnModel().getColumn(1).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        JScrollPane scrollPane = new JScrollPane(statisticsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

        Box controlPanel1 = new Box(BoxLayout.PAGE_AXIS);
        controlPanel1.setOpaque(false);

        JPanel panel3 = new JPanel(new BorderLayout(0, 0));
        panel3.setOpaque(false);
        rootingCheck = new JCheckBox("Best-fitting root");
        panel3.add(rootingCheck, BorderLayout.CENTER);

        controlPanel1.add(panel3);

        final JComboBox rootingFunctionCombo = new JComboBox(TemporalRooting.RootingFunction.values());

        JPanel panel4 = new JPanel(new BorderLayout(0,0));
        panel4.setOpaque(false);
        panel4.add(new JLabel("Function: "), BorderLayout.WEST);
        panel4.add(rootingFunctionCombo, BorderLayout.CENTER);
        controlPanel1.add(panel4);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));

        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.NORTH);

        // Set up tree panel

        Toolbar toolBar = new Toolbar();
        toolBar.setOpaque(false);
        toolBar.setBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, Color.darkGray));

        toolBar.setRollover(true);
        toolBar.setFloatable(false);

        Icon colourToolIcon = IconUtils.getIcon(this.getClass(), "images/coloursTool.png");

        final ToolbarAction colourToolbarAction = new ToolbarAction("Colour", COLOUR, colourToolIcon) {
            public void actionPerformed(ActionEvent e){
                colourSelected();
            }
        };
        ToolbarButton colourToolButton = new ToolbarButton(colourToolbarAction, true);
        colourToolButton.setFocusable(false);
        toolBar.addComponent(colourToolButton);

        toolBar.addFlexibleSpace();

        filterPopup = new JPopupMenu();

        final ButtonGroup bg = new ButtonGroup();
        boolean first = true;
        for (TreeViewer.TextSearchType searchType : TreeViewer.TextSearchType.values()) {
            final JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem(searchType.toString());
            if (first) {
                menuItem.setSelected(true);
                first = false;
            }
            filterPopup.add(menuItem);
            bg.add(menuItem);
        }
        filterPanel = new SearchPanel("Filter", filterPopup, true);
        filterPanel.setOpaque(false);
//        filterPanel.getSearchText().requestFocus();
        filterPanel.addSearchPanelListener(new SearchPanelListener() {

            /**
             * Called when the user requests a search by pressing return having
             * typed a search string into the text field. If the continuousUpdate
             * flag is true then this method is called when the user types into
             * the text field.
             *
             * @param searchString the user's search string
             */
            public void searchStarted(String searchString) {
                Enumeration e = bg.getElements();
                String value = null;
                while (e.hasMoreElements()) {
                    AbstractButton button = (AbstractButton)e.nextElement();
                    if (button.isSelected()) {
                        value = button.getText();
                    }
                }

                for (TreeViewer.TextSearchType searchType : TreeViewer.TextSearchType.values()) {
                    if (searchType.toString().equals(value)) {
                        treePanel.getTreeViewer().selectTaxa("!name", searchType, searchString, false);
                    }
                }
            }

            /**
             * Called when the user presses the cancel search button or presses
             * escape while the search is in focus.
             */
            public void searchStopped() {
//                treeViewer.clearSelectedTaxa();
            }
        });

        JPanel panel5 = new JPanel(new FlowLayout());
        panel5.setOpaque(false);
        panel5.add(filterPanel);
        toolBar.addComponent(panel5);

        treePanel = new FigTreePanel(FigTreePanel.Style.SIMPLE);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(treePanel, BorderLayout.CENTER);
        panel2.add(toolBar, BorderLayout.NORTH);

        tabbedPane.add("Tree", panel2);

        treePanel.getTreeViewer().setSelectionMode(TreePaneSelector.SelectionMode.TAXA);
        treePanel.getTreeViewer().addTreeSelectionListener(new TreeSelectionListener() {
            public void selectionChanged() {
                treeSelectionChanged();
            }
        });

        rootToTipChart = new JChart(new LinearAxis(), new LinearAxis(Axis.AT_ZERO, Axis.AT_MINOR_TICK));

        ChartSelector selector1 = new ChartSelector(rootToTipChart);

        rootToTipPanel = new JChartPanel(rootToTipChart, "", "time", "divergence");
        JPanel panel = new JPanel(new BorderLayout());
        panel.add(rootToTipPanel, BorderLayout.CENTER);
        panel.add(showMRCACheck, BorderLayout.SOUTH);
        panel.setOpaque(false);

        tabbedPane.add("Root-to-tip", panel);

        residualChart = new JChart(new LinearAxis(), new LinearAxis(Axis.AT_ZERO, Axis.AT_MINOR_TICK));

        ChartSelector selector2 = new ChartSelector(residualChart);

        residualPanel = new JChartPanel(residualChart, "", "time", "residual");
        residualPanel.setOpaque(false);

        tabbedPane.add("Residuals", residualPanel);

//        textArea.setEditable(false);

        JPanel panel6 = new JPanel(new BorderLayout(0, 0));
        panel6.setOpaque(false);
        panel6.add(tabbedPane, BorderLayout.CENTER);
//        panel6.add(textArea, BorderLayout.SOUTH);

        if (SHOW_NODE_DENSITY) {
            nodeDensityChart = new JChart(new LinearAxis(), new LinearAxis(Axis.AT_ZERO, Axis.AT_MINOR_TICK));
            nodeDensityPanel = new JChartPanel(nodeDensityChart, "", "time", "node density");
            JPanel panel7 = new JPanel(new BorderLayout());
            panel7.add(nodeDensityPanel, BorderLayout.CENTER);
            panel7.setOpaque(false);

            ChartSelector selector3 = new ChartSelector(nodeDensityChart);

            tabbedPane.add("Node density", panel7);
        }


        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, tabbedPane);
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
                setBestFittingRoot(rootingCheck.isSelected(), (TemporalRooting.RootingFunction) rootingFunctionCombo.getSelectedItem());
            }
        });

        rootingFunctionCombo.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setBestFittingRoot(rootingCheck.isSelected(), (TemporalRooting.RootingFunction) rootingFunctionCombo.getSelectedItem());
            }
        });

        showMRCACheck.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setupPanel();
            }
        });


        setTree(tree);
    }

    public List<String> getSelectedTips() {
        List<String> tips = new ArrayList<String>();
        jebl.evolution.trees.Tree tree = treePanel.getTreeViewer().getTrees().get(0);

        for (Node node : treePanel.getTreeViewer().getSelectedTips()) {
            tips.add(tree.getTaxon(node).getName());
        }
        return tips;
    }

    private static Color lastColor = Color.GRAY;

    private void colourSelected() {
        Color color = JColorChooser.showDialog(this, "Select Colour", lastColor);
        if (color != null) {
            treePanel.getTreeViewer().annotateSelectedTips("!color", color);
            lastColor = color;
        }
        setupPanel();
    }

    private void treeSelectionChanged() {
        Set<Node> selectedTips = treePanel.getTreeViewer().getSelectedTips();
        frame.getCopyAction().setEnabled(selectedTips != null && selectedTips.size() > 0);
        selectedPoints = new HashSet<Integer>();
        for (Node node : selectedTips) {
            selectedPoints.add(pointMap.get(node));
        }
        if (rootToTipPlot != null) {
            rootToTipPlot.setSelectedPoints(selectedPoints);
        }
        if (residualPlot != null) {
            residualPlot.setSelectedPoints(selectedPoints);
        }
        if (SHOW_NODE_DENSITY && nodeDensityPlot != null) {
            nodeDensityPlot.setSelectedPoints(selectedPoints);
        }

        selectMRCA();
    }

    private void plotSelectionChanged(final Set<Integer> selectedPoints) {
        this.selectedPoints = selectedPoints;
        Set<String> selectedTaxa = new HashSet<String>();
        for (Integer i : selectedPoints) {
            selectedTaxa.add(tree.getTaxon(i).toString());
        }

        treePanel.getTreeViewer().selectTaxa(selectedTaxa);

        selectMRCA();
    }

    private void selectMRCA() {
        if (mrcaPlot == null) return;

        if (selectedPoints != null && selectedPoints.size() > 0) {

            Set<String> selectedTaxa = new HashSet<String>();
            for (Integer i : selectedPoints) {
                selectedTaxa.add(tree.getTaxon(i).toString());
            }

            Regression r = temporalRooting.getRootToTipRegression(currentTree);
            NodeRef mrca = TreeUtils.getCommonAncestorNode(currentTree, selectedTaxa);
            double mrcaDistance1 = temporalRooting.getRootToTipDistance(currentTree, mrca);
            double mrcaTime1 = r.getX(mrcaDistance1);
            if (tree.isExternal(mrca)) {
                mrca = tree.getParent(mrca);
            }
            double mrcaDistance = temporalRooting.getRootToTipDistance(currentTree, mrca);
            double mrcaTime = r.getX(mrcaDistance);

            mrcaPlot.setSelectedPoints(selectedPoints, mrcaTime, mrcaDistance);
        } else {
            mrcaPlot.clearSelection();
        }
        repaint();
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
        return (JComponent) tabbedPane.getSelectedComponent();
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

            RootedTree jtree = dr.evolution.tree.TreeUtils.asJeblTree(currentTree);

            List<Color> colours = new ArrayList<Color>();
            for (Node tip : jtree.getExternalNodes()) {
                Taxon taxon = jtree.getTaxon(tip);
                colours.add((Color)taxon.getAttribute("!color"));
            }

            if (temporalRooting.isContemporaneous()) {
                double[] dv = temporalRooting.getRootToTipDistances(currentTree);

                List<Double> values = new ArrayList<Double>();
                for (double d : dv) {
                    values.add(d);
                }

                rootToTipChart.removeAllPlots();
                NumericalDensityPlot dp = new NumericalDensityPlot(values, 20);
                dp.setLineColor(new Color(9, 70, 15));

                double yOffset = (Double) dp.getYData().getMax() / 2;
                List<Double> dummyValues = new ArrayList<Double>();
                for (int i = 0; i < values.size(); i++) {
                    // add a random y offset to give some visual spread
                    double y = MathUtils.nextGaussian() * ((Double) dp.getYData().getMax() * 0.05);
                    dummyValues.add(yOffset + y);
                }

                rootToTipPlot = new ScatterPlot(values, dummyValues);
                rootToTipPlot.setColours(colours);
                rootToTipPlot.setMarkStyle(Plot.CIRCLE_MARK, 8, new BasicStroke(0.0F), new Color(44, 44, 44), new Color(129, 149, 149));
                rootToTipPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44, 44, 44), UIManager.getColor("List.selectionBackground"));
                rootToTipPlot.addListener(new Plot.Adaptor() {
                    @Override
                    public void markClicked(int index, double x, double y, boolean isShiftDown) {
                        rootToTipPlot.selectPoint(index, isShiftDown);
                    }

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
                sb.append(", mean root-tip distance: " + nf.format(DiscreteStatistics.mean(dv)));
                sb.append(", coefficient of variation: " + nf.format(DiscreteStatistics.stdev(dv) / DiscreteStatistics.mean(dv)));
                sb.append(", stdev: " + nf.format(DiscreteStatistics.stdev(dv)));
                sb.append(", variance: " + nf.format(DiscreteStatistics.variance(dv)));

                showMRCACheck.setVisible(false);
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

                if (showMRCACheck.isSelected()) {
                    double[] dv = temporalRooting.getParentRootToTipDistances(currentTree);

                    List<Double> parentDistances = new ArrayList<Double>();
                    for (int i = 0; i < dv.length; i++) {
                        parentDistances.add(i, dv[i]);
                    }

                    List<Double> parentTimes = new ArrayList<Double>();
                    for (int i = 0; i < parentDistances.size(); i++) {
                        parentTimes.add(i, r.getX(parentDistances.get(i)));
                    }
                    mrcaPlot = new ParentPlot(r.getXData(), r.getYData(), parentTimes, parentDistances);
                    mrcaPlot.setLineColor(new Color(105, 202, 105));
                    mrcaPlot.setLineStroke(new BasicStroke(0.5F));

                    rootToTipChart.addPlot(mrcaPlot);
                }

                if (true) {
                    double[] datePrecisions = temporalRooting.getTipDatePrecisions(currentTree);

                    Variate.D ed = new Variate.D();

                    for (int i = 0; i < datePrecisions.length; i++) {
                        ed.add(datePrecisions[i]);
                    }

                    errorBarPlot = new ErrorBarPlot(ErrorBarPlot.Orientation.HORIZONTAL, r.getXData(), r.getYData(), ed);
                    errorBarPlot.setLineColor(new Color(44, 44, 44));
                    errorBarPlot.setLineStroke(new BasicStroke(1.0F));

                    rootToTipChart.addPlot(errorBarPlot);
                }

                rootToTipPlot = new ScatterPlot(r.getXData(), r.getYData());
                rootToTipPlot.addListener(new Plot.Adaptor() {
                    public void selectionChanged(final Set<Integer> selectedPoints) {
                        plotSelectionChanged(selectedPoints);
                    }
                });

                rootToTipPlot.setColours(colours);

                rootToTipPlot.setMarkStyle(Plot.CIRCLE_MARK, 8, new BasicStroke(0.0F), new Color(44, 44, 44), new Color(129, 149, 149));
                rootToTipPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44, 44, 44), UIManager.getColor("List.selectionBackground"));

                rootToTipChart.addPlot(rootToTipPlot);

                rootToTipChart.addPlot(new RegressionPlot(r));

                rootToTipChart.getXAxis().addRange(r.getXIntercept(), (Double) r.getXData().getMax());
                rootToTipPanel.setXAxisTitle("time");
                rootToTipPanel.setYAxisTitle("root-to-tip divergence");

                residualChart.removeAllPlots();
                Variate.D values = (Variate.D) r.getYResidualData();
                NumericalDensityPlot dp = new NumericalDensityPlot(values, 20);
                dp.setLineColor(new Color(103, 128, 144));

                double yOffset = (Double) dp.getYData().getMax() / 2;
                Double[] dummyValues = new Double[values.getCount()];
                for (int i = 0; i < dummyValues.length; i++) {
                    // add a random y offset to give some visual spread
                    double y = MathUtils.nextGaussian() * ((Double) dp.getYData().getMax() * 0.05);
                    dummyValues[i] = yOffset + y;
                }
                Variate.D yOffsetValues = new Variate.D(dummyValues);
                residualPlot = new ScatterPlot(values, yOffsetValues);
                residualPlot.addListener(new Plot.Adaptor() {
                    @Override
                    public void markClicked(int index, double x, double y, boolean isShiftDown) {
                        rootToTipPlot.selectPoint(index, isShiftDown);
                    }

                    @Override
                    public void selectionChanged(final Set<Integer> selectedPoints) {
                        plotSelectionChanged(selectedPoints);
                    }
                });
                residualPlot.setColours(colours);
                residualPlot.setMarkStyle(Plot.CIRCLE_MARK, 8, new BasicStroke(0.0F), new Color(44, 44, 44), new Color(129, 149, 149));
                residualPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44, 44, 44), UIManager.getColor("List.selectionBackground"));

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

                if (SHOW_NODE_DENSITY) {
                    Regression r2 = temporalRooting.getNodeDensityRegression(currentTree);
                    nodeDensityChart.removeAllPlots();
                    nodeDensityPlot = new ScatterPlot(r2.getXData(), r2.getYData());
                    nodeDensityPlot.addListener(new Plot.Adaptor() {
                        public void selectionChanged(final Set<Integer> selectedPoints) {
                            plotSelectionChanged(selectedPoints);
                        }
                    });
                    nodeDensityPlot.setColours(colours);
                    nodeDensityPlot.setMarkStyle(Plot.CIRCLE_MARK, 8, new BasicStroke(0.0F), new Color(44, 44, 44), new Color(129, 149, 149));
                    nodeDensityPlot.setHilightedMarkStyle(new BasicStroke(0.5F), new Color(44, 44, 44), UIManager.getColor("List.selectionBackground"));

                    nodeDensityChart.addPlot(nodeDensityPlot);

                    nodeDensityChart.addPlot(new RegressionPlot(r2));

                    nodeDensityChart.getXAxis().addRange(r2.getXIntercept(), (Double) r2.getXData().getMax());
                    nodeDensityPanel.setXAxisTitle("time");
                    nodeDensityPanel.setYAxisTitle("node density");
                }

                sb.append(", dated tips");
                sb.append(", date range: " + nf.format(temporalRooting.getDateRange()));
                sb.append(", slope (rate): " + nf.format(r.getGradient()));
                sb.append(", x-intercept (TMRCA): " + nf.format(r.getXIntercept()));
                sb.append(", corr. coeff: " + nf.format(r.getCorrelationCoefficient()));
                sb.append(", R^2: " + nf.format(r.getRSquared()));

                showMRCACheck.setVisible(true);
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
                        value = DiscreteStatistics.mean(values);
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
