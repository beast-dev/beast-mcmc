/*
 * BayesFactorsFrame.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.app.tracer.analysis;

import dr.inference.trace.MarginalLikelihoodAnalysis;
import jam.framework.AuxilaryFrame;
import jam.framework.DocumentFrame;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

public class BayesFactorsFrame extends AuxilaryFrame {

    private List<MarginalLikelihoodAnalysis> marginalLikelihoods = new ArrayList<MarginalLikelihoodAnalysis>();
    private JPanel contentPanel;
    private JComboBox transformCombo;

    private BayesFactorsModel bayesFactorsModel;
    private JTable bayesFactorsTable;

    public BayesFactorsFrame(DocumentFrame frame, String title, boolean hasErrors, boolean isAICM) {

        super(frame);

        setTitle(title);

        bayesFactorsModel = new BayesFactorsModel(hasErrors, isAICM);
        bayesFactorsTable = new JTable(bayesFactorsModel);

        JScrollPane scrollPane1 = new JScrollPane(bayesFactorsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

        contentPanel = new JPanel(new BorderLayout(0, 0));
        contentPanel.setOpaque(false);
        contentPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
                new java.awt.Insets(0, 0, 6, 0)));

        JToolBar toolBar = new JToolBar();
        toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
        toolBar.setFloatable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);

        JLabel label = new JLabel("<html>Models compared by marginal likelihood (S.E. estimated from bootstrap " +
                "replicates). Differences between log marginal likelihoods (specifically, log Bayes factors) are " +
                "reported.  Positive values indicate better relative model fit of the rows's model compared to the " +
                "column's model.</html>");
        if (isAICM) {
            label = new JLabel("<html>Models compared by AICM (S.E. estimated from bootstrap replicates). " +
                    "Lower AICM values indicate better model fit.  Differences between AICM estimates are reported. " +
                    "Positive values indicate better relative model fit of the row's model compared to the column's " +
                    "model.</html>");
             label.setFont(UIManager.getFont("SmallSystemFont"));
        }
        label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setBorder(new EmptyBorder(5, 5, 5, 5));

        panel1.add(label, BorderLayout.NORTH);
        panel1.add(scrollPane1, BorderLayout.CENTER);

        contentPanel.add(panel1, BorderLayout.CENTER);

        label = new JLabel("<html>Marginal likelihood estimated using the smoothed harmonic mean estimator following " +
                "Newton & Raftery (1994, <i>J. Roy. Statist. Soc. B.</i> <b>56</b>: 3-48) " +
                "with the modifications proprosed by Suchard et al. (2001, <i>MBE</i> <b>18</b>: 1001-1013)</html>");
        if (isAICM) {
            label = new JLabel("<html>AICM estimated using the method-of-moments estimator. Please cite: Baele, " +
                    "Lemey, Bedford, Rambaut, Suchard and Alekseyenko. Improving the accuracy of demographic and " +
                    "molecular clock model comparison while accommodating phylogenetic uncertainty. In prep.</html>");
             label.setFont(UIManager.getFont("SmallSystemFont"));
        }
        label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setBorder(new EmptyBorder(5, 5, 5, 5));

        contentPanel.add(label, BorderLayout.SOUTH);

        setContentsPanel(contentPanel);

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

    public void addMarginalLikelihood(MarginalLikelihoodAnalysis marginalLikelihood) {
        this.marginalLikelihoods.add(marginalLikelihood);
        bayesFactorsModel.fireTableStructureChanged();
        bayesFactorsTable.repaint();
    }

    public void initializeComponents() {

        setSize(new Dimension(640, 480));
    }

    public boolean useExportAction() {
        return true;
    }

    public JComponent getExportableComponent() {
        return contentPanel;
    }

    public void doCopy() {
        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(bayesFactorsModel.toString());

        clipboard.setContents(selection, selection);
    }


    class BayesFactorsModel extends AbstractTableModel {

        String[] columnNames = {"Trace", "ln P(data | model)", "S.E."};
        boolean isAICM = false;

        private DecimalFormat formatter2 = new DecimalFormat("####0.###");

        private int columnCount;

        public BayesFactorsModel(boolean hasErrors, boolean isAICM) {
            this.columnCount = (hasErrors ? 3 : 2);
            this.isAICM = isAICM;
            if (isAICM) {
                columnNames[1] = "AICM";
            }
        }

        public int getColumnCount() {
            return columnCount + marginalLikelihoods.size();
        }

        public int getRowCount() {
            return marginalLikelihoods.size();
        }

        public Object getValueAt(int row, int col) {

            if (col == 0) {
                return marginalLikelihoods.get(row).getTraceName();
            }

            if (col == 1) {
                return formatter2.format(marginalLikelihoods.get(row).getLogMarginalLikelihood());
            } else if (columnCount > 2 && col == 2) {
                return " +/- " + formatter2.format(marginalLikelihoods.get(row).getBootstrappedSE());
            } else {
                if (col - columnCount != row) {
                    double lnML1 = marginalLikelihoods.get(row).getLogMarginalLikelihood();
                    double lnML2 = marginalLikelihoods.get(col - columnCount).getLogMarginalLikelihood();
                    double lnRatio = lnML1 - lnML2;
                    if (isAICM) {
                        lnRatio = lnML2 - lnML1;
                    }
                    return formatter2.format(lnRatio);
                } else {
                    return "-";
                }
            }
        }

        public String getColumnName(int column) {
            if (column < columnCount) {
                return columnNames[column];
            }
            return marginalLikelihoods.get(column - columnCount).getTraceName();
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
}