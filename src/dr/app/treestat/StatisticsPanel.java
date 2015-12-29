/*
 * StatisticsPanel.java
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

package dr.app.treestat;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import dr.app.treestat.statistics.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;
import jam.table.TableRenderer;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;


public class StatisticsPanel extends OptionsPanel implements Exportable {


    /**
     *
     */
    private static final long serialVersionUID = -8026203872020056264L;
    TreeStatFrame frame;
    ArrayList<TreeSummaryStatistic.Factory> availableStatistics = new ArrayList<TreeSummaryStatistic.Factory>();
    TreeStatData treeStatData = null;

    JScrollPane scrollPane1 = null;
    JScrollPane scrollPane2 = null;
    JTable includedStatisticsTable = null;
    JTable availableStatisticsTable = null;
    AvailableStatisticsTableModel availableStatisticsTableModel = null;
    IncludedStatisticsTableModel includedStatisticsTableModel = null;
    TreeSummaryStatisticLabel statisticLabel = new TreeSummaryStatisticLabel(null);

    public StatisticsPanel(TreeStatFrame frame, TreeStatData treeStatData) {

        this.frame = frame;
        this.treeStatData = treeStatData;

        // default
        //treeStatDatastics.add(TreeSummaryStatistic.Utils.createTMRCAStatistic());

        // add generic tree statistics here
        availableStatistics.add(TreeLength.FACTORY);
        availableStatistics.add(TreeHeight.FACTORY);
        availableStatistics.add(NodeHeights.FACTORY);
        availableStatistics.add(InternalBranchLengths.FACTORY);
        availableStatistics.add(InternalBranchRates.FACTORY);
        availableStatistics.add(ExternalBranchRates.FACTORY);
        availableStatistics.add(InternalNodeAttribute.FACTORY);
        availableStatistics.add(RootToTipLengths.FACTORY);
        availableStatistics.add(TMRCASummaryStatistic.FACTORY);
        availableStatistics.add(CladeMRCAAttributeStatistic.FACTORY);
        availableStatistics.add(CladeMeanAttributeStatistic.FACTORY);
        availableStatistics.add(BetaTreeDiversityStatistic.FACTORY);
//        availableStatistics.add(MeanRootToTipLength.FACTORY);
//        availableStatistics.add(MedianRootToTipLength.FACTORY);

        availableStatistics.add(B1Statistic.FACTORY);
        availableStatistics.add(CollessIndex.FACTORY);
        availableStatistics.add(CherryStatistic.FACTORY);
        availableStatistics.add(SingleChildCountStatistic.FACTORY);
        availableStatistics.add(Nbar.FACTORY);
        availableStatistics.add(TreenessStatistic.FACTORY);
        availableStatistics.add(GammaStatistic.FACTORY);
        availableStatistics.add(DeltaStatistic.FACTORY);
//        availableStatistics.add(MonophylySummaryStatistic.FACTORY);
//        availableStatistics.add(ParsimonySummaryStatistic.FACTORY);
        availableStatistics.add(ExternalInternalRatio.FACTORY);
        availableStatistics.add(FuLiD.FACTORY);
        availableStatistics.add(RankProportionStatistic.FACTORY);
        availableStatistics.add(IntervalKStatistic.FACTORY);
        availableStatistics.add(LineageCountStatistic.FACTORY);
        availableStatistics.add(LineageProportionStatistic.FACTORY);
        availableStatistics.add(MRCAOlderThanStatistic.FACTORY);


//		if (treeStatDataeeStatDatanull) {
//			for (int i = 0; i < treeStatDataSets.size(); i++) {
//				availableStatistics.add(TreeSummaryStatistic.Utils.createTMRCAStatistic((TaxonList)treeStatDataSets.get(i)));
//			}
//		}

        setOpaque(false);

        Icon includeIcon = IconUtils.getIcon(this.getClass(), "images/include.png");
        Icon excludeIcon = IconUtils.getIcon(this.getClass(), "images/exclude.png");

        // Available statistics
        availableStatisticsTableModel = new AvailableStatisticsTableModel();
        TableSorter sorter = new TableSorter(availableStatisticsTableModel);
        availableStatisticsTable = new JTable(sorter);

        availableStatisticsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        availableStatisticsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                statisticsTableSelectionChanged(false);
            }
        });

        scrollPane1 = new JScrollPane(availableStatisticsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        JPanel buttonPanel1 = createAddRemoveButtonPanel(
                includeStatisticAction, includeIcon, null,
                excludeStatisticAction, excludeIcon, null,
                javax.swing.BoxLayout.Y_AXIS);

        // Included statistics
        includedStatisticsTableModel = new IncludedStatisticsTableModel();
        sorter = new TableSorter(includedStatisticsTableModel);
        includedStatisticsTable = new JTable(sorter);

        includedStatisticsTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        includedStatisticsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                statisticsTableSelectionChanged(true);
            }
        });

        scrollPane2 = new JScrollPane(includedStatisticsTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);


        setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        c.weightx = 0.5;
        c.weighty = 0.75;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6, 6, 6, 6);
        c.gridx = 0;
        c.gridy = 0;
        add(scrollPane1, c);

        c.weightx = 0;
        c.weighty = 0.75;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6, 0, 6, 0);
        c.gridx = 1;
        c.gridy = 0;
        add(buttonPanel1, c);

        c.weightx = 0.5;
        c.weighty = 0.75;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6, 6, 6, 6);
        c.gridx = 2;
        c.gridy = 0;
        add(scrollPane2, c);

        c.weightx = 1.0;
        c.weighty = 0.25;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.CENTER;
        c.insets = new Insets(6, 6, 6, 6);
        c.gridx = 0;
        c.gridy = 1;
        c.gridwidth = 3;
        add(statisticLabel, c);

        setBorder(BorderFactory.createEmptyBorder(6, 6, 6, 6));
    }

    JPanel createAddRemoveButtonPanel(Action addAction, Icon addIcon, String addToolTip,
                                      Action removeAction, Icon removeIcon, String removeToolTip, int axis) {

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new BoxLayout(buttonPanel, axis));
        buttonPanel.setOpaque(false);
        JButton addButton = new JButton(addAction);
        if (addIcon != null) {
            addButton.setIcon(addIcon);
            addButton.setText(null);
        }
        addButton.setToolTipText(addToolTip);
        addButton.putClientProperty("JButton.buttonType", "textured");
        addButton.setOpaque(false);
        addAction.setEnabled(false);

        JButton removeButton = new JButton(removeAction);
        if (removeIcon != null) {
            removeButton.setIcon(removeIcon);
            removeButton.setText(null);
        }
        removeButton.setToolTipText(removeToolTip);
        removeButton.putClientProperty("JButton.buttonType", "textured");
        removeButton.setOpaque(false);
        removeAction.setEnabled(false);

        buttonPanel.add(addButton);
        buttonPanel.add(new JToolBar.Separator(new Dimension(6, 6)));
        buttonPanel.add(removeButton);

        return buttonPanel;
    }

    public JComponent getExportableComponent() {
        return this;
    }

    public void dataChanged() {
        availableStatisticsTableModel.fireTableDataChanged();
        includedStatisticsTableModel.fireTableDataChanged();
    }

    private void statisticsTableSelectionChanged(boolean includedTable) {
        if (includedTable) {
            int index = includedStatisticsTable.getSelectedRow();
            if (index != -1) {
                availableStatisticsTable.clearSelection();
                SummaryStatisticDescription ssd = treeStatData.statistics.get(index);
                statisticLabel.setSummaryStatisticDescription(ssd);
                excludeStatisticAction.setEnabled(true);
            } else {
                excludeStatisticAction.setEnabled(false);
            }
        } else {
            int index = availableStatisticsTable.getSelectedRow();
            if (index != -1) {
                includedStatisticsTable.clearSelection();
                SummaryStatisticDescription ssd = availableStatistics.get(index);
                statisticLabel.setSummaryStatisticDescription(ssd);
                includeStatisticAction.setEnabled(true);
            } else {
                includeStatisticAction.setEnabled(false);
            }
        }
    }

    Action includeStatisticAction = new AbstractAction("->") {
        /**
         *
         */
        private static final long serialVersionUID = -7179224487959650620L;

        public void actionPerformed(ActionEvent ae) {
            int[] indices = availableStatisticsTable.getSelectedRows();
            for (int i = indices.length - 1; i >= 0; i--) {
                TreeSummaryStatistic.Factory ssd = availableStatistics.get(indices[i]);
                TreeSummaryStatistic tss = createStatistic(ssd);
                if (tss != null) {
                    treeStatData.statistics.add(tss);
                }
            }
            frame.fireDataChanged();
            dataChanged();
        }
    };

    Action excludeStatisticAction = new AbstractAction("<-") {
        /**
         *
         */
        private static final long serialVersionUID = -3904236403703620633L;

        public void actionPerformed(ActionEvent ae) {
            int[] indices = includedStatisticsTable.getSelectedRows();
            for (int i = indices.length - 1; i >= 0; i--) {
                treeStatData.statistics.remove(indices[i]);
            }

            frame.fireDataChanged();
            dataChanged();
        }
    };

    public TreeSummaryStatistic createStatistic(TreeSummaryStatistic.Factory factory) {

        if (!factory.allowsTaxonList() &&
                !factory.allowsDouble() && !factory.allowsInteger() && !factory.allowsString()) {
            return factory.createStatistic();
        }

        OptionsPanel optionPanel = new OptionsPanel();

        optionPanel.addSpanningComponent(new JLabel(factory.getSummaryStatisticDescription()));

        final JRadioButton wholeTreeRadio = new JRadioButton("For the whole tree", false);
        final JRadioButton taxonSetRadio = new JRadioButton("Using a given taxon set", false);
        final JComboBox taxonSetCombo = new JComboBox();
        final JTextField valueField;

        if (factory.allowsTaxonList()) {

            for (Object taxonSet : treeStatData.taxonSets) {
                taxonSetCombo.addItem(taxonSet);
            }

            ButtonGroup group = new ButtonGroup();

            ItemListener listener = new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                    taxonSetCombo.setEnabled(taxonSetRadio.isSelected());
                }
            };

            if (factory.allowsWholeTree()) {
                group.add(wholeTreeRadio);
                wholeTreeRadio.addItemListener(listener);

                optionPanel.addSpanningComponent(wholeTreeRadio);
                optionPanel.addSeparator();
            }

            if (factory.allowsTaxonList()) {
                group.add(taxonSetRadio);
                taxonSetRadio.addItemListener(listener);

                optionPanel.addSpanningComponent(taxonSetRadio);
                optionPanel.addComponentWithLabel("Taxon Set: ", taxonSetCombo);
                optionPanel.addSeparator();
            }

            if (factory.allowsTaxonList()) {
                taxonSetRadio.setSelected(true);
            }
            if (factory.allowsWholeTree()) {
                wholeTreeRadio.setSelected(true);
            }
        }

        if (factory.allowsDouble() || factory.allowsInteger() || factory.allowsString()) {
            if (factory.allowsDouble()) {
                valueField = new RealNumberField();
                valueField.setColumns(12);
                optionPanel.addComponentWithLabel(factory.getValueName(), valueField);
            } else if (factory.allowsInteger()) {
                valueField = new WholeNumberField();
                valueField.setColumns(12);
                optionPanel.addComponentWithLabel(factory.getValueName(), valueField);
            } else {  // allowsString
                valueField = new JTextField();
                valueField.setColumns(24);
                optionPanel.addComponentWithLabel(factory.getValueName(), valueField);
            }
        } else {
            valueField = null;
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        JDialog dialog = optionPane.createDialog(frame, factory.getSummaryStatisticName());
        //		dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
        dialog.pack();
        dialog.setVisible(true);

        if (optionPane.getValue() == null) {
            return null;
        }

        int result = (Integer) optionPane.getValue();
        if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
            return null;
        }

        TreeSummaryStatistic statistic = factory.createStatistic();

        if (wholeTreeRadio.isSelected()) {
            statistic = factory.createStatistic();
        } else if (taxonSetRadio.isSelected()) {
            TreeStatData.TaxonSet t = (TreeStatData.TaxonSet) taxonSetCombo.getSelectedItem();
            Taxa taxa = new Taxa();
            taxa.setId(t.name);
            //Iterator iter = t.taxa.iterator();
            for (Object aTaxa : t.taxa) {
                String id = (String) aTaxa;
                Taxon taxon = new Taxon(id);
                taxa.addTaxon(taxon);
            }
            statistic.setTaxonList(taxa);
        } else {
            return null;
        }

        if (factory.allowsDouble()) {
            assert valueField instanceof RealNumberField;

            Double value = ((RealNumberField) valueField).getValue();
            statistic.setDouble(value);
        } else if (factory.allowsInteger()) {
            assert valueField instanceof WholeNumberField;

            Integer value = ((WholeNumberField) valueField).getValue();
            statistic.setInteger(value);
        } else if (factory.allowsString()) {
            String value = valueField.getText();
            statistic.setString(value);
        }

        return statistic;
    }

    class AvailableStatisticsTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = 86401307035717809L;

        public AvailableStatisticsTableModel() {
        }

        public int getColumnCount() {
            return 2;
        }

        public int getRowCount() {

            return availableStatistics.size();
        }

        public Object getValueAt(int row, int col) {
            if (col == 0) return availableStatistics.get(row).getSummaryStatisticName();
            return availableStatistics.get(row).getCategory();
        }

        public boolean isCellEditable(int row, int col) {
            return false;
        }

        public String getColumnName(int column) {
            if (column == 0) return "Statistic Name";
            return "Category";
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }
    }

    class IncludedStatisticsTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -7280629792388705376L;

        public IncludedStatisticsTableModel() {
        }

        public int getColumnCount() {
            return 1;
        }

        public int getRowCount() {
            if (treeStatData == null || treeStatData.statistics == null) return 0;

            return treeStatData.statistics.size();
        }

        public Object getValueAt(int row, int col) {
            if (treeStatData == null || treeStatData.statistics == null) return null;
            if (col == 0) return treeStatData.statistics.get(row).getSummaryStatisticName();
            return treeStatData.statistics.get(row).getSummaryStatisticDescription();
        }

        public boolean isCellEditable(int row, int col) {
            return false;
        }

        public String getColumnName(int column) {
            if (column == 0) return "Statistic Name";
            return "Description";
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }
    }

    class TreeSummaryStatisticLabel extends JLabel {

        /**
         *
         */
        private static final long serialVersionUID = -5204925491148650874L;

        public TreeSummaryStatisticLabel(SummaryStatisticDescription statistic) {
            setSummaryStatisticDescription(statistic);
            setVerticalAlignment(JLabel.TOP);
            setHorizontalAlignment(JLabel.LEFT);
        }

        public void setSummaryStatisticDescription(SummaryStatisticDescription statistic) {
            String html = "";
            if (statistic != null) {
                html = "<html><body><h3>" + statistic.getSummaryStatisticName() + "</h3>" +
                        "<em>" + statistic.getSummaryStatisticDescription() + "</em>";

                html += "<ul>";
                if (!statistic.allowsNonultrametricTrees()) {
                    html += "<li>Trees must be ultrametric.</li>";
                } else if (!statistic.allowsUnrootedTrees()) {
                    html += "<li>Trees must be rooted.</li>";
                }
                if (!statistic.allowsPolytomies()) {
                    html += "<li>Trees must be strictly bifurcating.</li>";
                }
                html += "</ul>";
                String ref = statistic.getSummaryStatisticReference();
                if (ref != null && !ref.equals("") && !ref.equals("-")) {
                    html += "Reference: " + ref;
                }
                html += "</body></html>";
            }

            setText(html);
        }
    }

}
