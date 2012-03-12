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

package dr.app.beauti.priorsPanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.hpm.HierarchicalModelComponentOptions;
import dr.app.beauti.components.hpm.HierarchicalPhylogeneticModel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockModelGroup;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.TableEditorStopper;
import dr.util.NumberFormatter;
import jam.framework.Exportable;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorsPanel extends BeautiPanel implements Exportable {
    private static final long serialVersionUID = -2936049032365493416L;
    JScrollPane scrollPane = new JScrollPane();
    JTable priorTable = null;
    PriorTableModel priorTableModel = null;
    JLabel messageLabel = new JLabel();
    JButton linkButton = null;

    public ArrayList<Parameter> parameters = new ArrayList<Parameter>();

    BeautiFrame frame = null;
    BeautiOptions options = null;

    boolean hasUndefinedPrior = false;
    boolean hasImproperPrior = false;
    boolean hasRate = false;

    private final boolean isDefaultOnly;

    private final static boolean HIERARCHICAL_ENABLED = true; // Change to false to disable

    public PriorsPanel(BeautiFrame parent, boolean isDefaultOnly) {
        this.frame = parent;
        this.isDefaultOnly = isDefaultOnly;

        priorTableModel = new PriorTableModel(this);
        priorTable = new JTable(priorTableModel);

        priorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_LAST_COLUMN);
        priorTable.getTableHeader().setReorderingAllowed(false);
//        priorTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        priorTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(0).setMinWidth(200);

        priorTable.getColumnModel().getColumn(1).setCellRenderer(
                new ButtonRenderer(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
        priorTable.getColumnModel().getColumn(1).setCellEditor(
                new ButtonEditor(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
        priorTable.getColumnModel().getColumn(1).setMinWidth(260);

        priorTable.getColumnModel().getColumn(2).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(2).setMinWidth(30);


        priorTable.getColumnModel().getColumn(3).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(3).setMinWidth(410);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(priorTable);

        scrollPane = new JScrollPane(priorTable,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

        scrollPane.setOpaque(false);

        Action setHierarchicalAction = new AbstractAction("Link parameters into a phylogenetic hierarchical model") {
            public void actionPerformed(ActionEvent actionEvent) {
                // Make list of selected parameters;
                int[] rows = priorTable.getSelectedRows();
                hierarchicalButtonPressed(rows);
            }
        };

        linkButton = new JButton(setHierarchicalAction);
        linkButton.setVisible(true);
        linkButton.setEnabled(false);
        linkButton.setToolTipText(HierarchicalPhylogeneticModel.TIP_TOOL);

        messageLabel.setText(getMessage());

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        if (!isDefaultOnly) {
            add(new JLabel("Priors for model parameters and statistics:"), BorderLayout.NORTH);
        }

        add(scrollPane, BorderLayout.CENTER);

        if (HIERARCHICAL_ENABLED && !isDefaultOnly) {
            JPanel southPanel = new JPanel();
            southPanel.setLayout(new BorderLayout(0, 0));
            JToolBar toolBar1 = new JToolBar();
            toolBar1.setFloatable(false);
            toolBar1.setOpaque(false);
            toolBar1.setLayout(new BoxLayout(toolBar1, BoxLayout.X_AXIS));

            PanelUtils.setupComponent(linkButton);
            toolBar1.add(linkButton);
            southPanel.add(toolBar1, BorderLayout.NORTH);
            southPanel.add(messageLabel, BorderLayout.SOUTH);
            add(southPanel, BorderLayout.SOUTH);
        } else {
            add(messageLabel, BorderLayout.SOUTH);
        }

        priorTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

    }

    public void selectionChanged() {
        int[] selRows = priorTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
        linkButton.setEnabled(hasSelection);
    }


    public void setOptions(BeautiOptions options) {
        this.options = options;

        parameters = options.selectParameters();

        messageLabel.setText(getMessage());

        priorTableModel.fireTableDataChanged();

        validate();
        repaint();
    }

    private String getMessage() {
        String message = "<html>";
        if (isDefaultOnly) {
            message += "<ul><li>These priors listed above are still set to the default values " +
                    "and need to be reviewed.</li>";

            hasUndefinedPrior = false;
            hasImproperPrior = false;
            hasRate = false;

            for (Parameter param : parameters) {
                if (param.priorType == PriorType.UNDEFINED) {
                    hasUndefinedPrior = true;
                }
                if (param.isPriorImproper()) {
                    hasImproperPrior = true;
                }
                if (param.getBaseName().endsWith("clock.rate") || param.getBaseName().endsWith(ClockType.UCED_MEAN)
                        || param.getBaseName().endsWith(ClockType.UCLD_MEAN)) {
                    hasRate = true;
                }
            }

            if (hasUndefinedPrior) {
                message += "<li><b><font color=\"#E42217\">These priors need to be defined by the user.</font></b></li>";
            }
            if (hasImproperPrior) {
                message += "<li><b><font color=\"#B4B417\">Warning: one or more parameters have improper priors.</font></b></li>";
            }
            if (hasRate) {
                message += "<li>" +
                        "Priors on clock rates in particular should be proper such as a high variance lognormal or gamma distribution " +
                        "with a mean appropriate for the organism and units of time employed.</li>";
            }
            message += "</ul>";

        } else {
            message += "* Marked parameters currently have a default prior distribution. " +
                    "You should check that these are appropriate.";
        }


        return message + "</html>";
    }

    public void setParametersList(BeautiOptions options) {
        this.options = options;

        parameters.clear();
        for (Parameter param : options.selectParameters()) {
            if (!param.isPriorEdited() || param.isPriorImproper()) {
                parameters.add(param);
            }
        }
        messageLabel.setText(getMessage());

        if (continueButton != null) {
            continueButton.setEnabled(!hasUndefinedPrior);
        }
    }

    private PriorDialog priorDialog = null;
    //    private DiscretePriorDialog discretePriorDialog = null;
    private HierarchicalPriorDialog hierarchicalPriorDialog = null;

    private JButton continueButton = null;

    public void setContinueButton(JButton continueButton) {
        this.continueButton = continueButton;
        if (continueButton != null) {
            continueButton.setEnabled(!hasUndefinedPrior);
        }
    }

    private void hierarchicalButtonPressed(int[] rows) {

        if (rows.length < 2) {
            JOptionPane.showMessageDialog(frame,
                    "Less than two parameters selected.",
                    "Parameter linking error",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }

        Double lowerBound = null;
        Double upperBound = null;

        List<Parameter> paramList = new ArrayList<Parameter>();
        for (int i = 0; i < rows.length; ++i) {
            Parameter parameter = parameters.get(rows[i]);
            if (parameter.isStatistic) {
                JOptionPane.showMessageDialog(frame,
                        "Statistics are not currently allowed.",
                        "HPM parameter linking error",
                        JOptionPane.WARNING_MESSAGE);
                return;
            }
            boolean sameBounds = true;
            if (lowerBound == null) {
                lowerBound = parameter.truncationLower;
            } else {
                if (lowerBound != parameter.truncationLower) {
                    sameBounds = false;
                }
            }
            if (upperBound == null) {
                upperBound = parameter.truncationUpper;
            } else {
                if (upperBound != parameter.truncationUpper) {
                    sameBounds = false;
                }
            }

            if (!sameBounds) {
                JOptionPane.showMessageDialog(frame,
                        "Only parameters that share the same bounds\n" +
                                "should be included in a HPM.",
                        "HPM parameter link error",
                        JOptionPane.WARNING_MESSAGE);
                return; // Bail out
            }

            paramList.add(parameter);
        }

        if (hierarchicalPriorDialog != null) { // Already called
            // Check to see if selected parameters are already in a HPM
            HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                    options.getComponentOptions(HierarchicalModelComponentOptions.class);
            boolean anyConflicts = false;
            for (Parameter parameter : paramList) {
                if (comp.isHierarchicalParameter(parameter)) {
                    anyConflicts = true;
                    break;
                }
            }
            if (anyConflicts) {
                int option = JOptionPane.showConfirmDialog(this,
                        "At one selected parameter already exists in a HPM.\n" +
                                "Constructing a new prior will remove these parameter\n" +
                                "from the original model. Continue?",
                        "HPM warning",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                if (option == JOptionPane.NO_OPTION) {
                    return;
                }
            }
        }

        if (hierarchicalPriorDialog == null) {
            hierarchicalPriorDialog = new HierarchicalPriorDialog(frame, options);
        }

        boolean done = false;

        while (!done) {
            int result = hierarchicalPriorDialog.showDialog(paramList);
            if (result == JOptionPane.OK_OPTION && hierarchicalPriorDialog.validateModelName()) {
                hierarchicalPriorDialog.getArguments();
                done = true;
            }
            if (result == JOptionPane.CANCEL_OPTION) {
                return;
            }
        }

        // Remove parameters from old list
        for (Parameter parameter : paramList) {
            HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                    options.getComponentOptions(HierarchicalModelComponentOptions.class);
            if (comp.isHierarchicalParameter(parameter)) {
                comp.removeParameter(this, parameter, false);
            }

        }

        // Add HPM to component manager
        hierarchicalPriorDialog.addHPM(paramList);

        for (Parameter parameter : paramList) {
            parameter.setPriorEdited(true);
        }
        priorTableModel.fireTableDataChanged();
    }

    private void priorButtonPressed(int row) {
        Parameter param = parameters.get(row);

        if (HIERARCHICAL_ENABLED && hierarchicalPriorDialog != null) {
            // Check that parameter is not already in a HPM
            HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                    options.getComponentOptions(HierarchicalModelComponentOptions.class);
            if (comp.isHierarchicalParameter(param)) {
                int option = JOptionPane.showConfirmDialog(this,
                        "Parameter already exists in a HPM. Selecting a\n" +
                                "new prior will remove the parameter. Continue?",
                        "HPM warning",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                if (option == JOptionPane.NO_OPTION) {
                    return;
                }
            }
        }

        int result;
        if (priorDialog == null) {
            priorDialog = new PriorDialog(frame);
        }

        priorDialog.setParameter(param);

        do {
            result = priorDialog.showDialog();
        } while (result == JOptionPane.OK_OPTION && priorDialog.hasInvalidInput(true));

        if (result == JOptionPane.OK_OPTION) {
            // move to individual Dialog, otherwise it will change if Cancel
            priorDialog.getArguments(param);
            // Only do this if OK button is pressed (not cancel):

            if (HIERARCHICAL_ENABLED && hierarchicalPriorDialog != null) {
                // Possibly remove parameter from its HPM
                HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                        options.getComponentOptions(HierarchicalModelComponentOptions.class);
                if (comp.isHierarchicalParameter(param)) {
                    if (comp.removeParameter(this, param, true) == JOptionPane.NO_OPTION) {
                        // Bail out
                        System.err.println("Bailing out of modification");
                        return;
                    }
                    System.err.println("Parameter removed from an HPM");
                }
            }

            param.setPriorEdited(true);

            if (isDefaultOnly) {
                setParametersList(options);
            }

            if (param.getBaseName().endsWith("treeModel.rootHeight") || param.taxaId != null) { // param.taxa != null is TMRCA

                if (options.treeModelOptions.isNodeCalibrated(param)) {
                    List<ClockModelGroup> groupList;
                    if (options.useStarBEAST) {
                        groupList = options.clockModelOptions.getClockModelGroups();
                    } else {
                        groupList = options.clockModelOptions.getClockModelGroups(options.getDataPartitions(param.getOptions()));
                    }

                    for (ClockModelGroup clockModelGroup : groupList) {
                        options.clockModelOptions.nodeCalibration(clockModelGroup);
                    }
                    frame.setAllOptions();
//        	} else {
//        		options.clockModelOptions.fixRateOfFirstClockPartition();
                }
            }

            priorTableModel.fireTableDataChanged();
        }

    }

    public void getOptions(BeautiOptions options) {
    }

    public JComponent getExportableComponent() {
        return priorTable;
    }

    NumberFormatter formatter = new NumberFormatter(4);

    class DoubleRenderer extends TableRenderer {

        private static final long serialVersionUID = -2614341608257369805L;

        public DoubleRenderer(int alignment, Insets insets) {

            super(true, alignment, insets);
        }

        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column) {

            String s;
            if (((Double) value).isNaN()) {
                s = "random";
            } else {
                s = formatter.format((Double) value);
            }
            return super.getTableCellRendererComponent(table, s, isSelected, hasFocus, row, column);

        }
    }

    public class ButtonRenderer extends JButton implements TableCellRenderer {
        protected Color undefinedColour = new Color(0xE4, 0x22, 0x17);
        protected Color improperColour = new Color(0xB4, 0xB4, 0x17);


        private static final long serialVersionUID = -2416184092883649169L;

        public ButtonRenderer(int alignment, Insets insets) {
            setOpaque(true);
            setHorizontalAlignment(alignment);
            setMargin(insets);

//            setFont(UIManager.getFont("SmallSystemFont"));
//            putClientProperty("JComponent.sizeVariant", "small");
//            putClientProperty("JButton.buttonType", "square");
        }

        public Component getTableCellRendererComponent(JTable table, Object value,
                                                       boolean isSelected, boolean hasFocus, int row, int column) {
            setEnabled(table.isEnabled());
            setFont(table.getFont());
            if (isSelected) {
                //setForeground(table.getSelectionForeground());
                setBackground(table.getSelectionBackground());
            } else {
                //setForeground(table.getForeground());
                setBackground(UIManager.getColor("Button.background"));
            }

            String text = (value == null) ? "" : value.toString();
            if (text.toString().startsWith("?")) {
                setForeground(undefinedColour);
            } else if (text.toString().startsWith("!")) {
                setForeground(improperColour);
            } else {
                setForeground(UIManager.getColor("Button.foreground"));
            }

            setText(text);

            return this;
        }
    }

    public class ButtonEditor extends DefaultCellEditor {

        private static final long serialVersionUID = 6372738480075411674L;
        protected JButton button;
        private String label;
        private boolean isPushed;
        private int row;

        public ButtonEditor(int alignment, Insets insets) {
            super(new JCheckBox());
            button = new JButton();
            button.setOpaque(true);
            button.setHorizontalAlignment(alignment);
            button.setMargin(insets);
            button.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    fireEditingStopped();
                }
            });
//            button.setFont(UIManager.getFont("SmallSystemFont"));
//            button.putClientProperty("JComponent.sizeVariant", "small");
//            button.putClientProperty("JButton.buttonType", "square");
        }

        public Component getTableCellEditorComponent(JTable table, Object value,
                                                     boolean isSelected, int row, int column) {
            button.setEnabled(table.isEnabled());
            button.setFont(table.getFont());
            if (isSelected) {
//                button.setForeground(table.getSelectionForeground());
                button.setBackground(table.getSelectionBackground());
            } else {
//                button.setForeground(table.getForeground());
                button.setBackground(UIManager.getColor("Button.background"));
            }
            label = (value == null) ? "" : value.toString();
            button.setText(label);
            isPushed = true;
            this.row = row;
            return button;
        }

        public Object getCellEditorValue() {
            if (isPushed) {
                priorButtonPressed(row);
            }
            isPushed = false;
            return label;
        }

        public boolean stopCellEditing() {
            isPushed = false;
            return super.stopCellEditing();
        }

        protected void fireEditingStopped() {
            super.fireEditingStopped();
        }
    }


}
