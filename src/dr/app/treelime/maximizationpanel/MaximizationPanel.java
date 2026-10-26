/*
 * MaximizationPanel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.treelime.maximizationpanel;

import dr.app.gui.components.WholeNumberField;
import dr.app.treelime.TreeLimeFrame;
import dr.app.treelime.TreeLimePanel;
import dr.app.treelime.options.GradientRegistry;
import dr.app.treelime.options.GradientRegistry.LikelihoodTerm;
import dr.app.treelime.options.MaximizationOptions;
import dr.app.treelime.options.MaximizationOptions.GradientMode;
import dr.app.treelime.options.MaximizationOptions.Objective;
import dr.app.treelime.options.MaximizationOptions.StepSize;
import dr.app.treelime.options.MaximizationOptions.Treatment;
import dr.app.treelime.options.MaximizationWarnings;
import dr.app.treelime.options.Parameter;
import dr.app.treelime.options.PartitionOptions;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.util.OSType;
import dr.evomodelxml.continuous.hmc.MaximizeWrtParameterParser;
import dr.xml.XMLParseException;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Replaces BEAUti's Operators and MCMC tabs: what to maximize, the L-BFGS budget, and which parameters to optimize.
 * Parameters are held fixed by unticking them here, or by fixing them in the Priors tab.
 */
public class MaximizationPanel extends TreeLimePanel {

    private static final Insets FIELD_INSETS = new Insets(3, 6, 3, 6);

    private final JComboBox<Objective> objectiveCombo = new JComboBox<Objective>(Objective.values());
    private final WholeNumberField iterationsField = new WholeNumberField(1, Integer.MAX_VALUE);
    private final JTextField maxRunTimeField = new JTextField(10);
    private final WholeNumberField maxEvaluationsField = new WholeNumberField(0, Integer.MAX_VALUE);
    private final JComboBox<GradientMode> gradientModeCombo = new JComboBox<GradientMode>(GradientMode.values());
    private final JComboBox<StepSize> stepSizeCombo = new JComboBox<StepSize>(StepSize.values());
    private final JCheckBox printToScreenCheck = new JCheckBox("Print each iteration to the screen");
    private final JTextField reportFileField = new JTextField(32);
    private final JLabel statusLabel = new JLabel();
    private final WholeNumberField replicatesField = new WholeNumberField(1, 100000);
    private final WholeNumberField seedField = new WholeNumberField(0, Integer.MAX_VALUE);

    private final ParameterTableModel tableModel = new ParameterTableModel();
    private final OptionsPanel optionsPanel;

    private final TreeLimeFrame frame;
    private TreeLimeOptions options;
    private boolean settingOptions = false;

    public MaximizationPanel(TreeLimeFrame parent) {
        this.frame = parent;
        setLayout(new BorderLayout());
        setOpaque(false);

        // Mac OS X components have more spacing round them already
        optionsPanel = new OptionsPanel(12, (OSType.isMac() ? 6 : 24));
        optionsPanel.setOpaque(false);

        objectiveCombo.setToolTipText("<html>MAP maximizes the joint density, likelihood times priors.<br>" +
                "ML maximizes the tree likelihoods only; priors then only fix parameters and give initial values.</html>");
        optionsPanel.addComponentWithLabel("Objective:", objectiveCombo);
        objectiveCombo.addActionListener(e -> {
            if (!settingOptions) {
                maximization().objective = (Objective) objectiveCombo.getSelectedItem();
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        iterationsField.setColumns(10);
        iterationsField.setBorder(new EmptyBorder(FIELD_INSETS));
        iterationsField.setToolTipText("The L-BFGS iteration limit (at least 1).");
        optionsPanel.addComponentWithLabel("Maximum iterations:", iterationsField);
        onEdit(iterationsField, () -> {
            Integer value = iterationsField.getValue();
            if (value != null && value >= 1) {
                maximization().nIterations = value;
            }
        });

        maxRunTimeField.setBorder(new EmptyBorder(FIELD_INSETS));
        maxRunTimeField.setToolTipText("<html>Stop after this much time, for example 90, 30min, 2h or 1d " +
                "(seconds by default).<br>Leave empty for no limit.</html>");
        optionsPanel.addComponentWithLabel("Maximum run time:", maxRunTimeField);
        onEdit(maxRunTimeField, () -> {
            maximization().maxRunTime = maxRunTimeField.getText().trim();
            updateStatus();
        });

        maxEvaluationsField.setColumns(10);
        maxEvaluationsField.setBorder(new EmptyBorder(FIELD_INSETS));
        maxEvaluationsField.setToolTipText("Stop after this many evaluations of the objective (0 for no limit).");
        optionsPanel.addComponentWithLabel("Maximum evaluations:", maxEvaluationsField);
        onEdit(maxEvaluationsField, () -> {
            Integer value = maxEvaluationsField.getValue();
            if (value != null && value >= 0) {
                maximization().maxEvaluations = value;
            }
        });

        optionsPanel.addSeparator();

        gradientModeCombo.setToolTipText("<html>Analytic: the tree-likelihood gradient is computed exactly where BEAST " +
                "can (substitution, site and branch-rate parameters),<br>and by finite differences otherwise; the " +
                "prior side is always taken by finite differences of the priors, which is cheap.</html>");
        optionsPanel.addComponentWithLabel("Gradient:", gradientModeCombo);
        gradientModeCombo.addActionListener(e -> {
            if (!settingOptions) {
                maximization().setGradientMode((GradientMode) gradientModeCombo.getSelectedItem());
                stepSizeCombo.setEnabled(!maximization().isAnalytic());
                tableModel.fireTableDataChanged();
                updateStatus();
                frame.setDirty();
            }
        });

        stepSizeCombo.setToolTipText("The finite-difference step of the all-numeric gradient, relative to each value.");
        optionsPanel.addComponentWithLabel("Numeric gradient step:", stepSizeCombo);
        stepSizeCombo.addActionListener(e -> {
            if (!settingOptions) {
                maximization().stepSize = (StepSize) stepSizeCombo.getSelectedItem();
                frame.setDirty();
            }
        });

        optionsPanel.addComponent(printToScreenCheck);
        printToScreenCheck.addActionListener(e -> {
            if (!settingOptions) {
                maximization().printToScreen = printToScreenCheck.isSelected();
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        reportFileField.setBorder(new EmptyBorder(FIELD_INSETS));
        reportFileField.setToolTipText("The maximizer's report: the optimum, the gradient there and why it stopped.");
        optionsPanel.addComponentWithLabel("Report file:", reportFileField);
        onEdit(reportFileField, () -> {
            String name = reportFileField.getText().trim();
            maximization().reportFileName = name.isEmpty() ? null : name;
        });

        optionsPanel.addSeparator();

        // bootstrap replicate files (plan 8 S3g): the maximization of each replicate runs in BEAST, like the original
        replicatesField.setColumns(6);
        replicatesField.setBorder(new EmptyBorder(FIELD_INSETS));
        replicatesField.setToolTipText("<html>Nonparametric bootstrap: each replicate resamples the columns of every " +
                "partition with replacement<br>and keeps the original starting tree (a fixed-topology bootstrap).</html>");
        seedField.setColumns(8);
        seedField.setBorder(new EmptyBorder(FIELD_INSETS));
        seedField.setToolTipText("The same seed always writes the same replicates.");
        JButton bootstrapButton = new JButton("Generate replicate files...");
        bootstrapButton.setToolTipText("<html>Writes one TreeLime XML per replicate, the counts that rebuild them " +
                "and a script that runs them in BEAST.</html>");
        bootstrapButton.addActionListener(e -> {
            Integer replicates = replicatesField.getValue();
            Integer seed = seedField.getValue();
            if (replicates != null && seed != null) {
                frame.doGenerateBootstrap(replicates, seed);
            }
        });
        JPanel bootstrapPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        bootstrapPanel.setOpaque(false);
        bootstrapPanel.add(replicatesField);
        bootstrapPanel.add(new JLabel("seed:"));
        bootstrapPanel.add(seedField);
        bootstrapPanel.add(bootstrapButton);
        optionsPanel.addComponentWithLabel("Bootstrap replicates:", bootstrapPanel);
        onEdit(replicatesField, () -> {
            Integer value = replicatesField.getValue();
            if (value != null && value >= 1) {
                maximization().bootstrapReplicates = value;
            }
        });
        onEdit(seedField, () -> {
            Integer value = seedField.getValue();
            if (value != null) {
                maximization().bootstrapSeed = value;
            }
        });

        JTable table = new JTable(tableModel) {
            // the reasons a parameter is held fixed and the gradient descriptions can be longer than their columns
            public String getToolTipText(java.awt.event.MouseEvent event) {
                int row = rowAtPoint(event.getPoint());
                int column = columnAtPoint(event.getPoint());
                if (row < 0 || column <= 0) {
                    return null;
                }
                Object value = getValueAt(row, column);
                return value == null || value.toString().isEmpty() ? null : value.toString();
            }
        };
        int[] widths = {60, 150, 200, 260, 90, 330};
        for (int column = 0; column < widths.length; column++) {
            table.getColumnModel().getColumn(column).setPreferredWidth(widths[column]);
        }
        table.getColumnModel().getColumn(0).setMaxWidth(70);
        table.getColumnModel().getColumn(4).setMinWidth(90);
        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setOpaque(false);

        JPanel tablePanel = new JPanel(new BorderLayout(0, 6));
        tablePanel.setOpaque(false);
        tablePanel.setBorder(new EmptyBorder(new Insets(6, 12, 12, 12)));
        tablePanel.add(new JLabel("<html>Parameters to optimize; fix one here or in the Priors tab. A tree's root " +
                "height row stands for all its internal node heights. The topology stays as in the starting tree. " +
                "A parameter that goes below its lower bound is held there (blank for none).</html>"),
                BorderLayout.NORTH);
        tablePanel.add(scrollPane, BorderLayout.CENTER);
        tablePanel.add(statusLabel, BorderLayout.SOUTH);

        add(optionsPanel, BorderLayout.NORTH);
        add(tablePanel, BorderLayout.CENTER);
    }

    private MaximizationOptions maximization() {
        return options.maximizationOptions;
    }

    // the generator also falls back to numeric when some partition is not in a tree data likelihood
    private LikelihoodTerm likelihoodGradient(Parameter parameter) {
        return maximization().isAnalytic() ? GradientRegistry.likelihoodTerm(options, parameter) : LikelihoodTerm.NUMERIC;
    }

    private void onEdit(JTextField field, Runnable update) {
        field.getDocument().addDocumentListener(new DocumentListener() {
            public void insertUpdate(DocumentEvent e) {
                changed();
            }

            public void removeUpdate(DocumentEvent e) {
                changed();
            }

            public void changedUpdate(DocumentEvent e) {
                changed();
            }

            private void changed() {
                if (!settingOptions && options != null) {
                    update.run();
                    frame.setDirty();
                }
            }
        });
    }

    @Override
    public void setOptions(TreeLimeOptions options) {
        this.options = options;
        MaximizationOptions maximization = options.maximizationOptions;

        settingOptions = true;
        objectiveCombo.setSelectedItem(maximization.objective);
        iterationsField.setValue(maximization.nIterations);
        maxRunTimeField.setText(maximization.maxRunTime);
        maxEvaluationsField.setValue((int) maximization.maxEvaluations);
        gradientModeCombo.setSelectedItem(maximization.getGradientMode());
        stepSizeCombo.setSelectedItem(maximization.stepSize);
        stepSizeCombo.setEnabled(!maximization.isAnalytic());
        printToScreenCheck.setSelected(maximization.printToScreen);
        reportFileField.setText(maximization.getReportFileName(options.fileNameStem));
        // files saved before these settings existed have zeros
        replicatesField.setValue(maximization.bootstrapReplicates >= 1 ? maximization.bootstrapReplicates : 100);
        seedField.setValue((int) Math.max(0, Math.min(Integer.MAX_VALUE, maximization.bootstrapSeed)));
        settingOptions = false;

        tableModel.setParameters(options.selectParameters());
        updateStatus();
    }

    @Override
    public void getOptions(TreeLimeOptions options) {
        // the listeners write every change as it is made
    }

    private void updateStatus() {
        if (options == null) {
            return;
        }
        MaximizationOptions maximization = maximization();
        int optimized = 0;
        int numeric = 0;
        int trees = 0;
        for (Parameter parameter : tableModel.parameters) {
            Treatment treatment = maximization.classify(parameter, false);
            if (treatment == Treatment.NODE_HEIGHTS) {
                trees++;
            } else if (treatment.isOptimized()) {
                optimized++;
                if (likelihoodGradient(parameter) == LikelihoodTerm.NUMERIC) {
                    numeric++;
                }
            }
        }
        String status = !maximization.isAnalytic() ?
                count(optimized, "parameter") + " optimized, all with numeric gradients (2 evaluations per " +
                        "dimension and step); node heights fixed." :
                count(optimized, "parameter") + (trees > 0 ? " and the node heights of " + count(trees, "tree") : "") +
                        " optimized: " + (optimized - numeric) + " parameters with analytic or no tree-likelihood " +
                        "gradients, " + numeric + " by finite differences of the tree likelihoods" +
                        (trees > 0 ? "; node heights analytic." : ".");
        StringBuilder problems = new StringBuilder();
        if (!maximization.maxRunTime.isEmpty()) {
            try {
                MaximizeWrtParameterParser.parseRunTime(maximization.maxRunTime);
            } catch (XMLParseException e) {
                problems.append("<font color=\"red\">").append(e.getMessage()).append("</font><br>");
            }
        }
        for (String warning : MaximizationWarnings.get(options)) {
            problems.append("<font color=\"#b35c00\">Warning: ").append(warning).append("</font><br>");
        }
        statusLabel.setText(problems.length() == 0 ? status : "<html>" + problems + status + "</html>");
    }

    private static String count(int n, String noun) {
        return n + " " + noun + (n == 1 ? "" : "s");
    }

    @Override
    public JComponent getExportableComponent() {
        return this;
    }

    private class ParameterTableModel extends AbstractTableModel {

        private final String[] columnNames = {"Optimize", "Parameter", "Description", "Transform", "Lower bound",
                "Likelihood gradient"};
        private List<Parameter> parameters = new ArrayList<Parameter>();

        void setParameters(List<Parameter> parameters) {
            this.parameters = new ArrayList<Parameter>(parameters);
            fireTableDataChanged();
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return parameters.size();
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class<?> getColumnClass(int column) {
            return column == 0 ? Boolean.class : String.class;
        }

        public Object getValueAt(int row, int column) {
            Parameter parameter = parameters.get(row);
            Treatment treatment = maximization().classify(parameter, false);
            switch (column) {
                case 0:
                    return treatment.isOptimized();
                case 1:
                    return parameter.getName();
                case 2:
                    return parameter.getDescription();
                case 3:
                    return treatment.isOptimized() ? treatment.getTransform() : "held fixed: " + treatment.getReason();
                case 4:
                    if (!hasLowerBound(treatment)) {
                        return "";
                    }
                    Double bound = maximization().getLowerBound(parameter);
                    return bound != null ? bound.toString() : "";
                case 5:
                    if (!treatment.isOptimized()) {
                        return "";
                    }
                    if (treatment == Treatment.NODE_HEIGHTS) {
                        PartitionOptions tree = GradientRegistry.findOwner(options, parameter);
                        return tree instanceof PartitionTreeModel ?
                                GradientRegistry.nodeHeightPlan(options, (PartitionTreeModel) tree).getDescription() : "";
                    }
                    LikelihoodTerm term = likelihoodGradient(parameter);
                    return term + (term == LikelihoodTerm.NUMERIC || term == LikelihoodTerm.NONE ? "" :
                            " (" + term.getDescription() + ")");
            }
            return null;
        }

        public boolean isCellEditable(int row, int column) {
            // only parameters that could be optimized can be ticked or unticked
            if (column == 4) {
                return hasLowerBound(maximization().classify(parameters.get(row), false));
            }
            return column == 0 && maximization().classify(parameters.get(row), true).isOptimized();
        }

        // element-by-element transforms only: the maximizer holds single elements at their bounds
        private boolean hasLowerBound(Treatment treatment) {
            return treatment == Treatment.LOG || treatment == Treatment.LOGIT || treatment == Treatment.NONE;
        }

        public void setValueAt(Object value, int row, int column) {
            if (column == 4) {
                String text = ((String) value).trim();
                try {
                    maximization().setLowerBound(parameters.get(row), text.isEmpty() ? null : Double.parseDouble(text));
                } catch (NumberFormatException e) {
                    Toolkit.getDefaultToolkit().beep();
                    return;
                }
                fireTableRowsUpdated(row, row);
                frame.setDirty();
                return;
            }
            String id = parameters.get(row).getName();
            if ((Boolean) value) {
                maximization().excluded.remove(id);
            } else {
                maximization().excluded.add(id);
            }
            fireTableRowsUpdated(row, row);
            updateStatus();
            frame.setDirty();
        }
    }
}
