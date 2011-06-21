/*
 * PriorDialog.java
 *
 * @author Marc A. Suchard
 */

package dr.app.beauti.priorsPanel;

import dr.app.beauti.components.hpm.HierarchicalModelComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;
import dr.app.gui.chart.Axis;
import dr.app.gui.chart.JChart;
import dr.app.gui.chart.LinearAxis;
import dr.app.gui.chart.PDFPlot;
import dr.app.gui.components.RealNumberField;
import dr.app.util.OSType;
import dr.math.distributions.Distribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.OffsetPositiveDistribution;
import dr.util.NumberFormatter;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Marc A. Suchard
 */
public class HierarchicalPriorDialog {

    private JFrame frame;

    private Map<PriorType, PriorOptionsPanel> optionsPanels = new HashMap<PriorType, PriorOptionsPanel>();

    private JComboBox priorCombo = new JComboBox(EnumSet.range(PriorType.NORMAL_HPM_PRIOR, PriorType.LOGNORMAL_HPM_PRIOR).toArray());
    private JComboBox nonPriorCombo = new JComboBox(EnumSet.range(PriorType.NONE_TREE_PRIOR, PriorType.BETA_PRIOR).toArray());
    private JCheckBox meanInRealSpaceCheck = new JCheckBox();
    private RealNumberField initialField = new RealNumberField();
    private RealNumberField selectedField;

    private JTextField nameField;
    private JPanel panel;

//    private final SpecialNumberPanel specialNumberPanel;
    private JChart chart;
    private JPanel quantilePanel;
    private JTextArea quantileText;

    private java.util.List<Parameter> parameterList;
    private Parameter parameter;

    final private BeautiOptions options;

    public HierarchicalPriorDialog(JFrame frame, BeautiOptions options) {
        this.frame = frame;
        this.options = options;

        initialField.setColumns(10);

//        optionsPanels.put(PriorType.UNIFORM_PRIOR, new UniformOptionsPanel());
//        optionsPanels.put(PriorType.LAPLACE_PRIOR, new LaplaceOptionsPanel());
        optionsPanels.put(PriorType.NORMAL_HPM_PRIOR, new NormalOptionsPanel());
//        optionsPanels.put(PriorType.EXPONENTIAL_PRIOR, new ExponentialOptionsPanel());
        optionsPanels.put(PriorType.LOGNORMAL_HPM_PRIOR, new LogNormalOptionsPanel());
//        optionsPanels.put(PriorType.GAMMA_PRIOR, new GammaOptionsPanel());
//        optionsPanels.put(PriorType.INVERSE_GAMMA_PRIOR, new InverseGammaOptionsPanel());
//        optionsPanels.put(PriorType.TRUNC_NORMAL_PRIOR, new TruncatedNormalOptionsPanel());
//        optionsPanels.put(PriorType.BETA_PRIOR, new BetaOptionsPanel());

        chart = new JChart(new LinearAxis(Axis.AT_MINOR_TICK, Axis.AT_MINOR_TICK),
                new LinearAxis(Axis.AT_ZERO, Axis.AT_DATA));

        JLabel quantileLabels = new JLabel();
        quantileLabels.setFont(quantileLabels.getFont().deriveFont(10.0f));
        quantileLabels.setOpaque(false);
        quantileLabels.setText("<html><p align=\"right\">Quantiles: 2.5%:<br>5%:<br>Median:<br>95%:<br>97.5%:</p></html>");

        quantileText = new JTextArea(0, 5);
        quantileText.setFont(quantileText.getFont().deriveFont(10.0f));
        quantileText.setOpaque(false);
        quantileText.setEditable(false);
        quantileLabels.setHorizontalAlignment(JLabel.LEFT);

        quantilePanel = new JPanel();
        quantilePanel.add(quantileLabels);
        quantilePanel.add(quantileText);

//        specialNumberPanel = new SpecialNumberPanel(this);
//        specialNumberPanel.setEnabled(false);
    }


    public void addHPM(java.util.List<Parameter> parameterList) {
        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);
        comp.addHPM(nameField.getText(), parameterList, parameterList.get(0).priorType);
    }

    public boolean validateModelName() {
        return validateModelName(nameField.getText());
    }

    private boolean validateModelName(String modelName) {
        System.err.println("Validating: " + modelName);
        // check that the name is valid
        if (modelName.trim().length() == 0) {
            Toolkit.getDefaultToolkit().beep();
            return false;
        }

//        // disallow a trait called 'date'
//        if (modelName.equalsIgnoreCase("date")) {
//            JOptionPane.showMessageDialog(frame,
//                    "This trait name has a special meaning. Use the 'Tip Date' panel\n" +
//                            " to set dates for taxa.",
//                    "Reserved trait name",
//                    JOptionPane.WARNING_MESSAGE);
//
//            return false;
//        }

        // check that the trait name doesn't exist
        if (modelExists(modelName)) {
            JOptionPane.showMessageDialog(frame,
                    "A model with this name already exists.",
                    "HPM name error",
//                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE);
            System.err.println("Model name exists");
            return false;

//            if (option == JOptionPane.NO_OPTION) {
//                return false;
//            }
        }

        return true;
    }

    private boolean modelExists(String modelName) {

        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);
        return comp.modelExists(modelName);
    }

    public int showDialog(final java.util.List<Parameter> parameterList) {

        this.parameterList = parameterList;
        this.parameter = parameterList.get(0);
        PriorType priorType = parameter.priorType;

        if (parameter.isNodeHeight || parameter.isStatistic) {
            nonPriorCombo.setSelectedItem(priorType);
        } else {
            priorCombo.setSelectedItem(priorType);
        }

        if (!parameter.isStatistic) {
            initialField.setRange(parameter.lower, parameter.upper);
            initialField.setValue(parameter.initial);
        }

        panel = new JPanel(new GridBagLayout());

        setupComponents();

        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        JOptionPane optionPane = new JOptionPane(scrollPane,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Phylogenetic Hierarchical Model Setup");

        priorCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                setupComponents();
                dialog.validate();
                dialog.repaint();
                dialog.pack();
            }
        });

//        nonPriorCombo.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent e) {
//                setupComponents();
//                dialog.validate();
//                dialog.repaint();
//                dialog.pack();
//            }
//        });


//        meanInRealSpaceCheck.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent ev) {
//                PriorOptionsPanel currentPanel = optionsPanels.get(PriorType.LOGNORMAL_PRIOR);
//
//                if (meanInRealSpaceCheck.isSelected()) {
//                    currentPanel.replaceFieldName(0, "Mean");
//                    if (currentPanel.getValue(0) <= 0) {
//                        currentPanel.getField(0).setValue(0.01);
//                    }
//                    currentPanel.getField(0).setRange(0.0, Double.POSITIVE_INFINITY);
//                } else {
//                    currentPanel.replaceFieldName(0, "Log(Mean)");
//                    currentPanel.getField(0).setRange(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//                }
//
//                setupChart();
//                dialog.validate();
//                dialog.repaint();
//                dialog.pack();
//            }
//        });

//        KeyListener listener = new KeyAdapter() {
//            public void keyReleased(KeyEvent e) {
//                if (e.getComponent() instanceof RealNumberField) {
//                    String number = ((RealNumberField) e.getComponent()).getText();
//                    if (!(number.equals("") || number.endsWith("e") || number.endsWith("E")
//                            || number.startsWith("-") || number.endsWith("-"))) {
////                        System.out.println(e.getID() + " = \"" + ((RealNumberField) e.getComponent()).getText() + "\"");
//                        setupChart();
//                        dialog.repaint();
//                    }
//                }
//            }
//        };
//        FocusListener flistener = new FocusAdapter() {
//            public void focusGained(FocusEvent e) {
//                if (e.getComponent() instanceof RealNumberField) {
//                    selectedField = (RealNumberField) e.getComponent();
//                    //specialNumberPanel.setEnabled(true);
//                }
//            }
//
//            public void focusLost(FocusEvent e) {
//                selectedField = null;
//                //specialNumberPanel.setEnabled(false);
//            }
//        };

//        for (PriorOptionsPanel optionsPanel : optionsPanels.values()) {
//            for (JComponent component : optionsPanel.getJComponents()) {
//                if (component instanceof RealNumberField) {
////                    component.addKeyListener(listener);
////                    component.addFocusListener(flistener);
//                }
//            }
//        }

        dialog.pack();
        if (OSType.isMac()) {
            dialog.setMinimumSize(new Dimension(dialog.getBounds().width, 300));
        } else {
            Toolkit tk = Toolkit.getDefaultToolkit();
            Dimension d = tk.getScreenSize();
            if (d.height < 700 && panel.getHeight() > 450) {
                dialog.setSize(new Dimension(panel.getWidth() + 100, 550));
            } else {
                // setSize because optionsPanel is shrunk in dialog
                dialog.setSize(new Dimension(panel.getWidth() + 100, panel.getHeight() + 100));
            }

//            System.out.println("panel width = " + panel.getWidth());
//            System.out.println("panel height = " + panel.getHeight());
        }

        dialog.setResizable(true);
        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

//        if (result == JOptionPane.OK_OPTION) {
//            getArguments();
//        }

        return result;
    }

    private void setArguments(PriorType priorType) {
        PriorOptionsPanel panel;
        switch (priorType) {
            case NORMAL_HPM_PRIOR:
                panel = optionsPanels.get(priorType);
                panel.getField(0).setValue(parameter.mean);
//                panel.getField(1).setRange(0.0, Double.MAX_VALUE);
                panel.getField(1).setValue(parameter.stdev);
                break;

            case LOGNORMAL_HPM_PRIOR:
                panel = optionsPanels.get(priorType);
                if (parameter.isMeanInRealSpace() && parameter.mean <= 0) {// if LOGNORMAL && meanInRealSpace = true, then mean > 0
                    panel.getField(0).setValue(0.01);
                } else {
                    panel.getField(0).setValue(parameter.mean);
                }
                panel.getField(1).setValue(parameter.stdev);
                panel.getField(2).setValue(parameter.offset);
                meanInRealSpaceCheck.setSelected(parameter.isMeanInRealSpace());
                break;
        }
    }

    public void getArguments() {        
        for (Parameter parameter : parameterList) {
            parameter.priorType = (PriorType) priorCombo.getSelectedItem();
            optionsPanels.get(parameter.priorType).setParameterPrior(parameter);
        }
    }

    private void setupComponents() {
        panel.removeAll();

        OptionsPanel optionsPanel = new OptionsPanel(12, (OSType.isMac() ? 6 : 24));

        JPanel panel1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel1.add(optionsPanel);

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 0.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;

        panel.add(panel1, gbc);

        optionsPanel.addSpanningComponent(new JLabel("Select HPM for parameters: "));
        for (Parameter p : parameterList) {
            // TODO Change background color to make these stand out
            optionsPanel.addSpanningComponent(new JLabel("\t" + p.getName()));
        }

        PriorType modelType;
        optionsPanel.addComponentWithLabel("Model Distribution: ", priorCombo);
        modelType = (PriorType) priorCombo.getSelectedItem();

        optionsPanel.addSeparator();

        String modelName = "untitled";
        nameField = new JTextField(modelName);
        nameField.setColumns(20);

//        optionsPanel.addSpanningComponent(new JLabel("Name: "));
//        optionsPanel.addSpanningComponent(nameField);
        
        optionsPanel.addComponentWithLabel("Name: ", nameField);

//        PriorType priorType;
//        if (parameter.isNodeHeight || parameter.isStatistic) {
//            optionsPanel.addComponentWithLabel("Prior Distribution: ", nonPriorCombo);
//            priorType = (PriorType) nonPriorCombo.getSelectedItem();
//            if (priorType == PriorType.NONE_TREE_PRIOR || priorType == PriorType.NONE_STATISTIC) {
//                return;
//            }
//        } else {
//            priorType = (PriorType) priorCombo.getSelectedItem();
//            if (!parameter.priorFixed) {
//                optionsPanel.addComponentWithLabel("Prior Distribution: ", priorCombo);
//            } else {
//                optionsPanel.addComponentWithLabel("Prior Distribution: ", new JLabel(priorType.toString()));
//            }
//        }

//        if (!parameter.isStatistic) {
//            optionsPanel.addSeparator();
//            optionsPanel.addComponentWithLabel("Initial Value: ", initialField);
//        }
//
//        if (parameter.getOptions() instanceof PartitionClockModel) {
//            PartitionClockModel pcm = (PartitionClockModel) parameter.getOptions();
//            initialField.setEnabled(!pcm.getClockModelGroup().isFixMean());
//        }
//
//        if (priorType != PriorType.ONE_OVER_X_PRIOR) {
//            optionsPanel.addSpanningComponent(optionsPanels.get(priorType));
//        }
//
//        if (priorType == PriorType.UNIFORM_PRIOR || priorType == PriorType.TRUNC_NORMAL_PRIOR) {
//            optionsPanel.addSeparator();
//            //optionsPanel.addSpanningComponent(specialNumberPanel);
//        }
//
//        // UNIFORM_PRIOR and ONE_OVER_X_PRIOR have no chart
//        if (priorType != PriorType.UNIFORM_PRIOR && priorType != PriorType.ONE_OVER_X_PRIOR) {
//            optionsPanel.addSeparator();
//
//            setupChart();
//            chart.setPreferredSize(new Dimension(300, 200));
//            chart.setFontSize(8);
//
//            gbc.gridy = 1;
//            gbc.weighty = 1.0;
//            gbc.fill = GridBagConstraints.BOTH;
//
//            panel.add(chart, gbc);
//
//            gbc.gridy = 2;
//            gbc.weighty = 0.0;
//            gbc.anchor = GridBagConstraints.PAGE_END;
//            gbc.fill = GridBagConstraints.HORIZONTAL;
//
//            panel.add(quantilePanel, gbc);
//
//        }
        setArguments(modelType);
        panel.repaint();
    }

    public void setSelectedField(RealNumberField selectedField) {
        this.selectedField = selectedField;
    }

    public RealNumberField getSelectedField() {
        return selectedField;
    }

    NumberFormatter formatter = new NumberFormatter(4);

    private void setupChart() {
        chart.removeAllPlots();

        PriorType priorType;
        if (parameter.isNodeHeight || parameter.isStatistic) {
            priorType = (PriorType) nonPriorCombo.getSelectedItem();
            if (priorType == PriorType.NONE_TREE_PRIOR || priorType == PriorType.NONE_STATISTIC) {
                return;
            }
        } else {
            priorType = (PriorType) priorCombo.getSelectedItem();
        }
        // ExponentialDistribution(1.0 / mean)
//        if (priorType == PriorType.EXPONENTIAL_PRIOR && parameter.mean == 0) parameter.mean = 1;

        double offset = 0.0;
        Distribution distribution = optionsPanels.get(priorType).getDistribution();

        chart.addPlot(new PDFPlot(distribution, offset));
        if (distribution != null) {
            quantileText.setText(formatter.format(distribution.quantile(0.025)) +
                    "\n" + formatter.format(distribution.quantile(0.05)) +
                    "\n" + formatter.format(distribution.quantile(0.5)) +
                    "\n" + formatter.format(distribution.quantile(0.95)) +
                    "\n" + formatter.format(distribution.quantile(0.975)));
        }
    }

    // options panels

    class NormalOptionsPanel extends PriorOptionsPanel {

        public NormalOptionsPanel() {

            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Stdev", 1.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new NormalDistribution(getValue(0), getValue(1));
        }

        public void setParameterPrior(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.stdev = getValue(1);
        }
    }

    class LogNormalOptionsPanel extends PriorOptionsPanel {

        public LogNormalOptionsPanel() {
            if (meanInRealSpaceCheck.isSelected()) {
                addField("Mean", 0.01, 0.0, Double.POSITIVE_INFINITY);
            } else {
                addField("Log(Mean)", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            }
            addField("Log(Stdev)", 1.0, 0.0, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
            addCheckBox("Mean In Real Space", meanInRealSpaceCheck);
        }

        public Distribution getDistribution() {
            double mean = getValue(0);
            if (meanInRealSpaceCheck.isSelected()) {
                if (mean <= 0) {
                    throw new IllegalArgumentException("meanInRealSpace works only for a positive mean");
                }
                mean = Math.log(getValue(0)) - 0.5 * getValue(1) * getValue(1);
            }
            return new OffsetPositiveDistribution(
                    new LogNormalDistribution(mean, getValue(1)), getValue(2));
        }

        public void setParameterPrior(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.stdev = getValue(1);
            parameter.offset = getValue(2);
            parameter.setMeanInRealSpace(meanInRealSpaceCheck.isSelected());
        }

    }
}