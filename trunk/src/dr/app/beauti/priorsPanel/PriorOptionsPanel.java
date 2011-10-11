package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.util.OSType;
import dr.math.distributions.*;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
abstract class PriorOptionsPanel extends OptionsPanel {

    interface Listener {
        void optionsPanelChanged();
    }

    private List<JComponent> argumentFields = new ArrayList<JComponent>();
    private List<String> argumentNames = new ArrayList<String>();

    private boolean isInitializable = true;
    private final boolean isTruncatable;

    private final RealNumberField initialField = new RealNumberField();
    private final JButton negativeInfinityButton;
    private final JButton positiveInfinityButton;

    private final JCheckBox isTruncatedCheck = new JCheckBox("Truncate to:");
    private final RealNumberField lowerField = new RealNumberField();
    private final JLabel lowerLabel = new JLabel("Lower: ");
    private final RealNumberField upperField = new RealNumberField();
    private final JLabel upperLabel = new JLabel("Upper: ");

    protected final Set<Listener> listeners = new HashSet<Listener>();

    PriorOptionsPanel(boolean isTruncatable) {
        super(12, (OSType.isMac() ? 6 : 24));

        this.isTruncatable = isTruncatable;

        negativeInfinityButton = new JButton(NumberFormat.getNumberInstance().format(Double.NEGATIVE_INFINITY));
        PanelUtils.setupComponent(negativeInfinityButton);
        negativeInfinityButton.setFocusable(false);
        negativeInfinityButton.setActionCommand(RealNumberField.NEGATIVE_INFINITY);
        negativeInfinityButton.addActionListener(new ActionListener() {
            public void actionPerformed(final ActionEvent e) {
                lowerField.setText(e.getActionCommand());
            }
        });
        negativeInfinityButton.setToolTipText("Click to set 'Positive Infinity' in the numerical field.");

        positiveInfinityButton = new JButton(NumberFormat.getNumberInstance().format(Double.POSITIVE_INFINITY));
        PanelUtils.setupComponent(positiveInfinityButton);
        positiveInfinityButton.setFocusable(false);
        positiveInfinityButton.setActionCommand(RealNumberField.POSITIVE_INFINITY);
        positiveInfinityButton.addActionListener(new ActionListener() {
            public void actionPerformed(final ActionEvent e) {
                upperField.setText(e.getActionCommand());
            }
        });
        positiveInfinityButton.setToolTipText("Click to set 'Negative Infinity' in the numerical field.");

        initialField.setColumns(10);
        lowerField.setColumns(10);
        upperField.setColumns(10);

        setup();

        isTruncatedCheck.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                lowerField.setEnabled(isTruncatedCheck.isSelected());
                lowerLabel.setEnabled(isTruncatedCheck.isSelected());
                negativeInfinityButton.setEnabled(isTruncatedCheck.isSelected());
                upperField.setEnabled(isTruncatedCheck.isSelected());
                upperLabel.setEnabled(isTruncatedCheck.isSelected());
                positiveInfinityButton.setEnabled(isTruncatedCheck.isSelected());
            }
        });

        KeyListener listener = new KeyAdapter() {
            public void keyReleased(KeyEvent e) {
                if (e.getComponent() instanceof RealNumberField) {
                    String number = ((RealNumberField) e.getComponent()).getText();
                    if (!(number.equals("") || number.endsWith("e") || number.endsWith("E")
                            || number.startsWith("-") || number.endsWith("-"))) {
//                        System.out.println(e.getID() + " = \"" + ((RealNumberField) e.getComponent()).getText() + "\"");
//                        setupChart();
//                        dialog.repaint();
                        for (Listener listener : listeners) {
                            listener.optionsPanelChanged();
                        }
                    }
                }
            }
        };

        initialField.addKeyListener(listener);

        for (JComponent component : argumentFields) {
            if (component instanceof RealNumberField) {
                component.addKeyListener(listener);
            }
        }

        lowerField.addKeyListener(listener);
        upperField.addKeyListener(listener);
    }

    protected RealNumberField getInitialField() {
        return initialField;
    }

    void addListener(Listener listener) {
        listeners.add(listener);
    }

    void removeAllListeners() {
        listeners.clear();
    }
    protected void setFieldRange(RealNumberField field, boolean isNonNegative, boolean isZeroOne) {
        double lower = Double.NEGATIVE_INFINITY;
        double upper = Double.POSITIVE_INFINITY;

        if (isZeroOne) {
            lower = 0.0;
            upper = 1.0;
        } else if (isNonNegative) {
            lower = 0.0;
        }

        field.setRange(lower, upper);
    }

    protected void setFieldRange(RealNumberField field, boolean isNonNegative, boolean isZeroOne, double truncationLower, double truncationUpper) {
        double lower = Double.NEGATIVE_INFINITY;
        double upper = Double.POSITIVE_INFINITY;

        if (isZeroOne) {
            lower = 0.0;
            upper = 1.0;
        } else if (isNonNegative) {
            lower = 0.0;
        }

        if (lower < truncationLower) {
            lower = truncationLower;
        }
        if (upper > truncationUpper) {
            upper = truncationUpper;
        }

        field.setRange(lower, upper);
    }

    protected void addField(String name, double initialValue, double min, double max) {

        RealNumberField field = new RealNumberField(min, max);
        field.setValue(initialValue);
        addField(name, field);
    }

    protected void addField(String name, RealNumberField field) {
        argumentNames.add(name);

        field.setColumns(10);
        argumentFields.add(field);
        setupComponents();
    }

    protected void addCheckBox(String name, JCheckBox jCheckBox) {
        argumentNames.add(name);

        argumentFields.add(jCheckBox);
        setupComponents();
    }

    protected void replaceFieldName(int i, String name) {
        argumentNames.set(i, name);
        setupComponents();
    }

    protected double getValue(int i) {
        return ((RealNumberField) argumentFields.get(i)).getValue();
    }

    private void setupComponents() {
        removeAll();

        if (isInitializable) {
            addComponentWithLabel("Initial value: ", initialField);
        }

        for (int i = 0; i < argumentFields.size(); i++) {
            addComponentWithLabel(argumentNames.get(i) + ":", argumentFields.get(i));
        }

        if (isTruncatable) {
            addSpanningComponent(isTruncatedCheck);
            JPanel panel = new JPanel();
            panel.add(upperField);
            panel.add(positiveInfinityButton);
            addComponents(upperLabel, panel);
             panel = new JPanel();
            panel.add(lowerField);
            panel.add(negativeInfinityButton);
            addComponents(lowerLabel, panel);

            positiveInfinityButton.setMinimumSize(new Dimension(negativeInfinityButton.getWidth(), negativeInfinityButton.getHeight()));

            lowerField.setEnabled(isTruncatedCheck.isSelected());
            lowerLabel.setEnabled(isTruncatedCheck.isSelected());
            negativeInfinityButton.setEnabled(isTruncatedCheck.isSelected());
            upperField.setEnabled(isTruncatedCheck.isSelected());
            upperLabel.setEnabled(isTruncatedCheck.isSelected());
            positiveInfinityButton.setEnabled(isTruncatedCheck.isSelected());
        }
    }

    RealNumberField getField(int i) {
        return (RealNumberField) argumentFields.get(i);
    }

    Distribution getDistribution(Parameter parameter) {
        Distribution dist = getDistribution();

        boolean isBounded = isTruncatedCheck.isSelected();

        double lower = Double.NEGATIVE_INFINITY;
        double upper = Double.POSITIVE_INFINITY;

        if (parameter.isZeroOne) {
            lower = 0.0;
            upper = 1.0;
            isBounded = true;
        } else if (parameter.isNonNegative) {
            lower = 0.0;

            isBounded = true;
        }

        if (dist != null && isBounded) {
            if (isTruncatedCheck.isSelected()) {
                lower = lowerField.getValue();
                upper = upperField.getValue();
            }
            dist = new TruncatedDistribution(dist, lower, upper);
        }
        return dist;
    }

    void setArguments(Parameter parameter, PriorType priorType) {
        if (!parameter.isStatistic) {
            setFieldRange(initialField, parameter.isNonNegative, parameter.isZeroOne);
            initialField.setValue(parameter.initial);
        }
        isTruncatedCheck.setSelected(parameter.isTruncated);
        lowerField.setValue(parameter.getLowerBound());
        upperField.setValue(parameter.getUpperBound());

        setArguments(parameter);

        setupComponents();
    }

    void getArguments(Parameter parameter, PriorType priorType) {
        if (!parameter.isStatistic) {
            parameter.initial = initialField.getValue();
        }
        parameter.isTruncated = isTruncatedCheck.isSelected();
        if (parameter.isTruncated) {
            parameter.truncationLower = lowerField.getValue();
            parameter.truncationUpper = upperField.getValue();
        }

        getArguments(parameter);
    }

    abstract void setup();

    abstract Distribution getDistribution();

    abstract void setArguments(Parameter parameter);

    abstract void getArguments(Parameter parameter);

    static final PriorOptionsPanel UNIFORM = new PriorOptionsPanel(false) {
        void setup() {
            addField("Lower", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Upper", 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        }

        Distribution getDistribution() {
            return new UniformDistribution(getValue(0), getValue(1));
        }

        void setArguments(Parameter parameter) {
            super.setFieldRange(getField(0), parameter.isNonNegative, parameter.isZeroOne);
            super.setFieldRange(getField(1), parameter.isNonNegative, parameter.isZeroOne);

            getField(0).setValue(parameter.getLowerBound());
            getField(1).setValue(parameter.getUpperBound());
        }

        void getArguments(Parameter parameter) {
            parameter.isTruncated = true;
            parameter.truncationLower = getValue(0);
            parameter.truncationUpper = getValue(1);
        }
    };

    static final PriorOptionsPanel EXPONENTIAL = new PriorOptionsPanel(true) {

        void setup() {
            addField("Mean", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new ExponentialDistribution(1.0 / getValue(0)), getValue(1));
        }

        public void setArguments(Parameter parameter) {
            setFieldRange(getField(0), true, parameter.isZeroOne);
            getField(0).setValue(parameter.mean != 0.0 ? parameter.mean : 1.0);
            getField(1).setValue(parameter.offset);
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.offset = getValue(1);
        }
    };

    static final PriorOptionsPanel LAPLACE = new PriorOptionsPanel(true) {
        void setup() {
            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new LaplaceDistribution(getValue(0), getValue(1));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.mean);
            setFieldRange(getField(0), true, false);
            getField(1).setValue(parameter.scale);
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.scale = getValue(1);
        }
    };

    static final PriorOptionsPanel NORMAL = new PriorOptionsPanel(true) {

        void setup() {
            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Stdev", 1.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new NormalDistribution(getValue(0), getValue(1));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.mean);
            getField(1).setValue(parameter.stdev);
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.stdev = getValue(1);
        }
    };

    static final PriorOptionsPanel LOG_NORMAL = new PriorOptionsPanel(true) {
        private JCheckBox meanInRealSpaceCheck;

        void setup() {
            meanInRealSpaceCheck = new JCheckBox();
            if (meanInRealSpaceCheck.isSelected()) {
                addField("Mean", 0.01, 0.0, Double.POSITIVE_INFINITY);
            } else {
                addField("Log(Mean)", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            }
            addField("Log(Stdev)", 1.0, 0.0, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
            addCheckBox("Mean In Real Space", meanInRealSpaceCheck);

            meanInRealSpaceCheck.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent ev) {

                    if (meanInRealSpaceCheck.isSelected()) {
                        replaceFieldName(0, "Mean");
                        if (getValue(0) <= 0) {
                            getField(0).setValue(0.01);
                        }
                        getField(0).setRange(0.0, Double.POSITIVE_INFINITY);
                    } else {
                        replaceFieldName(0, "Log(Mean)");
                        getField(0).setRange(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                    }

                    for (Listener listener : listeners) {
                        listener.optionsPanelChanged();
                    }
                }
            });
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

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.mean);
            getField(1).setValue(parameter.stdev);
            getField(2).setValue(parameter.offset);
            meanInRealSpaceCheck.setSelected(parameter.isMeanInRealSpace());
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.stdev = getValue(1);
            parameter.offset = getValue(2);
            parameter.setMeanInRealSpace(meanInRealSpaceCheck.isSelected());
        }

    };

    static final PriorOptionsPanel GAMMA = new PriorOptionsPanel(true) {

        void setup() {
            addField("Shape", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new GammaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.shape);
            getField(1).setValue(parameter.scale);
            getField(2).setValue(parameter.offset);
        }

        public void getArguments(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.scale = getValue(1);
            parameter.offset = getValue(2);
        }
    };

    static final PriorOptionsPanel INVERSE_GAMMA = new PriorOptionsPanel(true) {

        void setup() {
            addField("Shape", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new InverseGammaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.shape);
            getField(1).setValue(parameter.scale);
            getField(2).setValue(parameter.offset);
        }

        public void getArguments(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.scale = getValue(1);
            parameter.offset = getValue(2);
        }
    };

//    class TruncatedNormalOptionsPanel extends PriorOptionsPanel {
//
//        public TruncatedNormalOptionsPanel() {
//
//            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//            addField("Stdev", 1.0, 0.0, Double.MAX_VALUE);
//            addField("Lower", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//            addField("Upper", 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//        }
//
//        public Distribution getDistribution() {
//            return new TruncatedNormalDistribution(getValue(0), getValue(1), getValue(2), getValue(3));
//        }
//
//        public void setParameterPrior(Parameter parameter) {
//            parameter.mean = getValue(0);
//            parameter.stdev = getValue(1);
//            parameter.isTruncated = true;
//            parameter.truncationLower = getValue(2);
//            parameter.truncationUpper = getValue(3);
//        }
//    }

    static final PriorOptionsPanel BETA = new PriorOptionsPanel(true) {

        void setup() {
            addField("Shape", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("ShapeB", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new BetaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.shape);
            getField(1).setValue(parameter.shapeB);
            getField(2).setValue(parameter.offset);
        }

        public void getArguments(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.shapeB = getValue(1);
            parameter.offset = getValue(2);
        }
    };

    static final PriorOptionsPanel CTMC_RATE_REFERENCE = new PriorOptionsPanel(false) {

        void setup() {
        }

        public Distribution getDistribution() {
            return null;
        }

        public void setArguments(Parameter parameter) {
        }

        public void getArguments(Parameter parameter) {
        }
    };

    static final PriorOptionsPanel ONE_OVER_X = new PriorOptionsPanel(false) {

        void setup() {
        }

        public Distribution getDistribution() {
            return null;
        }

        public void setArguments(Parameter parameter) {
        }

        public void getArguments(Parameter parameter) {
        }
    };

}
