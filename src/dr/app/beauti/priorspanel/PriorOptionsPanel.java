/*
 * PriorOptionsPanel.java
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

package dr.app.beauti.priorspanel;

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
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
abstract class PriorOptionsPanel extends OptionsPanel {

    // this not throw exception, exception thrown by FocusListener in RealNumberField
    public boolean hasInvalidInput(PriorType priorType) { // TODO move all validation here
        for (JComponent component : argumentFields) {
            if (component instanceof RealNumberField && !((RealNumberField) component).isValueValid()) {
                error = ((RealNumberField) component).getErrorMessage();
                return true;
            }
        }
        if (priorType == PriorType.UNIFORM_PRIOR && !isInputValid()) {
            error = "Invalid uniform bound !";
            return true;
        }
        if (isTruncatable && isTruncatedCheck.isSelected()) {
            if (!lowerField.isValueValid()) {
                error = lowerField.getErrorMessage();
                return true;
            } else if (!upperField.isValueValid()) {
                error = upperField.getErrorMessage();
                return true;
            } else if (lowerField.getValue() >= upperField.getValue()) {
                error = "Invalid truncation bound !";
                return true;
            } else if (getValue(OFFSET) > -1 && lowerField.getValue() < getValue(OFFSET)) {
                error = "Offset cannot be smaller than truncation lower !";
                return true;
            } else {
                error = "";
                return false;
            }
        }
        error = "";
        return false;
    }

    public String error = "";

    interface Listener {
        void optionsPanelChanged();
    }

    private List<JComponent> argumentFields = new ArrayList<JComponent>();
    private List<String> argumentNames = new ArrayList<String>();

    private boolean isCalibratedYule = false;
    private boolean isInitializable = true;
    private final boolean isTruncatable;

    private final RealNumberField initialField = new RealNumberField();
    private final JButton negativeInfinityButton;
    private final JButton positiveInfinityButton;

    private final JCheckBox isTruncatedCheck = new JCheckBox("Truncate to:");
    // only this RealNumberField constructor adds FocusListener
    private final RealNumberField lowerField = new RealNumberField(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "truncate lower");
    private final JLabel lowerLabel = new JLabel("Lower: ");
    private final RealNumberField upperField = new RealNumberField(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "truncate upper");
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
                for (Listener listener : listeners) {
                    listener.optionsPanelChanged();
                }
            }
        });

        KeyListener listener = new KeyAdapter() {
            public void keyReleased(KeyEvent e) {
                if (e.getComponent() instanceof RealNumberField) {
                    String number = ((RealNumberField) e.getComponent()).getText();
                    if (!(number.equals("") || number.endsWith("e") || number.endsWith("E") || number.endsWith("-"))) {
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

    protected void addField(String name, double initialValue, double min, boolean includeMin, double max, boolean includeMax) {
        RealNumberField field = new RealNumberField(min, includeMin, max, includeMax, name);
        field.setValue(initialValue);
        addField(name, field);
    }

    protected void addField(String name, double initialValue, double min, double max) {
        RealNumberField field = new RealNumberField(min, max, name);
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
        ((RealNumberField) argumentFields.get(i)).setLabel(name);
        setupComponents();
    }

    protected double getValue(int i) {
        return ((RealNumberField) argumentFields.get(i)).getValue();
    }

    protected double getValue(String fieldName) {
        int i = argumentNames.indexOf(fieldName);
        if (i < 0) return -1; // has no offset field
        return ((RealNumberField) argumentFields.get(i)).getValue();
    }

    protected String getArguName(int i) {
        return argumentNames.get(i);
    }

    private void setupComponents() {
        removeAll();

        if (isInitializable && !isCalibratedYule) {
            addComponentWithLabel("Initial value: ", initialField);
        }

        for (int i = 0; i < argumentFields.size(); i++) {
            addComponentWithLabel(argumentNames.get(i) + ":", argumentFields.get(i));
        }

        if (isTruncatable && !isCalibratedYule) {
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
//            isBounded = true;
        } else if (parameter.isNonNegative) {
            lower = 0.0;

//            isBounded = true;
        }

        if (dist != null && isTruncatable && isBounded) {
            if (isTruncatedCheck.isSelected()) {
                lower = lowerField.getValue();
                upper = upperField.getValue();
            }
            dist = new TruncatedDistribution(dist, lower, upper);
        }
        return dist;
    }

    void setArguments(Parameter parameter, PriorType priorType) {
        this.isCalibratedYule = parameter.isCalibratedYule;
        this.isInitializable = priorType.isInitializable && !parameter.isStatistic && !parameter.isNodeHeight;
        if (!parameter.isStatistic && !parameter.isNodeHeight) {
            setFieldRange(initialField, parameter.isNonNegative, parameter.isZeroOne);
            initialField.setValue(parameter.getInitial());
        }
        isTruncatedCheck.setSelected(parameter.isTruncated);
        setFieldRange(lowerField, parameter.isNonNegative, parameter.isZeroOne, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        setFieldRange(upperField, parameter.isNonNegative, parameter.isZeroOne, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        lowerField.setLabel(parameter.getName() + " truncate lower");
        upperField.setLabel(parameter.getName() + " truncate upper");
        lowerField.setValue(parameter.getLowerBound());
        upperField.setValue(parameter.getUpperBound());

        setArguments(parameter);

        setupComponents();
    }

    void getArguments(Parameter parameter, PriorType priorType) {
        if (priorType.isInitializable && !parameter.isStatistic && !parameter.isNodeHeight) {
            parameter.setInitial(initialField.getValue());
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

    abstract boolean isInputValid();

    private static final String OFFSET = "Offset";

    static final PriorOptionsPanel INFINITE_UNIFORM = new PriorOptionsPanel(false) {
        void setup() {
        }

        Distribution getDistribution() {
            return null;
        }

        void setArguments(Parameter parameter) {
        }

        void getArguments(Parameter parameter) {
        }

        @Override
        boolean isInputValid() {
            return getValue(0) > getValue(1);
        }
    };

    static final PriorOptionsPanel UNIFORM = new PriorOptionsPanel(false) {
        void setup() {
            addField("Upper", 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Lower", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        }

        Distribution getDistribution() {
            return new UniformDistribution(
                    getValue(1), // lower
                    getValue(0) // upper
            );
        }

        void setArguments(Parameter parameter) {
            super.setFieldRange(getField(0), parameter.isNonNegative, parameter.isZeroOne, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            super.setFieldRange(getField(1), parameter.isNonNegative, parameter.isZeroOne, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);

            getField(0).setValue(parameter.uniformUpper);
            getField(1).setValue(parameter.uniformLower);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
        }

        void getArguments(Parameter parameter) {
            parameter.isTruncated = false;
            parameter.uniformUpper = getValue(0);
            parameter.uniformLower = getValue(1);
        }

        @Override
        boolean isInputValid() {
            return getValue(0) > getValue(1);
        }
    };

    static final PriorOptionsPanel EXPONENTIAL = new PriorOptionsPanel(true) {

        void setup() {
            addField("Mean", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField(OFFSET, 0.0, 0.0, true, Double.POSITIVE_INFINITY, true);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new ExponentialDistribution(1.0 / getValue(0)), getValue(1));
        }

        public void setArguments(Parameter parameter) {
            setFieldRange(getField(0), true, parameter.isZeroOne);
            getField(0).setValue(parameter.mean != 0.0 ? parameter.mean : 1.0);
            getField(1).setValue(parameter.offset);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.offset = getValue(1);
        }

        @Override
        boolean isInputValid() {
            return true;
        }
    };

    static final PriorOptionsPanel LAPLACE = new PriorOptionsPanel(true) {
        void setup() {
            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Scale", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
        }

        public Distribution getDistribution() {
            return new LaplaceDistribution(getValue(0), getValue(1));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.mean);
            setFieldRange(getField(0), true, false);
            getField(1).setValue(parameter.scale);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.scale = getValue(1);
        }

        @Override
        boolean isInputValid() {
            return true;
        }
    };

    static final PriorOptionsPanel NORMAL = new PriorOptionsPanel(true) {

        void setup() {
            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            addField("Stdev", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
        }

        public Distribution getDistribution() {
            return new NormalDistribution(getValue(0), getValue(1));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.mean);
            getField(1).setValue(parameter.stdev);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(0);
            parameter.stdev = getValue(1);
        }

        @Override
        boolean isInputValid() {
            return true;
        }
    };

    static final PriorOptionsPanel LOG_NORMAL = new PriorOptionsPanel(true) {
        private JCheckBox parametersInRealSpaceCheck;

        void setup() {
            parametersInRealSpaceCheck = new JCheckBox();
            addCheckBox("Mean/Stdev in real space", parametersInRealSpaceCheck);
            if (parametersInRealSpaceCheck.isSelected()) {
                addField("Mean", 0.01, 0.0, false, Double.POSITIVE_INFINITY, true);
            } else {
                addField("mu", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            }
            if (parametersInRealSpaceCheck.isSelected()) {
                addField("Stdev", 1.0, 0.0, Double.POSITIVE_INFINITY);
            } else {
                addField("sigma", 1.0, 0.0, Double.POSITIVE_INFINITY);
            }
            addField(OFFSET, 0.0, 0.0, Double.POSITIVE_INFINITY);

            parametersInRealSpaceCheck.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent ev) {

                    if (parametersInRealSpaceCheck.isSelected()) {
                        replaceFieldName(1, "Mean");
                        if (getValue(1) <= 0) {
                            getField(1).setValue(0.01);
                        }
                        getField(1).setRange(0.0, Double.POSITIVE_INFINITY);
                        replaceFieldName(2, "Stdev");
                    } else {
                        replaceFieldName(1, "mu");
                        getField(1).setRange(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                        replaceFieldName(2, "sigma");
                    }

                    for (Listener listener : listeners) {
                        listener.optionsPanelChanged();
                    }
                }
            });
        }

        public Distribution getDistribution() {
            double mu = getValue(1);
            double sigma = getValue(2);
            if (parametersInRealSpaceCheck.isSelected()) {
                double mean = getValue(1);
                double stdev = getValue(2);
                if (mean <= 0) {
                    throw new IllegalArgumentException("'Mean in real space' works only for a positive mean");
                }
                mu = Math.log(mean/Math.sqrt(1 + (stdev * stdev) / (mean * mean)));
                sigma = Math.sqrt(Math.log(1 + (stdev * stdev) / (mean * mean)));
            }
            return new OffsetPositiveDistribution(
                    new LogNormalDistribution(mu, sigma),
                    getValue(3));
        }

        public void setArguments(Parameter parameter) {
            getField(1).setValue(parameter.mean);
            getField(2).setValue(parameter.stdev);
            getField(3).setValue(parameter.offset);
            parametersInRealSpaceCheck.setSelected(parameter.isInRealSpace());

            getField(1).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(2).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
            getField(3).setLabel(parameter.getName() + " " + getArguName(2).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.mean = getValue(1);
            parameter.stdev = getValue(2);
            parameter.offset = getValue(3);
            parameter.setMeanInRealSpace(parametersInRealSpaceCheck.isSelected());
        }

        @Override
        boolean isInputValid() {
            return true;
        }

    };

    static final PriorOptionsPanel GAMMA = new PriorOptionsPanel(true) {

        void setup() {
            addField("Shape", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField("Scale", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField(OFFSET, 0.0, 0.0, Double.POSITIVE_INFINITY);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new GammaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.shape);
            getField(1).setValue(parameter.scale);
            getField(2).setValue(parameter.offset);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
            getField(2).setLabel(parameter.getName() + " " + getArguName(2).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.scale = getValue(1);
            parameter.offset = getValue(2);
        }

        @Override
        boolean isInputValid() {
            return true;
        }
    };

    static final PriorOptionsPanel INVERSE_GAMMA = new PriorOptionsPanel(true) {

        void setup() {
            addField("Shape", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField("Scale", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField(OFFSET, 0.0, 0.0, Double.POSITIVE_INFINITY);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new InverseGammaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.shape);
            getField(1).setValue(parameter.scale);
            getField(2).setValue(parameter.offset);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
            getField(2).setLabel(parameter.getName() + " " + getArguName(2).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.scale = getValue(1);
            parameter.offset = getValue(2);
        }

        @Override
        boolean isInputValid() {
            return true;
        }
    };

//    class TruncatedNormalOptionsPanel extends PriorOptionsPanel {
//
//        public TruncatedNormalOptionsPanel() {
//
//            addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//            addField("Stdev", 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
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
            addField("Shape", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField("ShapeB", 1.0, 0.0, false, Double.POSITIVE_INFINITY, true);
            addField(OFFSET, 0.0, 0.0, Double.POSITIVE_INFINITY);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new BetaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setArguments(Parameter parameter) {
            getField(0).setValue(parameter.shape);
            getField(1).setValue(parameter.shapeB);
            getField(2).setValue(parameter.offset);

            getField(0).setLabel(parameter.getName() + " " + getArguName(0).toLowerCase());
            getField(1).setLabel(parameter.getName() + " " + getArguName(1).toLowerCase());
            getField(2).setLabel(parameter.getName() + " " + getArguName(2).toLowerCase());
        }

        public void getArguments(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.shapeB = getValue(1);
            parameter.offset = getValue(2);
        }

        @Override
        boolean isInputValid() {
            return true;
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

        @Override
        boolean isInputValid() {
            return true;
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

        @Override
        boolean isInputValid() {
            return true;
        }
    };

}
