/*
 * RestoreInvariantTest.java
 *
 * Copyright (c) 2002-2024 the BEAST Development Team
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
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.xml.unittest;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.util.Identifiable;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.Reportable;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * XML unittest helper that verifies observable values survive a
 * store/propose/restore cycle.
 */
public class RestoreInvariantTest implements Reportable {

    private final Model target;
    private final List<Probe> observables;
    private final List<Probe> touches;
    private final List<Perturbation> perturbations;

    public RestoreInvariantTest(Model target,
                                List<Probe> observables,
                                List<Probe> touches,
                                List<Perturbation> perturbations) {
        this.target = target;
        this.observables = observables;
        this.touches = touches;
        this.perturbations = perturbations;
    }

    @Override
    public String getReport() {

        final List<Parameter> perturbedParameters = uniquePerturbedParameters();
        final Map<Parameter, double[]> originalParameterValues = captureParameterValues(perturbedParameters);

        boolean stored = false;
        Snapshot[] baseline = null;
        Snapshot[] proposal = null;
        Snapshot[] restored = null;
        double maxParameterRestoreDifference = Double.NaN;

        try {
            baseline = snapshotAll(observables);

            target.storeModelState();
            stored = true;

            applyPerturbations(originalParameterValues);
            proposal = snapshotAll(observables);
            snapshotAll(touches);

            target.restoreModelState();
            stored = false;

            restored = snapshotAll(observables);
            maxParameterRestoreDifference = maxParameterDifference(originalParameterValues);

        } finally {
            if (stored && !target.isValidState()) {
                target.restoreModelState();
            }
            if (!parametersMatch(originalParameterValues)) {
                restoreParameterValues(originalParameterValues);
            }
        }

        final double maxAbsDifference = maxDifference(baseline, restored);
        final double proposalMaxAbsDifference = maxDifference(baseline, proposal);

        StringBuilder sb = new StringBuilder();
        sb.append("RestoreInvariantTest\n");
        sb.append("target: ").append(labelFor(target)).append('\n');
        sb.append("observableCount: ").append(observables.size()).append('\n');
        sb.append("touchCount: ").append(touches.size()).append('\n');
        sb.append("perturbationCount: ").append(perturbations.size()).append('\n');
        sb.append("maxAbsDifference: ").append(maxAbsDifference).append('\n');
        sb.append("proposalMaxAbsDifference: ").append(proposalMaxAbsDifference).append('\n');
        sb.append("maxParameterRestoreDifference: ").append(maxParameterRestoreDifference).append('\n');

        for (int i = 0; i < observables.size(); ++i) {
            sb.append("observable[").append(i).append("]: ").append(observables.get(i).getLabel()).append('\n');
            appendSnapshot(sb, "baseline", baseline[i]);
            appendSnapshot(sb, "proposal", proposal[i]);
            appendSnapshot(sb, "restored", restored[i]);
        }

        return sb.toString();
    }

    private Snapshot[] snapshotAll(List<Probe> probes) {
        Snapshot[] snapshots = new Snapshot[probes.size()];
        for (int i = 0; i < probes.size(); ++i) {
            snapshots[i] = probes.get(i).snapshot();
        }
        return snapshots;
    }

    private List<Parameter> uniquePerturbedParameters() {
        List<Parameter> parameters = new ArrayList<Parameter>();
        for (Perturbation perturbation : perturbations) {
            if (!containsIdentity(parameters, perturbation.parameter)) {
                parameters.add(perturbation.parameter);
            }
        }
        return parameters;
    }

    private Map<Parameter, double[]> captureParameterValues(List<Parameter> parameters) {
        Map<Parameter, double[]> values = new IdentityHashMap<Parameter, double[]>();
        for (Parameter parameter : parameters) {
            values.put(parameter, parameter.getParameterValues());
        }
        return values;
    }

    private void applyPerturbations(Map<Parameter, double[]> originalParameterValues) {
        List<Parameter> changedParameters = new ArrayList<Parameter>();

        for (Perturbation perturbation : perturbations) {
            double[] originalValues = originalParameterValues.get(perturbation.parameter);
            if (perturbation.dimension < 0 || perturbation.dimension >= originalValues.length) {
                throw new IllegalArgumentException("Perturbation dimension out of range");
            }

            double value = perturbation.hasAbsoluteValue ?
                    perturbation.value :
                    originalValues[perturbation.dimension] + perturbation.offset;

            perturbation.parameter.setParameterValueQuietly(perturbation.dimension, value);
            if (!containsIdentity(changedParameters, perturbation.parameter)) {
                changedParameters.add(perturbation.parameter);
            }
        }

        for (Parameter parameter : changedParameters) {
            parameter.fireParameterChangedEvent();
        }
    }

    private boolean parametersMatch(Map<Parameter, double[]> originalParameterValues) {
        return maxParameterDifference(originalParameterValues) == 0.0;
    }

    private double maxParameterDifference(Map<Parameter, double[]> originalParameterValues) {
        double max = 0.0;
        for (Map.Entry<Parameter, double[]> entry : originalParameterValues.entrySet()) {
            Parameter parameter = entry.getKey();
            double[] originalValues = entry.getValue();
            for (int i = 0; i < originalValues.length; ++i) {
                max = Math.max(max, difference(originalValues[i], parameter.getParameterValue(i)));
            }
        }
        return max;
    }

    private void restoreParameterValues(Map<Parameter, double[]> originalParameterValues) {
        for (Map.Entry<Parameter, double[]> entry : originalParameterValues.entrySet()) {
            Parameter parameter = entry.getKey();
            double[] values = entry.getValue();
            for (int i = 0; i < values.length; ++i) {
                parameter.setParameterValueQuietly(i, values[i]);
            }
            parameter.fireParameterChangedEvent();
        }
    }

    private static boolean containsIdentity(List<Parameter> parameters, Parameter parameter) {
        for (Parameter candidate : parameters) {
            if (candidate == parameter) {
                return true;
            }
        }
        return false;
    }

    private double maxDifference(Snapshot[] lhs, Snapshot[] rhs) {
        if (lhs == null || rhs == null || lhs.length != rhs.length) {
            return Double.POSITIVE_INFINITY;
        }

        double max = 0.0;
        for (int i = 0; i < lhs.length; ++i) {
            max = Math.max(max, lhs[i].maxDifference(rhs[i]));
        }
        return max;
    }

    private static void appendSnapshot(StringBuilder sb, String label, Snapshot snapshot) {
        sb.append(label).append(": ");
        snapshot.appendTo(sb);
        sb.append('\n');
    }

    private static double difference(double lhs, double rhs) {
        if (Double.isNaN(lhs) || Double.isNaN(rhs)) {
            return Double.isNaN(lhs) && Double.isNaN(rhs) ? 0.0 : Double.POSITIVE_INFINITY;
        }
        if (Double.isInfinite(lhs) || Double.isInfinite(rhs)) {
            return lhs == rhs ? 0.0 : Double.POSITIVE_INFINITY;
        }
        return Math.abs(lhs - rhs);
    }

    private static String labelFor(Object object) {
        if (object instanceof Identifiable) {
            String id = ((Identifiable) object).getId();
            if (id != null) {
                return id;
            }
        }
        return object.getClass().getName();
    }

    private interface Probe {
        Snapshot snapshot();

        String getLabel();
    }

    private static class NumericProbe implements Probe {

        private final Object object;
        private final String label;
        private final ProbeType type;

        private NumericProbe(Object object, String label, ProbeType type) {
            this.object = object;
            this.label = label;
            this.type = type;
        }

        @Override
        public Snapshot snapshot() {
            return new Snapshot(type.snapshot(object));
        }

        @Override
        public String getLabel() {
            return label;
        }
    }

    private static class ReportProbe implements Probe {

        private final Reportable reportable;
        private final String label;
        private final Pattern pattern;

        private ReportProbe(Reportable reportable, String label, String regex) {
            this.reportable = reportable;
            this.label = label;
            this.pattern = regex == null ? null : Pattern.compile(regex, Pattern.DOTALL);
        }

        @Override
        public Snapshot snapshot() {
            String report = reportable.getReport();
            if (pattern != null) {
                Matcher matcher = pattern.matcher(report);
                if (!matcher.find()) {
                    throw new IllegalArgumentException("Report did not match regex for restore invariant probe");
                }
                report = matcher.group(1);
            }
            return new Snapshot(report);
        }

        @Override
        public String getLabel() {
            return label;
        }
    }

    private enum ProbeType {
        PARAMETER("parameter") {
            @Override
            double[] snapshot(Object object) {
                return ((Parameter) object).getParameterValues();
            }
        },
        STATISTIC("statistic") {
            @Override
            double[] snapshot(Object object) {
                Statistic statistic = (Statistic) object;
                double[] values = new double[statistic.getDimension()];
                for (int i = 0; i < values.length; ++i) {
                    values[i] = statistic.getStatisticValue(i);
                }
                return values;
            }
        },
        LIKELIHOOD("likelihood") {
            @Override
            double[] snapshot(Object object) {
                return new double[]{((Likelihood) object).getLogLikelihood()};
            }
        },
        GRADIENT("gradient") {
            @Override
            double[] snapshot(Object object) {
                return ((GradientWrtParameterProvider) object).getGradientLogDensity();
            }
        };

        private final String label;

        ProbeType(String label) {
            this.label = label;
        }

        abstract double[] snapshot(Object object);

        static ProbeType parse(String label) {
            for (ProbeType type : ProbeType.values()) {
                if (type.label.equalsIgnoreCase(label)) {
                    return type;
                }
            }
            throw new IllegalArgumentException("Unknown restore invariant probe type: " + label);
        }
    }

    private static class Snapshot {

        private final double[] values;
        private final String text;

        private Snapshot(double[] values) {
            this.values = new double[values.length];
            System.arraycopy(values, 0, this.values, 0, values.length);
            this.text = null;
        }

        private Snapshot(String text) {
            this.values = null;
            this.text = text;
        }

        double maxDifference(Snapshot other) {
            if (values == null || other.values == null) {
                if (values == null && other.values == null && text.equals(other.text)) {
                    return 0.0;
                }
                return Double.POSITIVE_INFINITY;
            }

            if (values.length != other.values.length) {
                return Double.POSITIVE_INFINITY;
            }

            double max = 0.0;
            for (int i = 0; i < values.length; ++i) {
                max = Math.max(max, difference(values[i], other.values[i]));
            }
            return max;
        }

        void appendTo(StringBuilder sb) {
            if (values == null) {
                sb.append(text);
            } else {
                for (int i = 0; i < values.length; ++i) {
                    if (i > 0) {
                        sb.append(' ');
                    }
                    sb.append(values[i]);
                }
            }
        }
    }

    private static class Perturbation {

        private final Parameter parameter;
        private final int dimension;
        private final double offset;
        private final double value;
        private final boolean hasAbsoluteValue;

        private Perturbation(Parameter parameter,
                             int dimension,
                             double offset,
                             double value,
                             boolean hasAbsoluteValue) {
            this.parameter = parameter;
            this.dimension = dimension;
            this.offset = offset;
            this.value = value;
            this.hasAbsoluteValue = hasAbsoluteValue;
        }
    }

    private static final String PARSER_NAME = "restoreInvariantTest";
    private static final String TARGET = "target";
    private static final String OBSERVABLE = "observable";
    private static final String TOUCH = "touch";
    private static final String PERTURBATION = "perturbation";
    private static final String TYPE = "type";
    private static final String LABEL = "label";
    private static final String REGEX = "regex";
    private static final String DIMENSION = "dim";
    private static final String OFFSET = "offset";
    private static final String VALUE = "value";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Model target = parseTarget(xo.getChild(TARGET));
            List<Probe> observables = parseProbes(xo.getAllChildren(OBSERVABLE));
            List<Probe> touches = parseProbes(xo.getAllChildren(TOUCH));
            List<Perturbation> perturbations = parsePerturbations(xo.getAllChildren(PERTURBATION));

            if (observables.isEmpty()) {
                throw new XMLParseException(PARSER_NAME + " requires at least one observable");
            }
            if (perturbations.isEmpty()) {
                throw new XMLParseException(PARSER_NAME + " requires at least one perturbation");
            }
            validatePerturbationsBelongToTarget(target, perturbations);

            return new RestoreInvariantTest(target, observables, touches, perturbations);
        }

        private Model parseTarget(XMLObject xo) throws XMLParseException {
            Object object = getSingleChild(xo);
            if (object instanceof Model) {
                return (Model) object;
            }
            if (object instanceof Likelihood) {
                return ((Likelihood) object).getModel();
            }
            if (object instanceof GradientWrtParameterProvider) {
                return ((GradientWrtParameterProvider) object).getLikelihood().getModel();
            }
            throw new XMLParseException(PARSER_NAME + " target must be a Model, Likelihood, or GradientWrtParameterProvider");
        }

        private List<Probe> parseProbes(List<XMLObject> wrappers) throws XMLParseException {
            List<Probe> probes = new ArrayList<Probe>();
            for (XMLObject wrapper : wrappers) {
                Object object = getSingleChild(wrapper);
                String label = wrapper.getAttribute(LABEL, labelFor(object));
                String type = wrapper.getAttribute(TYPE, null);
                String regex = wrapper.getAttribute(REGEX, null);

                probes.add(parseProbe(object, label, type, regex));
            }
            return probes;
        }

        private Probe parseProbe(Object object, String label, String type, String regex) throws XMLParseException {
            if (type != null && type.equalsIgnoreCase("report")) {
                if (!(object instanceof Reportable)) {
                    throw new XMLParseException(PARSER_NAME + " report probe requires a Reportable child");
                }
                return new ReportProbe((Reportable) object, label, regex);
            }

            ProbeType probeType = type == null ? inferProbeType(object) : ProbeType.parse(type);
            if (!isCompatibleProbe(object, probeType)) {
                throw new XMLParseException(PARSER_NAME + " " + probeType.label + " probe has incompatible child");
            }
            return new NumericProbe(object, label, probeType);
        }

        private ProbeType inferProbeType(Object object) throws XMLParseException {
            if (object instanceof GradientWrtParameterProvider) {
                return ProbeType.GRADIENT;
            }
            if (object instanceof Likelihood) {
                return ProbeType.LIKELIHOOD;
            }
            if (object instanceof Parameter) {
                return ProbeType.PARAMETER;
            }
            if (object instanceof Statistic) {
                return ProbeType.STATISTIC;
            }
            if (object instanceof Reportable) {
                throw new XMLParseException(PARSER_NAME + " Reportable probes require type=\"report\"");
            }
            throw new XMLParseException(PARSER_NAME + " observable child has unsupported type");
        }

        private boolean isCompatibleProbe(Object object, ProbeType type) {
            switch (type) {
                case PARAMETER:
                    return object instanceof Parameter;
                case STATISTIC:
                    return object instanceof Statistic;
                case LIKELIHOOD:
                    return object instanceof Likelihood;
                case GRADIENT:
                    return object instanceof GradientWrtParameterProvider;
                default:
                    return false;
            }
        }

        private List<Perturbation> parsePerturbations(List<XMLObject> wrappers) throws XMLParseException {
            List<Perturbation> perturbations = new ArrayList<Perturbation>();
            for (XMLObject wrapper : wrappers) {
                Parameter parameter = (Parameter) wrapper.getChild(Parameter.class);
                if (parameter == null) {
                    throw new XMLParseException(PARSER_NAME + " perturbation requires a Parameter child");
                }

                int dimension = wrapper.getIntegerAttribute(DIMENSION);
                boolean hasOffset = wrapper.hasAttribute(OFFSET);
                boolean hasValue = wrapper.hasAttribute(VALUE);
                if (hasOffset == hasValue) {
                    throw new XMLParseException(PARSER_NAME + " perturbation requires exactly one of offset or value");
                }

                double offset = hasOffset ? wrapper.getDoubleAttribute(OFFSET) : 0.0;
                double value = hasValue ? wrapper.getDoubleAttribute(VALUE) : 0.0;
                perturbations.add(new Perturbation(parameter, dimension, offset, value, hasValue));
            }
            return perturbations;
        }

        private void validatePerturbationsBelongToTarget(Model target, List<Perturbation> perturbations)
                throws XMLParseException {
            for (Perturbation perturbation : perturbations) {
                if (!modelContainsVariable(target, perturbation.parameter, new IdentityHashMap<Model, Boolean>())) {
                    throw new XMLParseException(PARSER_NAME + " target does not contain perturbed parameter " +
                            labelFor(perturbation.parameter));
                }
            }
        }

        private boolean modelContainsVariable(Model model, Variable variable, IdentityHashMap<Model, Boolean> visited) {
            if (visited.containsKey(model)) {
                return false;
            }
            visited.put(model, Boolean.TRUE);

            for (int i = 0; i < model.getVariableCount(); ++i) {
                if (model.getVariable(i) == variable) {
                    return true;
                }
            }

            for (int i = 0; i < model.getModelCount(); ++i) {
                if (modelContainsVariable(model.getModel(i), variable, visited)) {
                    return true;
                }
            }

            return false;
        }

        private Object getSingleChild(XMLObject xo) throws XMLParseException {
            if (xo == null) {
                throw new XMLParseException(PARSER_NAME + " has a missing wrapper element");
            }
            if (xo.getChildCount() != 1) {
                throw new XMLParseException(PARSER_NAME + " wrapper element must contain exactly one child");
            }
            return xo.getChild(0);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            XMLSyntaxRule[] singleObject = new XMLSyntaxRule[]{
                    new ElementRule(Object.class)
            };
            XMLSyntaxRule[] probe = new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(TYPE, true),
                    AttributeRule.newStringRule(LABEL, true),
                    AttributeRule.newStringRule(REGEX, true),
                    new ElementRule(Object.class)
            };
            XMLSyntaxRule[] perturbation = new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(DIMENSION),
                    new XORRule(
                            AttributeRule.newDoubleRule(OFFSET),
                            AttributeRule.newDoubleRule(VALUE)
                    ),
                    new ElementRule(Parameter.class)
            };
            return new XMLSyntaxRule[]{
                    new ElementRule(TARGET, singleObject),
                    new ElementRule(OBSERVABLE, probe, 1, Integer.MAX_VALUE),
                    new ElementRule(TOUCH, probe, 0, Integer.MAX_VALUE),
                    new ElementRule(PERTURBATION, perturbation, 1, Integer.MAX_VALUE)
            };
        }

        @Override
        public String getParserDescription() {
            return "Checks that observables return to their baseline values after a store/propose/restore cycle.";
        }

        @Override
        public Class getReturnType() {
            return RestoreInvariantTest.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };
}
