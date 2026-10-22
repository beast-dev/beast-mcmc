/*
 * SimpleOperatorScheduleParser.java
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

package dr.inferencexml.operators;

import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.xml.*;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 *
 */
public class SimpleOperatorScheduleParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_SCHEDULE = "operators";
    public static final String SEQUENTIAL = "sequential";
    public static final String OPTIMIZATION_SCHEDULE = "optimizationSchedule";

    public static final String ACCEPTANCE_THRESHOLD = "minAcceptance";
    public static final String USE_THRESHOLD = "minUsage";

    /**
     * System property (set by the -operator_weights command line option) holding a comma separated
     * list of key=factor pairs that multiply the XML weights. A key matches an operator's full name
     * (e.g. "scale(clock.rate)") or the part before the first '(' (e.g. "subtreeLeap"); "*" matches
     * any operator not otherwise matched. A factor of 0 removes the operator from the schedule.
     */
    public static final String OPERATOR_WEIGHTS_PROPERTY = "mcmc.operator_weights";

    public String getParserName() {
        return OPERATOR_SCHEDULE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int useThreshold = xo.getAttribute(USE_THRESHOLD, 1000);
        double acceptanceThreshold = xo.getAttribute(ACCEPTANCE_THRESHOLD, 0.0);

        SimpleOperatorSchedule schedule = new SimpleOperatorSchedule(useThreshold, acceptanceThreshold);

        if (xo.hasAttribute(SEQUENTIAL)) {
            schedule.setSequential(xo.getBooleanAttribute(SEQUENTIAL));
        }

        Logger.getLogger("dr.inference").info("\nCreating operator scheduler");

        if (xo.hasAttribute(OPTIMIZATION_SCHEDULE)) {
            String type = xo.getStringAttribute(OPTIMIZATION_SCHEDULE);
            Logger.getLogger("dr.inference").info("  Optimization schedule: " + type);

            try {
                if (type.equalsIgnoreCase("default")) {
                    schedule.setOptimizationTransform(OperatorSchedule.DEFAULT_TRANSFORM);
                } else {
                    schedule.setOptimizationTransform(OperatorSchedule.OptimizationTransform.valueOf(type.toUpperCase()));
                }
            } catch (IllegalArgumentException iae) {
                throw new RuntimeException("Unsupported optimization schedule");
            }
        }

        Map<String, Double> weightFactors = parseWeightFactors(System.getProperty(OPERATOR_WEIGHTS_PROPERTY));

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof MCMCOperator) {
                MCMCOperator operator = (MCMCOperator) child;
                if (weightFactors != null && !applyWeightFactor(operator, weightFactors)) {
                    continue;
                }
                schedule.addOperator(operator);
            }
        }
        return schedule;
    }

    private static Map<String, Double> parseWeightFactors(String spec) {
        if (spec == null || spec.trim().isEmpty()) {
            return null;
        }
        Map<String, Double> factors = new LinkedHashMap<>();
        for (String entry : spec.split(",")) {
            String[] parts = entry.split("=");
            if (parts.length != 2) {
                throw new IllegalArgumentException("Operator weight override should be key=factor: " + entry);
            }
            double factor = Double.parseDouble(parts[1].trim());
            if (factor < 0.0) {
                throw new IllegalArgumentException("Operator weight factor must be non-negative: " + entry);
            }
            factors.put(parts[0].trim(), factor);
        }
        return factors;
    }

    /**
     * Multiplies the operator's weight by the matching factor.
     * @return false if the factor is zero and the operator should be left out of the schedule
     */
    private static boolean applyWeightFactor(MCMCOperator operator, Map<String, Double> factors) {
        String name = operator.getOperatorName();
        int bracket = name.indexOf('(');
        String prefix = (bracket > 0 ? name.substring(0, bracket) : name).trim();

        Double factor = factors.get(name);
        if (factor == null) {
            factor = factors.get(prefix);
        }
        if (factor == null) {
            factor = factors.get("*");
        }
        if (factor == null || factor == 1.0) {
            return true;
        }

        Logger logger = Logger.getLogger("dr.inference");
        if (factor == 0.0) {
            logger.info("  Operator weight override: removing " + name);
            return false;
        }
        double weight = operator.getWeight() * factor;
        logger.info("  Operator weight override: " + name + " weight " + operator.getWeight() + " -> " + weight);
        operator.setWeight(weight);
        return true;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SEQUENTIAL, true),
            new ElementRule(MCMCOperator.class, 1, Integer.MAX_VALUE),
            AttributeRule.newStringRule(OPTIMIZATION_SCHEDULE, true),
            AttributeRule.newDoubleRule(ACCEPTANCE_THRESHOLD, true, "Acceptance rate below which an operator will be switched off"),
            AttributeRule.newIntegerRule(USE_THRESHOLD, true, "Minimum number of usage before testing acceptance threshold")
    };

    public String getParserDescription() {
        return "A simple operator scheduler";
    }

    public Class getReturnType() {
        return SimpleOperatorSchedule.class;
    }
    
}
