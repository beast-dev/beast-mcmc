/*
 * OperatorAnalysisPrinter.java
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

package dr.inference.operators;

import dr.util.NumberFormatter;

import java.io.PrintStream;

/**
 * Package: OperatorAnalysisPrinter
 * Description:
 * Factoring out OperatorAnalysis from MCMC class. Maybe someone will depricate use in MCMC
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Jan 12, 2010
 *         Time: 6:02:30 PM
 */
public class OperatorAnalysisPrinter {

    private static final NumberFormatter formatter = new NumberFormatter(8);

    /**
     * Writes ano operator analysis to the provided print stream
     *
     * @param out the print stream to write operator analysis to
     */
    public static void showOperatorAnalysis(PrintStream out, OperatorSchedule schedule, boolean useAdaptation) {
        out.println();
        out.println("Operator analysis");
        out.println(formatter.formatToFieldWidth("Operator", 50) +
                formatter.formatToFieldWidth("Tuning", 9) +
                formatter.formatToFieldWidth("Count", 11) +
                formatter.formatToFieldWidth("Time", 9) +
                formatter.formatToFieldWidth("Time/Op", 9) +
                formatter.formatToFieldWidth("Pr(accept)", 11) +
                (useAdaptation ? "" : " Performance suggestion"));

        for (int i = 0; i < schedule.getOperatorCount(); i++) {

            final MCMCOperator op = schedule.getOperator(i);
            if (op instanceof JointOperator) {
                JointOperator jointOp = (JointOperator) op;
                for (int k = 0; k < jointOp.getNumberOfSubOperators(); k++) {
                    out.println(formattedOperatorName(jointOp.getSubOperatorName(k))
                                    + formattedParameterString(jointOp.getSubOperator(k))
                                    + formattedCountString(op)
                                    + formattedTimeString(op)
                                    + formattedTimePerOpString(op)
                                    + formattedProbString(jointOp)
                                    + (useAdaptation ? "" : formattedDiagnostics(jointOp, jointOp.getAcceptanceProbability()))
                    );
                }
            } else {
                out.println(formattedOperatorName(op.getOperatorName())
                                + formattedParameterString(op)
                                + formattedCountString(op)
                                + formattedTimeString(op)
                                + formattedTimePerOpString(op)
                                + formattedProbString(op)
                                + (useAdaptation ? "" : formattedDiagnostics(op, op.getAcceptanceProbability()))
                );
            }

        }
        out.println();
    }

    private static String formattedOperatorName(String operatorName) {
        return formatter.formatToFieldWidth(operatorName, 50);
    }

    private static String formattedParameterString(MCMCOperator op) {
        String pString = "        ";
        if (op instanceof AdaptableMCMCOperator && ((AdaptableMCMCOperator) op).getMode() != AdaptationMode.ADAPTATION_OFF) {
            pString = formatter.formatToFieldWidth(formatter.formatDecimal(((AdaptableMCMCOperator) op).getRawParameter(), 3), 8);
        }
        return pString;
    }

    private static String formattedCountString(MCMCOperator op) {
        final long count = op.getCount();
        return formatter.formatToFieldWidth(Long.toString(count), 10) + " ";
    }

    private static String formattedTimeString(MCMCOperator op) {
        final long time = op.getTotalEvaluationTime();
        return formatter.formatToFieldWidth(Long.toString(time), 8) + " ";
    }

    private static String formattedTimePerOpString(MCMCOperator op) {
        final double time = op.getMeanEvaluationTime();
        return formatter.formatToFieldWidth(formatter.formatDecimal(time, 2), 8) + " ";
    }

    private static String formattedProbString(MCMCOperator op) {
        final double acceptanceProb = op.getAcceptanceProbability();
        return formatter.formatToFieldWidth(formatter.formatDecimal(acceptanceProb, 4), 11) + " ";
    }

    private static String formattedDiagnostics(MCMCOperator op, double acceptanceProb) {

        String message = "good";
        if (acceptanceProb < AbstractAdaptableOperator.MINIMUM_GOOD_ACCEPTANCE_LEVEL) {
            if (acceptanceProb < (AbstractAdaptableOperator.MINIMUM_ACCEPTANCE_LEVEL / 10.0)) {
                message = "very low";
            } else if (acceptanceProb < AbstractAdaptableOperator.MINIMUM_ACCEPTANCE_LEVEL) {
                message = "low";
            } else message = "slightly low";

        } else if (acceptanceProb > AbstractAdaptableOperator.MAXIMUM_GOOD_ACCEPTANCE_LEVEL) {
            double reallyHigh = 1.0 - ((1.0 - AbstractAdaptableOperator.MAXIMUM_ACCEPTANCE_LEVEL) / 10.0);
            if (acceptanceProb > reallyHigh) {
                message = "very high";
            } else if (acceptanceProb > AbstractAdaptableOperator.MAXIMUM_ACCEPTANCE_LEVEL) {
                message = "high";
            } else message = "slightly high";
        }

        String performancesMsg = "";
        if (op instanceof GibbsOperator) {
            performancesMsg = "none (Gibbs operator)";
        } else if (op instanceof AdaptableMCMCOperator){
            final String suggestion = ((AdaptableMCMCOperator)op).getPerformanceSuggestion();
            performancesMsg = message + "\t" + suggestion;
        }

        return performancesMsg;
    }
}
