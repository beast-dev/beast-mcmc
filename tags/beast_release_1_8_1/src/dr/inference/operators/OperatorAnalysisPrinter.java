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

    private final OperatorSchedule schedule;
    private final NumberFormatter formatter = new NumberFormatter(8);

    public OperatorAnalysisPrinter(OperatorSchedule schedule) {
        this.schedule = schedule;
    }

    public void showOperatorAnalysis(PrintStream out) {
        out.println();
        out.println("Operator analysis");
        out.println(formatter.formatToFieldWidth("Operator", 50) +
                formatter.formatToFieldWidth("Tuning", 9) +
                formatter.formatToFieldWidth("Count", 11) +
                formatter.formatToFieldWidth("Time", 9) +
                formatter.formatToFieldWidth("Time/Op", 9) +
                formatter.formatToFieldWidth("Pr(accept)", 11) +
                " Performance suggestion");

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
                            + formattedDiagnostics(jointOp, MCMCOperator.Utils.getAcceptanceProbability(jointOp)));
                }
            } else {
                out.println(formattedOperatorName(op.getOperatorName())
                        + formattedParameterString(op)
                        + formattedCountString(op)
                        + formattedTimeString(op)
                        + formattedTimePerOpString(op)
                        + formattedProbString(op)
                        + formattedDiagnostics(op, MCMCOperator.Utils.getAcceptanceProbability(op)));
            }

        }
        out.println();
    }

    private String formattedOperatorName(String operatorName) {
        return formatter.formatToFieldWidth(operatorName, 50);
    }

    private String formattedParameterString(MCMCOperator op) {
        String pString = "        ";
        if (op instanceof CoercableMCMCOperator && ((CoercableMCMCOperator) op).getMode() != CoercionMode.COERCION_OFF) {
            pString = formatter.formatToFieldWidth(formatter.formatDecimal(((CoercableMCMCOperator) op).getRawParameter(), 3), 8);
        }
        return pString;
    }

    private String formattedCountString(MCMCOperator op) {
        final int count = op.getCount();
        return formatter.formatToFieldWidth(Integer.toString(count), 10) + " ";
    }

    private String formattedTimeString(MCMCOperator op) {
        final long time = op.getTotalEvaluationTime();
        return formatter.formatToFieldWidth(Long.toString(time), 8) + " ";
    }

    private String formattedTimePerOpString(MCMCOperator op) {
        final double time = op.getMeanEvaluationTime();
        return formatter.formatToFieldWidth(formatter.formatDecimal(time, 2), 8) + " ";
    }

    private String formattedProbString(MCMCOperator op) {
        final double acceptanceProb = MCMCOperator.Utils.getAcceptanceProbability(op);
        return formatter.formatToFieldWidth(formatter.formatDecimal(acceptanceProb, 4), 11) + " ";
    }

    private String formattedDiagnostics(MCMCOperator op, double acceptanceProb) {

        String message = "good";
        if (acceptanceProb < op.getMinimumGoodAcceptanceLevel()) {
            if (acceptanceProb < (op.getMinimumAcceptanceLevel() / 10.0)) {
                message = "very low";
            } else if (acceptanceProb < op.getMinimumAcceptanceLevel()) {
                message = "low";
            } else message = "slightly low";

        } else if (acceptanceProb > op.getMaximumGoodAcceptanceLevel()) {
            double reallyHigh = 1.0 - ((1.0 - op.getMaximumAcceptanceLevel()) / 10.0);
            if (acceptanceProb > reallyHigh) {
                message = "very high";
            } else if (acceptanceProb > op.getMaximumAcceptanceLevel()) {
                message = "high";
            } else message = "slightly high";
        }

        String performacsMsg;
        if (op instanceof GibbsOperator) {
            performacsMsg = "none (Gibbs operator)";
        } else {
            final String suggestion = op.getPerformanceSuggestion();
            performacsMsg = message + "\t" + suggestion;
        }

        return performacsMsg;
    }
}
