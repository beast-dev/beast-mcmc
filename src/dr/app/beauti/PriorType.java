package dr.app.beauti;

import java.text.NumberFormat;

/**
 * @author Alexei Drummond
 */
public enum PriorType {

    NONE,
    UNIFORM_PRIOR,
    EXPONENTIAL_PRIOR,
    NORMAL_PRIOR,
    LOGNORMAL_PRIOR,
    GAMMA_PRIOR,
    JEFFREYS_PRIOR,
    POISSON_PRIOR;

    public String toString() {

        switch (this) {
            case NONE: return "none";
            case UNIFORM_PRIOR: return "Uniform";
            case EXPONENTIAL_PRIOR: return "Exponential";
            case NORMAL_PRIOR: return "Normal";
            case LOGNORMAL_PRIOR: return "Lognormal";
            case GAMMA_PRIOR: return "Gamma";
            case JEFFREYS_PRIOR: return "Jeffreys'";
            case POISSON_PRIOR: return "Poisson";
            default: return "";
        }
    }

    public String getPriorString(BeastGenerator.Parameter param) {

        NumberFormat formatter = NumberFormat.getNumberInstance();
        StringBuffer buffer = new StringBuffer();

        if (!param.priorEdited) {
            buffer.append("* ");
        }
        switch (param.priorType) {
            case NONE:
                buffer.append("Using Tree Prior");
                break;
            case UNIFORM_PRIOR:
                if (!param.isDiscrete && !param.isStatistic) {
                    buffer.append("Uniform [");
                    buffer.append(formatter.format(param.uniformLower));
                    buffer.append(", ");
                    buffer.append(formatter.format(param.uniformUpper));
                    buffer.append("]");
                } else {
                    buffer.append("Uniform");
                }
                break;
            case EXPONENTIAL_PRIOR:
                buffer.append("Exponential [");
                buffer.append(formatter.format(param.exponentialMean));
                buffer.append("]");
                break;
            case NORMAL_PRIOR:
                buffer.append("Normal [");
                buffer.append(formatter.format(param.normalMean));
                buffer.append(", ");
                buffer.append(formatter.format(param.normalStdev));
                buffer.append("]");
                break;
            case LOGNORMAL_PRIOR:
                buffer.append("LogNormal [");
                buffer.append(formatter.format(param.logNormalMean));
                buffer.append(", ");
                buffer.append(formatter.format(param.logNormalStdev));
                buffer.append("]");
                break;
            case GAMMA_PRIOR:
                buffer.append("Gamma [");
                buffer.append(formatter.format(param.gammaAlpha));
                buffer.append(", ");
                buffer.append(formatter.format(param.gammaBeta));
                buffer.append("]");
                break;
            case JEFFREYS_PRIOR:
                buffer.append("Jeffreys");
                break;
            case POISSON_PRIOR:
                buffer.append("Poisson [");
                buffer.append(formatter.format(param.poissonMean));
                buffer.append("]");
                break;
            default:
                throw new IllegalArgumentException("Unknown prior type");
        }
        if (param.priorType != PriorType.NONE && !param.isStatistic) {
            buffer.append(", initial=").append(param.initial);
        }

        return buffer.toString();
    }
}
