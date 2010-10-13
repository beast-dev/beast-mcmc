package dr.app.beauti.enumTypes;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.util.NumberUtil;
import dr.math.distributions.*;

/**
 * @author Alexei Drummond
 */
public enum PriorType {

    UNDEFINED,
    NONE_TREE_PRIOR,
    NONE_STATISTIC,
    UNIFORM_PRIOR,
    EXPONENTIAL_PRIOR,
    LAPLACE_PRIOR,
    NORMAL_PRIOR,
    LOGNORMAL_PRIOR,
    GAMMA_PRIOR,
    INVERSE_GAMMA_PRIOR,
    ONE_OVER_X_PRIOR,
    TRUNC_NORMAL_PRIOR,
    POISSON_PRIOR;


    public String toString() {

        switch (this) {
            case UNDEFINED:
                return "undefined";
            case NONE_TREE_PRIOR:
                return "None (Tree Prior Only)";
            case NONE_STATISTIC:
                return "None (Statistic)";
            case UNIFORM_PRIOR:
                return "Uniform";
            case EXPONENTIAL_PRIOR:
                return "Exponential";
            case LAPLACE_PRIOR:
                return "Laplace";
            case NORMAL_PRIOR:
                return "Normal";
            case LOGNORMAL_PRIOR:
                return "Lognormal";
            case GAMMA_PRIOR:
                return "Gamma";
            case INVERSE_GAMMA_PRIOR:
                return "Inverse Gamma";
            case ONE_OVER_X_PRIOR:
                return "1/x"; //rename Jeffreys prior to 1/x prior everywhere in Beauti
            case POISSON_PRIOR:
                return "Poisson";
            case TRUNC_NORMAL_PRIOR:
                return "Truncated Normal";
            default:
                return "";
        }
    }

    public Distribution getDistributionClass(Parameter param) {
        Distribution dist = null;
        switch (this) {
            case UNIFORM_PRIOR:
                dist = new UniformDistribution(param.lower, param.upper);
                break;
            case EXPONENTIAL_PRIOR:
                if (param.mean == 0) throw new IllegalArgumentException("The mean of exponential prior cannot be 0.");
                dist = new OffsetPositiveDistribution(new ExponentialDistribution(1/param.mean), param.offset);
                break;
            case LAPLACE_PRIOR:
                dist = new LaplaceDistribution(param.mean, param.stdev);
                break;
            case NORMAL_PRIOR:
                dist = new NormalDistribution(param.mean, param.stdev);
                break;
            case LOGNORMAL_PRIOR:
                dist = new OffsetPositiveDistribution(new LogNormalDistribution(param.mean, param.stdev), param.offset);
                break;
            case GAMMA_PRIOR:
                dist = new OffsetPositiveDistribution(new GammaDistribution(param.shape, param.scale), param.offset);
                break;
            case INVERSE_GAMMA_PRIOR:
                dist = new OffsetPositiveDistribution(new InverseGammaDistribution(param.shape, param.scale), param.offset);
                break;
//            case ONE_OVER_X_PRIOR:
//                return ;
            case POISSON_PRIOR:
                dist = new OffsetPositiveDistribution(new PoissonDistribution(param.mean), param.offset);
                break;
            case TRUNC_NORMAL_PRIOR:
                dist = new TruncatedNormalDistribution(param.mean, param.stdev, param.lower, param.upper);
                break;
//            default: // wrong Exception for other priors without distribution implementation  
//                throw new IllegalArgumentException("Distribution class not available for this prior");
        }
        return dist;
    }

    public String getPriorString(Parameter param) {

//        NumberFormat formatter = NumberFormat.getNumberInstance();
        StringBuffer buffer = new StringBuffer();

        if (param.priorType == PriorType.UNDEFINED) {
            buffer.append("? ");
        } else if (param.isPriorImproper()) {
            buffer.append("! ");
        } else if (!param.isPriorEdited()) {
            buffer.append("* ");
        } else {
            buffer.append("  ");
        }

        switch (param.priorType) {
            case NONE_TREE_PRIOR:
                buffer.append("Using Tree Prior");
                break;
            case NONE_STATISTIC:
                buffer.append("Indirectly Specified Through Other Parameter");
                break;
            case UNDEFINED:
                buffer.append("Not yet specified");
                break;
            case UNIFORM_PRIOR:
                if (!param.isDiscrete && !param.isStatistic) {
                    buffer.append("Uniform [");
                    buffer.append(NumberUtil.formatDecimal(param.lower, 10, 6));
                    buffer.append(", ");
                    buffer.append(NumberUtil.formatDecimal(param.upper, 10, 6));
                    buffer.append("]");
                } else {
                    buffer.append("Uniform");
                }
                break;
            case EXPONENTIAL_PRIOR:
                buffer.append("Exponential [");
                buffer.append(NumberUtil.formatDecimal(param.mean, 10, 6));
                buffer.append("]");
                break;
            case LAPLACE_PRIOR:
                buffer.append("Laplace [");
                buffer.append(NumberUtil.formatDecimal(param.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.stdev, 10, 6));
                buffer.append("]");
                break;
            case NORMAL_PRIOR:
                buffer.append("Normal [");
                buffer.append(NumberUtil.formatDecimal(param.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.stdev, 10, 6));
                buffer.append("]");
                break;
            case LOGNORMAL_PRIOR:
                buffer.append("LogNormal [");
                if (param.isMeanInRealSpace()) buffer.append("R");
                buffer.append(NumberUtil.formatDecimal(param.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.stdev, 10, 6));
                buffer.append("]");
                break;
            case GAMMA_PRIOR:
                buffer.append("Gamma [");
                buffer.append(NumberUtil.formatDecimal(param.shape, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.scale, 10, 6));
                buffer.append("]");
                break;
            case INVERSE_GAMMA_PRIOR:
                buffer.append("Inverse Gamma [");
                buffer.append(NumberUtil.formatDecimal(param.shape, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.scale, 10, 6));
                buffer.append("]");
                break;
            case ONE_OVER_X_PRIOR:
                buffer.append("1/x"); // rename Jeffreys prior to 1/x prior everywhere in Beauti
                break;
            case POISSON_PRIOR:
                buffer.append("Poisson [");
                buffer.append(NumberUtil.formatDecimal(param.mean, 10, 6));
                buffer.append("]");
                break;
            case TRUNC_NORMAL_PRIOR:
                buffer.append("Truncated Normal [");
                buffer.append(NumberUtil.formatDecimal(param.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.stdev, 10, 6));
                buffer.append("]");
                buffer.append(" in [");
                buffer.append(NumberUtil.formatDecimal(param.lower, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.upper, 10, 6));
                buffer.append("]");

                break;
            default:
                throw new IllegalArgumentException("Unknown prior type");
        }
        if (param.priorType != PriorType.NONE_TREE_PRIOR && !param.isStatistic) {
            buffer.append(", initial=").append(NumberUtil.formatDecimal(param.initial, 10, 6));
        }

        return buffer.toString();
    }

    public String getPriorBoundString(Parameter param) {

        if (param.isStatistic) {
            return "n/a";
        }

//        NumberFormat formatter = NumberFormat.getNumberInstance();
        StringBuffer buffer = new StringBuffer();

        switch (param.priorType) {
            case NONE_TREE_PRIOR:
            case NONE_STATISTIC:
//                buffer.append("None");
//                break;
            case UNDEFINED:
            case UNIFORM_PRIOR:
            case EXPONENTIAL_PRIOR:
            case LAPLACE_PRIOR:
            case NORMAL_PRIOR:
            case LOGNORMAL_PRIOR:
            case GAMMA_PRIOR:
            case INVERSE_GAMMA_PRIOR:
            case ONE_OVER_X_PRIOR:
            case POISSON_PRIOR:
            case TRUNC_NORMAL_PRIOR:
                buffer.append("[");
                buffer.append(NumberUtil.formatDecimal(param.lower, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.upper, 10, 6));
                buffer.append("]");
                break;
            default:
                throw new IllegalArgumentException("Unknown prior type");
        }
        return buffer.toString();
    }
}
