package dr.app.beauti.types;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.util.NumberUtil;
import dr.math.distributions.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public enum PriorType {

    UNDEFINED("undefined", false, true, false, false),
    NONE_TREE_PRIOR("None (Tree Prior Only)", false, true, false, false),
    NONE_STATISTIC("None (Statistic)", false, true, false, false),
    UNIFORM_PRIOR("Uniform", false, false, true, false),
    EXPONENTIAL_PRIOR("Exponential", false, false, true, true),
    LAPLACE_PRIOR("Laplace", false, false, true, true),
    NORMAL_PRIOR("Normal", false, false, true, true),
    LOGNORMAL_PRIOR("Lognormal", false, false, true, true),
    GAMMA_PRIOR("Gamma", false, false, true, true),
    INVERSE_GAMMA_PRIOR("Inverse Gamma", false, false, true, true),
    BETA_PRIOR("Beta", false, false, true, true),
    ONE_OVER_X_PRIOR("1/x", false, false, false, false),
    TRUNC_NORMAL_PRIOR("Truncated Normal", false, false, true, true),
    SUBSTITUTION_REFERENCE_PRIOR("Subst Reference", false, false, false, false),
    NORMAL_HPM_PRIOR("Normal HPM", false, false, false, false),
    LOGNORMAL_HPM_PRIOR("Lognormal HPM", false, false, false, false),
    POISSON_PRIOR("Poisson", true, false, false, false);

    PriorType(final String name, final boolean isDiscrete, final boolean isSpecial, final boolean hasBounds, final boolean hasChart) {
        this.name = name;
        this.isDiscrete = isDiscrete;
        this.isSpecial = isSpecial;
        this.hasBounds = hasBounds;
        this.hasChart = hasChart;
    }

    public String toString() {
        return name;
    }

    public Distribution getDistributionClass(Parameter param) {
        Distribution dist = null;
        switch (this) {
            case UNIFORM_PRIOR:
                dist = new UniformDistribution(param.uniformLower, param.uniformUpper);
                break;
            case EXPONENTIAL_PRIOR:
                if (param.mean == 0) throw new IllegalArgumentException("The mean of exponential prior cannot be 0.");
                dist = new OffsetPositiveDistribution(new ExponentialDistribution(1 / param.mean), param.offset);
                break;
            case LAPLACE_PRIOR:
                dist = new LaplaceDistribution(param.mean, param.scale);
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
            case POISSON_PRIOR:
                dist = new OffsetPositiveDistribution(new PoissonDistribution(param.mean), param.offset);
                break;
            case TRUNC_NORMAL_PRIOR:
                dist = new TruncatedNormalDistribution(param.mean, param.stdev, param.uniformLower, param.uniformUpper);
                break;
            case BETA_PRIOR:
                dist = new OffsetPositiveDistribution(new BetaDistribution(param.shape, param.shapeB), param.offset);
                break;
            // The rest do not give a distribution class (return null)
            case UNDEFINED:
                break;
            case NONE_TREE_PRIOR:
                break;
            case NONE_STATISTIC:
                break;
            case ONE_OVER_X_PRIOR:
                break;
            case SUBSTITUTION_REFERENCE_PRIOR:
                break;
            case NORMAL_HPM_PRIOR:
                break;
            case LOGNORMAL_HPM_PRIOR:
                break;
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
                if (!param.isDiscrete) { // && !param.isStatistic) {
                    buffer.append("Uniform [");
                    buffer.append(NumberUtil.formatDecimal(param.uniformLower, 10, 6));
                    buffer.append(", ");
                    buffer.append(NumberUtil.formatDecimal(param.uniformUpper, 10, 6));
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
                buffer.append(NumberUtil.formatDecimal(param.scale, 10, 6));
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
                buffer.append(NumberUtil.formatDecimal(param.uniformLower, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.uniformUpper, 10, 6));
                buffer.append("]");

                break;
            case BETA_PRIOR:
                buffer.append("Beta [");
                buffer.append(NumberUtil.formatDecimal(param.shape, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(param.shapeB, 10, 6));
                buffer.append("]");
                break;
            case SUBSTITUTION_REFERENCE_PRIOR:
                buffer.append("Approx. Reference Prior");
                break;
            case NORMAL_HPM_PRIOR:
                buffer.append("Normal HPM [mean, precision]");
                break;
            case LOGNORMAL_HPM_PRIOR:
                buffer.append("Lognormal HPM [mean, precision]");
                break;
            default:
                throw new IllegalArgumentException("Unknown prior type");
        }
        if (param.priorType != PriorType.NONE_TREE_PRIOR && (!param.isStatistic) && param.initial != Double.NaN) {
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
            case BETA_PRIOR:
            case TRUNC_NORMAL_PRIOR:
            case SUBSTITUTION_REFERENCE_PRIOR:
            case NORMAL_HPM_PRIOR:
            case LOGNORMAL_HPM_PRIOR:
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

    public String getName() {
        return name;
    }

    public boolean isDiscrete() {
        return isDiscrete;
    }

    public boolean isSpecial() {
        return isSpecial;
    }

    public boolean hasBounds() {
        return hasBounds;
    }

    public boolean hasChart() {
        return hasChart;
    }

    private final String name;
    private final boolean isDiscrete;
    private final boolean isSpecial;
    private final boolean hasBounds;
    private final boolean hasChart;

}
