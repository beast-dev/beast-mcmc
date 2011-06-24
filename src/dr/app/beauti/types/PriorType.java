package dr.app.beauti.types;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.util.NumberUtil;
import dr.math.distributions.*;
import dr.util.DataTable;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public enum PriorType {

    UNDEFINED("undefined", false, true, false, false, true),
    NONE_TREE_PRIOR("None (Tree Prior Only)", false, true, false, false, true),
    NONE_STATISTIC("None (Statistic)", false, true, false, false, true),
    UNIFORM_PRIOR("Uniform", false, false, true, false, true),
    EXPONENTIAL_PRIOR("Exponential", false, false, true, true, true),
    LAPLACE_PRIOR("Laplace", false, false, true, true, true),
    NORMAL_PRIOR("Normal", false, false, true, true, true),
    LOGNORMAL_PRIOR("Lognormal", false, false, true, true, true),
    GAMMA_PRIOR("Gamma", false, false, true, true, true),
    INVERSE_GAMMA_PRIOR("Inverse Gamma", false, false, true, true, true),
    BETA_PRIOR("Beta", false, false, true, true, true),
    ONE_OVER_X_PRIOR("1/x", false, false, false, false, true),
    TRUNC_NORMAL_PRIOR("Truncated Normal", false, false, true, true, true),
    CMTC_RATE_REFERENCE_PRIOR("Subst Reference", false, false, false, false, false),
    LOGNORMAL_HPM_PRIOR("Lognormal HPM", false, false, false, false, false),
    NORMAL_HPM_PRIOR("Normal HPM", false, false, false, false, false),
    POISSON_PRIOR("Poisson", true, false, false, false, true);

    PriorType(final String name, final boolean isDiscrete, final boolean isSpecial,
              final boolean hasBounds, final boolean hasChart, final boolean displayByDefault) {
        this.name = name;
        this.isDiscrete = isDiscrete;
        this.isSpecial = isSpecial;
        this.hasBounds = hasBounds;
        this.hasChart = hasChart;
        this.displayByDefault = displayByDefault;
    }

    public String toString() {
        return name;
    }

    public Distribution getDistributionClass(Parameter param) {
        Distribution dist = null;
        switch (this) {
            case UNIFORM_PRIOR:
                dist = new UniformDistribution(param.truncationLower, param.truncationUpper);
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
                dist = new TruncatedNormalDistribution(param.mean, param.stdev, param.truncationLower, param.truncationUpper);
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
            case CMTC_RATE_REFERENCE_PRIOR:
                break;
            case NORMAL_HPM_PRIOR:
                break;
            case LOGNORMAL_HPM_PRIOR:
                break;
        }
        return dist;
    }

    public String getPriorString(Parameter parameter) {

//        NumberFormat formatter = NumberFormat.getNumberInstance();
        StringBuffer buffer = new StringBuffer();

        if (parameter.priorType == PriorType.UNDEFINED) {
            buffer.append("? ");
        } else if (parameter.isPriorImproper()) {
            buffer.append("! ");
        } else if (!parameter.isPriorEdited()) {
            buffer.append("* ");
        } else {
            buffer.append("  ");
        }

        double lower = parameter.getLowerBound();
        double upper = parameter.getUpperBound();

        switch (parameter.priorType) {
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
                if (!parameter.isDiscrete) { // && !param.isStatistic) {
                    buffer.append("Uniform [");
                    buffer.append(NumberUtil.formatDecimal(lower, 10, 6));
                    buffer.append(", ");
                    buffer.append(NumberUtil.formatDecimal(upper, 10, 6));
                    buffer.append("]");
                } else {
                    buffer.append("Uniform");
                }
                break;
            case EXPONENTIAL_PRIOR:
                buffer.append("Exponential [");
                buffer.append(NumberUtil.formatDecimal(parameter.mean, 10, 6));
                buffer.append("]");
                break;
            case LAPLACE_PRIOR:
                buffer.append("Laplace [");
                buffer.append(NumberUtil.formatDecimal(parameter.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(parameter.scale, 10, 6));
                buffer.append("]");
                break;
            case NORMAL_PRIOR:
            case TRUNC_NORMAL_PRIOR:
                buffer.append("Normal [");
                buffer.append(NumberUtil.formatDecimal(parameter.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(parameter.stdev, 10, 6));
                buffer.append("]");
                break;
            case LOGNORMAL_PRIOR:
                buffer.append("LogNormal [");
                if (parameter.isMeanInRealSpace()) buffer.append("R");
                buffer.append(NumberUtil.formatDecimal(parameter.mean, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(parameter.stdev, 10, 6));
                buffer.append("]");
                break;
            case GAMMA_PRIOR:
                buffer.append("Gamma [");
                buffer.append(NumberUtil.formatDecimal(parameter.shape, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(parameter.scale, 10, 6));
                buffer.append("]");
                break;
            case INVERSE_GAMMA_PRIOR:
                buffer.append("Inverse Gamma [");
                buffer.append(NumberUtil.formatDecimal(parameter.shape, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(parameter.scale, 10, 6));
                buffer.append("]");
                break;
            case ONE_OVER_X_PRIOR:
                buffer.append("1/x"); // rename Jeffreys prior to 1/x prior everywhere in Beauti
                break;
            case POISSON_PRIOR:
                buffer.append("Poisson [");
                buffer.append(NumberUtil.formatDecimal(parameter.mean, 10, 6));
                buffer.append("]");
                break;
            case BETA_PRIOR:
                buffer.append("Beta [");
                buffer.append(NumberUtil.formatDecimal(parameter.shape, 10, 6));
                buffer.append(", ");
                buffer.append(NumberUtil.formatDecimal(parameter.shapeB, 10, 6));
                buffer.append("]");
                break;
            case CMTC_RATE_REFERENCE_PRIOR:
                buffer.append("Approx. Reference Prior");
                break;
            case NORMAL_HPM_PRIOR:
                buffer.append("Normal HPM [");
                buffer.append(parameter.hpmModelName);
                buffer.append("]");
                break;
            case LOGNORMAL_HPM_PRIOR:
                buffer.append("Lognormal HPM [");
                buffer.append(parameter.hpmModelName);
                buffer.append("]");
                break;
            default:
                throw new IllegalArgumentException("Unknown prior type");
        }
        if (parameter.priorType != UNIFORM_PRIOR && parameter.isTruncated) {
            buffer.append(" in [");
            buffer.append(NumberUtil.formatDecimal(parameter.truncationLower, 10, 6));
            buffer.append(", ");
            buffer.append(NumberUtil.formatDecimal(parameter.truncationUpper, 10, 6));
            buffer.append("]");
        }


        if (parameter.priorType != PriorType.NONE_TREE_PRIOR && (!parameter.isStatistic) && parameter.initial != Double.NaN) {
            buffer.append(", initial=").append(NumberUtil.formatDecimal(parameter.initial, 10, 6));
        }

        return buffer.toString();
    }

    public String getPriorBoundString(Parameter parameter) {

        if (parameter.isStatistic) {
            return "n/a";
        }

        double lower = parameter.getLowerBound();
        double upper = parameter.getUpperBound();

//        NumberFormat formatter = NumberFormat.getNumberInstance();
        StringBuffer buffer = new StringBuffer();
        buffer.append("[");
        buffer.append(NumberUtil.formatDecimal(lower, 10, 6));
        buffer.append(", ");
        buffer.append(NumberUtil.formatDecimal(upper, 10, 6));
        buffer.append("]");

        return buffer.toString();
    }

    public static PriorType[] getPriorTypes(Parameter parameter) {
        if (parameter.isDiscrete) {
            return new PriorType[] {UNIFORM_PRIOR, POISSON_PRIOR};
        } else if (parameter.isCMTCRate) {
            return new PriorType[] {
                    CMTC_RATE_REFERENCE_PRIOR,
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    NORMAL_PRIOR,
                    TRUNC_NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    INVERSE_GAMMA_PRIOR,
                    ONE_OVER_X_PRIOR};
        } else if (parameter.isHierarchical) {
            return new PriorType[] {
                    LOGNORMAL_HPM_PRIOR,
                    NORMAL_HPM_PRIOR};
        } else if (parameter.isZeroOne) {
            return new PriorType[] {
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    NORMAL_PRIOR,
                    TRUNC_NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    INVERSE_GAMMA_PRIOR,
                    BETA_PRIOR};
        } else if (parameter.isNonNegative) {
            return new PriorType[] {
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    LAPLACE_PRIOR,
                    NORMAL_PRIOR,
                    TRUNC_NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    INVERSE_GAMMA_PRIOR,
                    ONE_OVER_X_PRIOR};
        } else { // just a continuous parameter
            return new PriorType[] {
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    LAPLACE_PRIOR,
                    NORMAL_PRIOR,
                    TRUNC_NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    INVERSE_GAMMA_PRIOR};
        }
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

    public boolean displayByDefault() {
        return displayByDefault;
    }

    private final String name;
    private final boolean isDiscrete;
    private final boolean isSpecial;
    private final boolean hasBounds;
    private final boolean hasChart;
    private final boolean displayByDefault;

}
