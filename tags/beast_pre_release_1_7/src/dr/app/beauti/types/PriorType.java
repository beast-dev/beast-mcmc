package dr.app.beauti.types;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.util.NumberUtil;
import dr.math.distributions.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public enum PriorType {

    UNDEFINED("undefined", false, false),
    NONE_TREE_PRIOR("None (Tree Prior Only)", false, false),
    NONE_STATISTIC("None (Statistic)", false, false),
    NONE_IMPROPER("Infinite Uniform (Improper)", false, false),
    UNIFORM_PRIOR("Uniform", false, false),
    EXPONENTIAL_PRIOR("Exponential", true, true),
    LAPLACE_PRIOR("Laplace", true, true),
    NORMAL_PRIOR("Normal", true, true),
    LOGNORMAL_PRIOR("Lognormal", true, true),
    GAMMA_PRIOR("Gamma", true, true),
    INVERSE_GAMMA_PRIOR("Inverse Gamma", true, true),
    BETA_PRIOR("Beta", true, true),
    ONE_OVER_X_PRIOR("1/x", true, false),
    CMTC_RATE_REFERENCE_PRIOR("CTMC Rate Reference", false, false),
    LOGNORMAL_HPM_PRIOR("Lognormal HPM", false, false),
    NORMAL_HPM_PRIOR("Normal HPM", false, false),
    POISSON_PRIOR("Poisson", false, false);

    PriorType(final String name, final boolean isTruncatable, final boolean isPlottable) {
        this.name = name;
        this.isTruncatable = isTruncatable;
        this.isPlottable = isPlottable;
    }

    public String toString() {
        return name;
    }

    public Distribution getDistributionInstance(Parameter parameter) {
        Distribution dist = null;
        switch (this) {
            case UNIFORM_PRIOR:
                dist = new UniformDistribution(parameter.getLowerBound(), parameter.getUpperBound());
                break;
            case EXPONENTIAL_PRIOR:
                if (parameter.mean == 0)
                    throw new IllegalArgumentException("The mean of exponential prior cannot be 0.");
                dist = new OffsetPositiveDistribution(new ExponentialDistribution(1.0 / parameter.mean), parameter.offset);
                break;
            case LAPLACE_PRIOR:
                dist = new LaplaceDistribution(parameter.mean, parameter.scale);
                break;
            case NORMAL_PRIOR:
                dist = new NormalDistribution(parameter.mean, parameter.stdev);
                break;
            case LOGNORMAL_PRIOR:
                dist = new OffsetPositiveDistribution(new LogNormalDistribution(parameter.mean, parameter.stdev), parameter.offset);
                break;
            case GAMMA_PRIOR:
                dist = new OffsetPositiveDistribution(new GammaDistribution(parameter.shape, parameter.scale), parameter.offset);
                break;
            case INVERSE_GAMMA_PRIOR:
                dist = new OffsetPositiveDistribution(new InverseGammaDistribution(parameter.shape, parameter.scale), parameter.offset);
                break;
            case POISSON_PRIOR:
                dist = new OffsetPositiveDistribution(new PoissonDistribution(parameter.mean), parameter.offset);
                break;
            case BETA_PRIOR:
                dist = new OffsetPositiveDistribution(new BetaDistribution(parameter.shape, parameter.shapeB), parameter.offset);
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
        if (dist != null && parameter.isTruncated) {
            dist = new TruncatedDistribution(dist, parameter.getLowerBound(), parameter.getUpperBound());
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
            case NONE_IMPROPER:
                buffer.append("Uniform infinite bounds");
                break;
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
        if (parameter.isTruncated) {
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
            return new PriorType[]{
                    UNIFORM_PRIOR,
                    POISSON_PRIOR};
        }
        if (parameter.isNodeHeight) {
            if (parameter.isCalibratedYule) {
                return new PriorType[]{
                        NONE_TREE_PRIOR,
                        UNIFORM_PRIOR,
                        EXPONENTIAL_PRIOR,
                        NORMAL_PRIOR,
                        LOGNORMAL_PRIOR,
                        GAMMA_PRIOR};
            } else {
                return new PriorType[]{
                        NONE_TREE_PRIOR,
                        UNIFORM_PRIOR,
                        EXPONENTIAL_PRIOR,
                        LAPLACE_PRIOR,
                        NORMAL_PRIOR,
                        LOGNORMAL_PRIOR,
                        GAMMA_PRIOR,
                        INVERSE_GAMMA_PRIOR,
                        ONE_OVER_X_PRIOR};
            }
        }
        if (parameter.isStatistic) {
            if (parameter.isCalibratedYule) {
                return new PriorType[]{
                        NONE_STATISTIC,
                        UNIFORM_PRIOR,
                        EXPONENTIAL_PRIOR,
                        NORMAL_PRIOR,
                        LOGNORMAL_PRIOR,
                        GAMMA_PRIOR};
            } else {
                return new PriorType[]{
                        NONE_STATISTIC,
                        UNIFORM_PRIOR,
                        EXPONENTIAL_PRIOR,
                        LAPLACE_PRIOR,
                        NORMAL_PRIOR,
                        LOGNORMAL_PRIOR,
                        GAMMA_PRIOR,
                        INVERSE_GAMMA_PRIOR,
                        ONE_OVER_X_PRIOR};
            }
        }
        if (parameter.isCMTCRate) {
            return new PriorType[]{
                    NONE_IMPROPER,
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    CMTC_RATE_REFERENCE_PRIOR,
                    INVERSE_GAMMA_PRIOR,
                    ONE_OVER_X_PRIOR};
        }
        if (parameter.isHierarchical) {
            return new PriorType[]{
                    LOGNORMAL_HPM_PRIOR,
                    NORMAL_HPM_PRIOR};
        }
        if (parameter.isZeroOne) {
            return new PriorType[]{
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    INVERSE_GAMMA_PRIOR,
                    BETA_PRIOR};
        }
        if (parameter.isNonNegative) {
            return new PriorType[]{
                    NONE_IMPROPER,
                    UNIFORM_PRIOR,
                    EXPONENTIAL_PRIOR,
                    LAPLACE_PRIOR,
                    NORMAL_PRIOR,
                    LOGNORMAL_PRIOR,
                    GAMMA_PRIOR,
                    INVERSE_GAMMA_PRIOR,
                    ONE_OVER_X_PRIOR};
        }

        // just a continuous parameter
        return new PriorType[]{
                NONE_IMPROPER,
                UNIFORM_PRIOR,
                EXPONENTIAL_PRIOR,
                LAPLACE_PRIOR,
                NORMAL_PRIOR,
                LOGNORMAL_PRIOR,
                GAMMA_PRIOR,
                INVERSE_GAMMA_PRIOR};

    }

    public String getName() {
        return name;
    }

    public boolean isTruncatable() {
        return isTruncatable;
    }

    public boolean isPlottable() {
        return isPlottable;
    }

    private final String name;
    private final boolean isTruncatable;
    private final boolean isPlottable;

}
