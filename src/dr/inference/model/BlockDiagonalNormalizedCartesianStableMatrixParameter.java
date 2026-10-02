package dr.inference.model;

import dr.xml.*;

/**
 * 2x2 block chart for the normalized Cartesian stable parametrization
 *
 *     [ rate       rate * ratio - t ]
 *     [ rate * ratio + t   rate    ]
 *
 * with an optional leading 1x1 scalar block in odd dimension.
 *
 * In the manuscript drift convention the corresponding block is negated, so
 * rate is the positive selection-strength counterpart of alpha < 0.
 */
public final class BlockDiagonalNormalizedCartesianStableMatrixParameter
        extends AbstractBlockDiagonalTwoByTwoMatrixParameter {

    public static final String NAME = "blockDiagonalNormalizedCartesianStableMatrixParameter";

    private static final String R_ELEMENT = "rotationMatrix";
    private static final String SCALAR_ELEMENT = "scalarBlock";
    private static final String RATE_ELEMENT = "blockRate";
    private static final String RATIO_ELEMENT = "blockRatio";
    private static final String T_ELEMENT = "blockT";
    private static final String RATE_ORDERING = "blockRateOrdering";

    private final RateOrdering rateOrdering;

    public BlockDiagonalNormalizedCartesianStableMatrixParameter(final String name,
                                                                 final MatrixParameter RParam,
                                                                 final Parameter scalarBlockParam,
                                                                 final Parameter rateParam,
                                                                 final Parameter ratioParam,
                                                                 final Parameter tParam) {
        this(name, RParam, scalarBlockParam, rateParam, ratioParam, tParam, RateOrdering.NONE);
    }

    public BlockDiagonalNormalizedCartesianStableMatrixParameter(final String name,
                                                                 final MatrixParameter RParam,
                                                                 final Parameter scalarBlockParam,
                                                                 final Parameter rateParam,
                                                                 final Parameter ratioParam,
                                                                 final Parameter tParam,
                                                                 final RateOrdering rateOrdering) {
        super(name, RParam, scalarBlockParam, rateParam, ratioParam, tParam);
        this.rateOrdering = rateOrdering == null ? RateOrdering.NONE : rateOrdering;
    }

    @Override
    protected int getTwoByTwoBlockParameterCount() {
        return 3;
    }

    @Override
    protected String getTwoByTwoBlockParameterName(final int k) {
        switch (k) {
            case 0: return "rate";
            case 1: return "ratio";
            case 2: return "t";
            default: throw new IllegalArgumentException("Invalid parameter index: " + k);
        }
    }

    @Override
    protected void fillTwoByTwoBlock(final int blockIndex,
                                     final double[] outBlock) {
        final double rate = getEffectiveRate(blockIndex);
        final double ratio = getTwoByTwoBlockParameterValue(1, blockIndex);
        final double t = getTwoByTwoBlockParameterValue(2, blockIndex);
        final double symmetricOffDiagonal = rate * ratio;

        outBlock[0] = rate;
        outBlock[1] = symmetricOffDiagonal - t;
        outBlock[2] = symmetricOffDiagonal + t;
        outBlock[3] = rate;
    }

    @Override
    protected void chainTwoByTwoBlockGradient(final int blockIndex,
                                              final double g00,
                                              final double g01,
                                              final double g10,
                                              final double g11,
                                              final double[] out,
                                              final int baseOffset) {
        final double rate = getEffectiveRate(blockIndex);
        final double ratio = getTwoByTwoBlockParameterValue(1, blockIndex);
        final double offDiagonalGradientSum = g01 + g10;

        final double gradRate = g00 + g11 + ratio * offDiagonalGradientSum;
        final double gradRatio = rate * offDiagonalGradientSum;
        final double gradT = -g01 + g10;

        addRateGradient(out, baseOffset, blockIndex, gradRate);
        out[twoByTwoGradientOffset(baseOffset, 1, blockIndex)] = gradRatio;
        out[twoByTwoGradientOffset(baseOffset, 2, blockIndex)] = gradT;
    }

    private double getEffectiveRate(final int blockIndex) {
        switch (rateOrdering) {
            case NONE:
                return getTwoByTwoBlockParameterValue(0, blockIndex);
            case ASCENDING:
                return sumRateIncrements(0, blockIndex + 1);
            case DESCENDING:
                return sumRateIncrements(blockIndex, num2x2Blocks);
            default:
                throw new IllegalArgumentException("Unhandled block rate ordering: " + rateOrdering);
        }
    }

    private double sumRateIncrements(final int fromInclusive,
                                     final int toExclusive) {
        double sum = 0.0;
        for (int i = fromInclusive; i < toExclusive; i++) {
            final double increment = getTwoByTwoBlockParameterValue(0, i);
            if (increment < 0.0) {
                throw new ArithmeticException("Ordered block rate increments must be non-negative");
            }
            sum += increment;
        }
        return sum;
    }

    private void addRateGradient(final double[] out,
                                 final int baseOffset,
                                 final int blockIndex,
                                 final double gradient) {
        switch (rateOrdering) {
            case NONE:
                out[twoByTwoGradientOffset(baseOffset, 0, blockIndex)] = gradient;
                return;
            case ASCENDING:
                for (int i = 0; i <= blockIndex; i++) {
                    out[twoByTwoGradientOffset(baseOffset, 0, i)] += gradient;
                }
                return;
            case DESCENDING:
                for (int i = blockIndex; i < num2x2Blocks; i++) {
                    out[twoByTwoGradientOffset(baseOffset, 0, i)] += gradient;
                }
                return;
            default:
                throw new IllegalArgumentException("Unhandled block rate ordering: " + rateOrdering);
        }
    }

    public Parameter getScalarBlockParameter() {
        return super.getScalarBlockParameter();
    }

    public Parameter getRateParameter() {
        return getTwoByTwoBlockParameter(0);
    }

    public RateOrdering getRateOrdering() {
        return rateOrdering;
    }

    public Parameter getRatioParameter() {
        return getTwoByTwoBlockParameter(1);
    }

    public Parameter getTParameter() {
        return getTwoByTwoBlockParameter(2);
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;
            final MatrixParameter R = (MatrixParameter) xo.getElementFirstChild(R_ELEMENT);
            final Parameter rateParam = (Parameter) xo.getElementFirstChild(RATE_ELEMENT);
            final Parameter ratioParam = (Parameter) xo.getElementFirstChild(RATIO_ELEMENT);
            final Parameter tParam = (Parameter) xo.getElementFirstChild(T_ELEMENT);
            final RateOrdering rateOrdering =
                    RateOrdering.parse(xo.getAttribute(RATE_ORDERING, RateOrdering.NONE.name));

            final boolean oddDimension = (R.getRowDimension() & 1) == 1;

            final Parameter scalarBlockParam;
            if (xo.hasChildNamed(SCALAR_ELEMENT)) {
                scalarBlockParam = (Parameter) xo.getElementFirstChild(SCALAR_ELEMENT);
            } else if (oddDimension) {
                throw new XMLParseException("Missing required element '" + SCALAR_ELEMENT + "' for odd matrix dimension.");
            } else {
                scalarBlockParam = new Parameter.Default(0);
            }

            return new BlockDiagonalNormalizedCartesianStableMatrixParameter(
                    name, R, scalarBlockParam, rateParam, ratioParam, tParam, rateOrdering
            );
        }

        @Override
        public String getParserDescription() {
            return "A matrix parameter A = R D R^{-1} with D composed of one optional 1x1 block " +
                    "and 2x2 stable blocks parametrized by normalized Cartesian coordinates (rate, ratio, t).";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return BlockDiagonalNormalizedCartesianStableMatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(RATE_ORDERING, true,
                    "Optional block rate ordering: none, ascending, or descending. " +
                            "Ascending/descending interpret blockRate values as non-negative cumulative increments."),
            new ElementRule(R_ELEMENT, new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
            new ElementRule(SCALAR_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(RATE_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RATIO_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(T_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };

    public enum RateOrdering {
        NONE("none"),
        ASCENDING("ascending"),
        DESCENDING("descending");

        private final String name;

        RateOrdering(final String name) {
            this.name = name;
        }

        private static RateOrdering parse(final String value) throws XMLParseException {
            for (final RateOrdering ordering : values()) {
                if (ordering.name.equalsIgnoreCase(value)) {
                    return ordering;
                }
            }
            throw new XMLParseException("Unknown " + RATE_ORDERING + " value '" + value +
                    "'. Expected none, ascending, or descending.");
        }
    }
}
