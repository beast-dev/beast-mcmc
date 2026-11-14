package dr.inference.model;

import dr.xml.*;

/**
 * 2x2 block chart for the polar stable parametrization
 *
 *     [ rho cos(theta)   rho sin(theta) - t ]
 *     [ rho sin(theta) + t   rho cos(theta) ]
 *
 * with an optional leading 1x1 scalar block in odd dimension.
 *
 * Under the new superclass:
 * - scalarBlockParam stores the optional leading 1x1 block
 *   (dimension 1 if matrix dimension is odd, 0 if even),
 * - rhoParam, thetaParam, tParam each have dimension num2x2Blocks.
 */
public final class BlockDiagonalPolarStableMatrixParameter extends AbstractBlockDiagonalTwoByTwoMatrixParameter {

    public static final String NAME = "blockDiagonalPolarStableMatrixParameter";

    private static final String R_ELEMENT = "rotationMatrix";
    private static final String SCALAR_ELEMENT = "scalarBlock";
    private static final String RHO_ELEMENT = "blockRho";
    private static final String THETA_ELEMENT = "blockTheta";
    private static final String T_ELEMENT = "blockT";
    private static final String RHO_ORDERING = "blockRhoOrdering";

    private final RhoOrdering rhoOrdering;

    public BlockDiagonalPolarStableMatrixParameter(final String name,
                                                   final MatrixParameter RParam,
                                                   final Parameter scalarBlockParam,
                                                   final Parameter rhoParam,
                                                   final Parameter thetaParam,
                                                   final Parameter tParam) {
        this(name, RParam, scalarBlockParam, rhoParam, thetaParam, tParam, RhoOrdering.NONE);
    }

    public BlockDiagonalPolarStableMatrixParameter(final String name,
                                                   final MatrixParameter RParam,
                                                   final Parameter scalarBlockParam,
                                                   final Parameter rhoParam,
                                                   final Parameter thetaParam,
                                                   final Parameter tParam,
                                                   final RhoOrdering rhoOrdering) {
        super(name, RParam, scalarBlockParam, rhoParam, thetaParam, tParam);
        this.rhoOrdering = rhoOrdering == null ? RhoOrdering.NONE : rhoOrdering;
    }

    @Override
    protected int getTwoByTwoBlockParameterCount() {
        return 3;
    }

    @Override
    protected String getTwoByTwoBlockParameterName(final int k) {
        switch (k) {
            case 0: return "rho";
            case 1: return "theta";
            case 2: return "t";
            default: throw new IllegalArgumentException("Invalid parameter index: " + k);
        }
    }

    @Override
    protected void fillTwoByTwoBlock(final int blockIndex,
                                     final double[] outBlock) {
        final double rho = getEffectiveRho(blockIndex);
        final double theta = getTwoByTwoBlockParameterValue(1, blockIndex);
        final double t = getTwoByTwoBlockParameterValue(2, blockIndex);

        final double cos = Math.cos(theta);
        final double sin = Math.sin(theta);

        outBlock[0] = rho * cos;
        outBlock[1] = rho * sin - t;
        outBlock[2] = rho * sin + t;
        outBlock[3] = rho * cos;
    }

    @Override
    protected void chainTwoByTwoBlockGradient(final int blockIndex,
                                              final double g00,
                                              final double g01,
                                              final double g10,
                                              final double g11,
                                              final double[] out,
                                              final int baseOffset) {
        final double rho = getEffectiveRho(blockIndex);
        final double theta = getTwoByTwoBlockParameterValue(1, blockIndex);

        final double cos = Math.cos(theta);
        final double sin = Math.sin(theta);

        // d00 = rho cos(theta)
        // d01 = rho sin(theta) - t
        // d10 = rho sin(theta) + t
        // d11 = rho cos(theta)

        final double gradRho =
                g00 * cos +
                        g01 * sin +
                        g10 * sin +
                        g11 * cos;

        final double gradTheta =
                g00 * (-rho * sin) +
                        g01 * ( rho * cos) +
                        g10 * ( rho * cos) +
                        g11 * (-rho * sin);

        final double gradT =
                -g01 + g10;

        addRhoGradient(out, baseOffset, blockIndex, gradRho);
        out[twoByTwoGradientOffset(baseOffset, 1, blockIndex)] = gradTheta;
        out[twoByTwoGradientOffset(baseOffset, 2, blockIndex)] = gradT;
    }

    private double getEffectiveRho(final int blockIndex) {
        switch (rhoOrdering) {
            case NONE:
                return getTwoByTwoBlockParameterValue(0, blockIndex);
            case ASCENDING:
                return sumRhoIncrements(0, blockIndex + 1);
            case DESCENDING:
                return sumRhoIncrements(blockIndex, num2x2Blocks);
            default:
                throw new IllegalArgumentException("Unhandled block rho ordering: " + rhoOrdering);
        }
    }

    private double sumRhoIncrements(final int fromInclusive,
                                    final int toExclusive) {
        double sum = 0.0;
        for (int i = fromInclusive; i < toExclusive; i++) {
            final double increment = getTwoByTwoBlockParameterValue(0, i);
            if (increment < 0.0) {
                throw new ArithmeticException("Ordered block rho increments must be non-negative");
            }
            sum += increment;
        }
        return sum;
    }

    private void addRhoGradient(final double[] out,
                                final int baseOffset,
                                final int blockIndex,
                                final double gradient) {
        switch (rhoOrdering) {
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
                throw new IllegalArgumentException("Unhandled block rho ordering: " + rhoOrdering);
        }
    }

    public Parameter getScalarBlockParameter() {
        return super.getScalarBlockParameter();
    }

    public Parameter getRhoParameter() {
        return getTwoByTwoBlockParameter(0);
    }

    public RhoOrdering getRhoOrdering() {
        return rhoOrdering;
    }

    public Parameter getThetaParameter() {
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
            final Parameter rhoParam = (Parameter) xo.getElementFirstChild(RHO_ELEMENT);
            final Parameter thetaParam = (Parameter) xo.getElementFirstChild(THETA_ELEMENT);
            final Parameter tParam = (Parameter) xo.getElementFirstChild(T_ELEMENT);
            final RhoOrdering rhoOrdering = RhoOrdering.parse(xo.getAttribute(RHO_ORDERING, RhoOrdering.NONE.name));

            final boolean oddDimension = (R.getRowDimension() & 1) == 1;

            final Parameter scalarBlockParam;
            if (xo.hasChildNamed(SCALAR_ELEMENT)) {
                scalarBlockParam = (Parameter) xo.getElementFirstChild(SCALAR_ELEMENT);
            } else if (oddDimension) {
                throw new XMLParseException("Missing required element '" + SCALAR_ELEMENT + "' for odd matrix dimension.");
            } else {
                scalarBlockParam = new Parameter.Default(0);
            }

            return new BlockDiagonalPolarStableMatrixParameter(
                    name, R, scalarBlockParam, rhoParam, thetaParam, tParam, rhoOrdering
            );
        }

        @Override
        public String getParserDescription() {
            return "A matrix parameter A = R D R^{-1} with D composed of one optional 1x1 block " +
                    "and 2x2 stable blocks parametrized by (rho, theta, t).";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return BlockDiagonalPolarStableMatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(RHO_ORDERING, true,
                    "Optional block rho ordering: none, ascending, or descending. " +
                            "Ascending/descending interpret blockRho values as non-negative cumulative increments."),
            new ElementRule(R_ELEMENT, new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
            new ElementRule(SCALAR_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(RHO_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(THETA_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(T_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };

    public enum RhoOrdering {
        NONE("none"),
        ASCENDING("ascending"),
        DESCENDING("descending");

        private final String name;

        RhoOrdering(final String name) {
            this.name = name;
        }

        private static RhoOrdering parse(final String value) throws XMLParseException {
            for (final RhoOrdering ordering : values()) {
                if (ordering.name.equalsIgnoreCase(value)) {
                    return ordering;
                }
            }
            throw new XMLParseException("Unknown " + RHO_ORDERING + " value '" + value +
                    "'. Expected none, ascending, or descending.");
        }
    }
}
