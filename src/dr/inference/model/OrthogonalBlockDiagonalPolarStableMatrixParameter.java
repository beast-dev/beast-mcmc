package dr.inference.model;

import dr.xml.*;

/**
 * True orthogonal variant of the polar stable block-diagonal parameterization.
 *
 * <p>The basis matrix is supplied through an {@link OrthogonalMatrixProvider}, so
 * the user parameterizes the change of basis through orthogonal chart coordinates
 * rather than free matrix entries.</p>
 */
public final class OrthogonalBlockDiagonalPolarStableMatrixParameter
        extends AbstractBlockDiagonalTwoByTwoMatrixParameter {

    public static final String NAME = "orthogonalBlockDiagonalPolarStableMatrixParameter";

    private static final String R_ELEMENT = "orthogonalRotationMatrix";
    private static final String SCALAR_ELEMENT = "scalarBlock";
    private static final String RHO_ELEMENT = "blockRho";
    private static final String THETA_ELEMENT = "blockTheta";
    private static final String T_ELEMENT = "blockT";
    private static final String RHO_ORDERING = "blockRhoOrdering";

    private final CompoundParameter nativeCompoundParameter;
    private final RhoOrdering rhoOrdering;

    public OrthogonalBlockDiagonalPolarStableMatrixParameter(final String name,
                                                             final MatrixParameter rotationParam,
                                                             final Parameter scalarBlockParam,
                                                             final Parameter rhoParam,
                                                             final Parameter thetaParam,
                                                             final Parameter tParam) {
        this(name, rotationParam, scalarBlockParam, rhoParam, thetaParam, tParam, RhoOrdering.NONE);
    }

    public OrthogonalBlockDiagonalPolarStableMatrixParameter(final String name,
                                                             final MatrixParameter rotationParam,
                                                             final Parameter scalarBlockParam,
                                                             final Parameter rhoParam,
                                                             final Parameter thetaParam,
                                                             final Parameter tParam,
                                                             final RhoOrdering rhoOrdering) {
        super(name, rotationParam, scalarBlockParam, rhoParam, thetaParam, tParam);
        if (!(rotationParam instanceof OrthogonalMatrixProvider)) {
            throw new IllegalArgumentException("rotationParam must implement OrthogonalMatrixProvider");
        }
        this.rhoOrdering = rhoOrdering == null ? RhoOrdering.NONE : rhoOrdering;
        this.nativeCompoundParameter = new CompoundParameter(getClass().getSimpleName() + ".native");
        this.nativeCompoundParameter.addParameter(scalarBlockParam);
        this.nativeCompoundParameter.addParameter(rhoParam);
        this.nativeCompoundParameter.addParameter(thetaParam);
        this.nativeCompoundParameter.addParameter(tParam);
        this.nativeCompoundParameter.addParameter(((OrthogonalMatrixProvider) rotationParam).getOrthogonalParameter());
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

        addRhoGradient(out, baseOffset, blockIndex,
                g00 * cos + g01 * sin + g10 * sin + g11 * cos);
        out[twoByTwoGradientOffset(baseOffset, 1, blockIndex)] =
                g00 * (-rho * sin) + g01 * (rho * cos) + g10 * (rho * cos) + g11 * (-rho * sin);
        out[twoByTwoGradientOffset(baseOffset, 2, blockIndex)] =
                -g01 + g10;
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

    @Override
    public CompoundParameter getParameter() {
        return nativeCompoundParameter;
    }

    public Parameter getRotationAngleParameter() {
        return ((OrthogonalMatrixProvider) getRotationMatrixParameter()).getOrthogonalParameter();
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
            final MatrixParameter rotation = (MatrixParameter) xo.getElementFirstChild(R_ELEMENT);
            final Parameter rhoParam = (Parameter) xo.getElementFirstChild(RHO_ELEMENT);
            final Parameter thetaParam = (Parameter) xo.getElementFirstChild(THETA_ELEMENT);
            final Parameter tParam = (Parameter) xo.getElementFirstChild(T_ELEMENT);
            final RhoOrdering rhoOrdering = RhoOrdering.parse(xo.getAttribute(RHO_ORDERING, RhoOrdering.NONE.name));

            final boolean oddDimension = (rotation.getRowDimension() & 1) == 1;
            final Parameter scalarBlockParam;
            if (xo.hasChildNamed(SCALAR_ELEMENT)) {
                scalarBlockParam = (Parameter) xo.getElementFirstChild(SCALAR_ELEMENT);
            } else if (oddDimension) {
                throw new XMLParseException("Missing required element '" + SCALAR_ELEMENT + "' for odd matrix dimension.");
            } else {
                scalarBlockParam = new Parameter.Default(0);
            }

            return new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                    name, rotation, scalarBlockParam, rhoParam, thetaParam, tParam, rhoOrdering
            );
        }

        @Override
        public String getParserDescription() {
            return "A matrix parameter A = R D R^T with orthogonal R and polar-stable 2x2 blocks.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return OrthogonalBlockDiagonalPolarStableMatrixParameter.class;
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
