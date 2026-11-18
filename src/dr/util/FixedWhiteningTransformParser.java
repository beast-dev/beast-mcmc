/*
 * FixedWhiteningTransformParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 */

package dr.util;

import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/** XML parser for {@link FixedWhiteningTransform}. */
public final class FixedWhiteningTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "fixedWhiteningTransform";
    private static final String PARAMETER = "parameter";
    private static final String INNER = "inner";
    private static final String REFERENCE = "reference";
    private static final String INVERSE_MAP = "inverseMap";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Parameter parameter = (Parameter) xo.getElementFirstChild(PARAMETER);
        final Transform.MultivariableTransform nativeTransform =
                (Transform.MultivariableTransform) xo.getElementFirstChild(INNER);
        final Parameter referenceParameter = (Parameter) xo.getElementFirstChild(REFERENCE);
        final Parameter mapParameter = (Parameter) xo.getElementFirstChild(INVERSE_MAP);
        final int dimension = parameter.getDimension();

        if (referenceParameter.getDimension() != dimension) {
            throw new XMLParseException("The fixed-whitening reference must have dimension " + dimension + ".");
        }
        if (mapParameter.getDimension() != dimension * dimension) {
            throw new XMLParseException("The fixed-whitening inverse map must contain " +
                    (dimension * dimension) + " row-major entries.");
        }

        final double[][] map = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                map[i][j] = mapParameter.getParameterValue(i * dimension + j);
            }
        }
        try {
            return new FixedWhiteningTransform(parameter, nativeTransform,
                    referenceParameter.getParameterValues(), map);
        } catch (IllegalArgumentException exception) {
            throw new XMLParseException(exception.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(PARAMETER, Parameter.class),
                new ElementRule(INNER, Transform.MultivariableTransform.class),
                new ElementRule(REFERENCE, Parameter.class),
                new ElementRule(INVERSE_MAP, Parameter.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Applies a fixed dense whitening after a native multivariable transform.";
    }

    @Override
    public Class getReturnType() {
        return FixedWhiteningTransform.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
