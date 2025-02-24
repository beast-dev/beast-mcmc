package dr.mathxml.geodesics;

import dr.math.geodesics.Euclidean;
import dr.xml.*;

public class EuclideanParser extends AbstractXMLObjectParser {

    private static final String EUCLIDEAN = "euclidean";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new Euclidean();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{};
    }

    @Override
    public String getParserDescription() {
        return "Spherical manifold";
    }

    @Override
    public Class getReturnType() {
        return Euclidean.class;
    }

    @Override
    public String getParserName() {
        return EUCLIDEAN;
    }
}
