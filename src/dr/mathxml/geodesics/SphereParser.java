package dr.mathxml.geodesics;

import dr.math.geodesics.Sphere;
import dr.xml.*;

public class SphereParser extends AbstractXMLObjectParser {

    private static final String SPHERE = "sphere";
    private static final String RADIUS = "radius";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        if (xo.hasAttribute(RADIUS)) {
            double radius = xo.getDoubleAttribute(RADIUS);
            return new Sphere(radius);
        }
        return new Sphere();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(RADIUS, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Spherical manifold";
    }

    @Override
    public Class getReturnType() {
        return Sphere.class;
    }

    @Override
    public String getParserName() {
        return SPHERE;
    }
}
