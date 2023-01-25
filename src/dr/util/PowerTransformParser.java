package dr.util;


import dr.xml.*;

/**
 * @author Gabriel Hassler
 * @author Marc Suchard
 */

public class PowerTransformParser extends AbstractXMLObjectParser {
    public static final String POWER_TRANSFORM = "powerTransform";
    private static final String POWER = "power";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Transform.ParsedTransform parsedTransform = (Transform.ParsedTransform) xo.getChild(Transform.ParsedTransform.class);
        if (!(parsedTransform.transform instanceof Transform.PowerTransform)) {
            throw new XMLParseException("The '" + TransformParsers.TYPE + "' attribute of the " +
                    TransformParsers.TRANSFORM + " xml element must be '" + Transform.Type.POWER.getName() + "'.");
        }

        double power = xo.getDoubleAttribute(POWER);

        parsedTransform.transform = new Transform.PowerTransform(power);
        return parsedTransform;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(POWER),
                new ElementRule(Transform.ParsedTransform.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Raises a parameter to a given power.";
    }

    @Override
    public Class getReturnType() {
        return Transform.ParsedTransform.class;
    }

    @Override
    public String getParserName() {
        return POWER_TRANSFORM;
    }
}
