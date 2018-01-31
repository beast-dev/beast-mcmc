package dr.inferencexml.distribution;

import dr.inference.distribution.IndependentNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class IndependentNormalDistributionModelParser extends AbstractXMLObjectParser {
    public static String INDEPENDENT_NORMAL_DISTRIBUTION_MODEL = "independentNormalDistributionModel";
    public static String MEAN = "mean";
    public static String VARIANCE = "variance";
    public static String PRECISION = "precision";
    public static String DATA = "data";
    public static String ID = "id";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String id = xo.getStringAttribute(ID);
        Parameter mean = (Parameter) xo.getChild(MEAN).getChild(Parameter.class);
        Parameter precision = null;
        if(xo.getChild(PRECISION) != null){
            precision = (Parameter) xo.getChild(PRECISION).getChild(Parameter.class);
        }
        Parameter variance = null;
        if(xo.getChild(VARIANCE) != null){
            variance = (Parameter) xo.getChild(VARIANCE).getChild(Parameter.class);
        }
        Parameter data = (Parameter) xo.getChild(DATA).getChild(Parameter.class);

        return new IndependentNormalDistributionModel(id, mean, variance, precision, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(MEAN, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(VARIANCE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(PRECISION, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            ),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            })
    };

    @Override
    public String getParserDescription() {
        return "A series of independent normal distribution models";
    }

    @Override
    public Class getReturnType() {
        return IndependentNormalDistributionModel.class;
    }

    @Override
    public String getParserName() {
        return INDEPENDENT_NORMAL_DISTRIBUTION_MODEL;
    }
}
