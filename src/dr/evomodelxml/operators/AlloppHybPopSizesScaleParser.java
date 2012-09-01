package dr.evomodelxml.operators;

import dr.evomodel.operators.AlloppHybPopSizesScale;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.inference.operators.MCMCOperator;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 03/08/12
 */
public class AlloppHybPopSizesScaleParser extends AbstractXMLObjectParser {

    public static final String HYB_POP_SIZES_SCALE = "hybPopSizesScaleOperator";


    public String getParserName() {
        return HYB_POP_SIZES_SCALE;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
        AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scalingfactor = xo.getDoubleAttribute(ScaleOperatorParser.SCALE_FACTOR);
        return new AlloppHybPopSizesScale(apspnet, apsp, scalingfactor, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(ScaleOperatorParser.SCALE_FACTOR),
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator which scales the population size of a newly formed hybrid.";

    }

    @Override
    public Class getReturnType() {
        return AlloppHybPopSizesScale.class;
    }

}
