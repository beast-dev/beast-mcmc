package dr.evomodelxml.operators;

import dr.evomodel.operators.AlloppChangeNumHybridizations;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 22/07/12
 * Time: 11:24
 */
public class AlloppChangeNumHybridizationsParser extends AbstractXMLObjectParser {

    public static final String CHANGE_NUM_HYBRIDIZATIONS = "changeNumHybridizations";


    public String getParserName() {
        return CHANGE_NUM_HYBRIDIZATIONS;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
        AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new AlloppChangeNumHybridizations(apspnet, apsp, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator which changes the number of tetraploid subtrees (hybridizations) by merging and splitting them.";

    }

    @Override
    public Class getReturnType() {
        return AlloppChangeNumHybridizations.class;
    }



}
