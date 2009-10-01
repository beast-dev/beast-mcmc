package dr.evomodelxml;

import dr.xml.*;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.treelikelihood.MicrosatelliteSamplerTreeLikelihood;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MicrosatelliteSamplerTreeLikelihood
 */
public class MicrosatelliteSamplerTreeLikelihoodParser extends AbstractXMLObjectParser {
    public static final String MUTATION_RATE = "mutationRate";
    public String getParserName(){
        return "microsatelliteSamplerTreeLikelihood";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MicrosatelliteSamplerTreeModel mstModel = (MicrosatelliteSamplerTreeModel)xo.getChild(MicrosatelliteSamplerTreeModel.class);
        MicrosatelliteModel microsatelliteModel = (MicrosatelliteModel)xo.getChild(MicrosatelliteModel.class);
        Parameter muRate = new Parameter.Default(1.0);
        if(xo.hasChildNamed(MUTATION_RATE)){
            muRate = (Parameter)xo.getElementFirstChild(MUTATION_RATE);
        }
        return new MicrosatelliteSamplerTreeLikelihood(mstModel,microsatelliteModel, muRate);
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MicrosatelliteSamplerTreeModel.class),
            new ElementRule(MicrosatelliteModel.class),
            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true)
    };

    public String getParserDescription(){
        return "this parser returns an object of the TreeMicrosatelliteSamplerLikelihood class";
    }

    public Class getReturnType(){
        return MicrosatelliteSamplerTreeLikelihood.class;
    }


}
