package dr.evomodel.speciation.ModelAveragingResearch;


import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;


/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class ModelAveragingIndexSpeciationLikelihoodParser extends AbstractXMLObjectParser {
    public static final String MODEL_AVE_Index_SPECIATION_LIKELIHOOD = "modelAveragingIndexSpeciationLikelihood";
    public static final String INDEX = "modelIndex";
    public static final String MAX_INDEX = "maxIndex";

    public String getParserName() {
        return MODEL_AVE_Index_SPECIATION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Variable<Integer> index;

//        cxo = xo.getChild(INDEX);
        index = (Variable<Integer>) xo.getElementFirstChild(INDEX); // integer index parameter size = real size - 1
        Parameter maxIndex = (Parameter) xo.getElementFirstChild(MAX_INDEX);

        return new ModelAveragingIndexSpeciationLikelihood(xo.getId(), index, maxIndex);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Model Averaging Speciation Likelihood.";
    }

    public Class getReturnType() {
        return ModelAveragingIndexSpeciationLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(INDEX, new XMLSyntaxRule[]{
                    new ElementRule(Variable.class)
            }),
            new ElementRule(MAX_INDEX, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };

}
