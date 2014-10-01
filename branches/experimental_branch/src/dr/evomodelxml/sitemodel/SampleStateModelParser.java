package dr.evomodelxml.sitemodel;

import dr.evomodel.sitemodel.SampleStateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Vector;

/**
 */
public class SampleStateModelParser extends AbstractXMLObjectParser {

    public static final String SAMPLE_STATE_MODEL = "sampleStateModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String PROPORTIONS = "proportions";

    public String getParserName() {
        return SAMPLE_STATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MUTATION_RATE);
        Parameter muParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PROPORTIONS);
        Parameter proportionParameter = (Parameter) cxo.getChild(Parameter.class);

        Vector<Object> subModels = new Vector<Object>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof SubstitutionModel) {
                subModels.addElement(xo.getChild(i));
            }
        }

        return new SampleStateModel(muParam, proportionParameter, subModels);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteModel that has a discrete distribution of substitution models over sites, " +
                "designed for sampling of internal states.";
    }

    public Class getReturnType() {
        return SampleStateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MUTATION_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PROPORTIONS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE)
    };
}
