package dr.evomodelxml.antigenic;

import dr.evomodel.antigenic.AntigenicGradientWrtParameter;
import dr.evomodel.antigenic.AntigenicLikelihoodGradient;
import dr.evomodel.antigenic.NewAntigenicLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class AntigenicLikelihoodGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "antigenicGradient";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        NewAntigenicLikelihood likelihood = (NewAntigenicLikelihood) xo.getChild(NewAntigenicLikelihood.class);
        List<Parameter> parameters = xo.getAllChildren(Parameter.class);

        List<AntigenicGradientWrtParameter> wrtList = new ArrayList<>();
        for (Parameter parameter : parameters) {
            AntigenicGradientWrtParameter wrt = likelihood.wrtFactory(parameter);
            wrtList.add(wrt);
        }

        return new AntigenicLikelihoodGradient(likelihood, wrtList);
    }

    public String getParserDescription() {
        return "Provides the gradient of the likelihood of immunological assay data such as " +
                "hemagglutinin inhibition (HI) given vectors of coordinates " +
                "for viruses and sera/antisera in some multidimensional 'antigenic' space.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(NewAntigenicLikelihood.class),
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
    };

    public Class getReturnType() {
        return AntigenicLikelihoodGradient.class;
    }
}

