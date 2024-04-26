package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.InfinitesimalRatesLogger;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

public class InfinitesimalRatesLoggerParser extends AbstractXMLObjectParser {

    private static final String NAME = "infinitesimalRatesLogger";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        return new InfinitesimalRatesLogger(substitutionModel);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(SubstitutionModel.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Logger to report infinitesimal rates of a substitution model";
    }

    @Override
    public Class getReturnType() {
        return InfinitesimalRatesLogger.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
