package dr.evomodelxml.sitemodel;

import dr.evomodel.sitemodel.CategorySiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class CategorySiteModelParser extends AbstractXMLObjectParser {

    public static final String SITE_MODEL = "categorySiteModel";
    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String RATE_PARAMETER = "rates";
    public static final String CATEGORIES = "categories";
    public static final String CATEGORY_STATES = "states";
    public static final String CATEGORY_STRING = "values";
    public static final String RELATIVE_TO = "relativeTo";

    public String getParserName() {
        return SITE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(SUBSTITUTION_MODEL);
        SubstitutionModel substitutionModel = (SubstitutionModel) cxo.getChild(SubstitutionModel.class);

        cxo = xo.getChild(MUTATION_RATE);
        Parameter muParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RATE_PARAMETER);
        Parameter rateParam = null;
        int relativeTo = 0;
        if (cxo != null) {
            rateParam = (Parameter) cxo.getChild(Parameter.class);
            relativeTo = cxo.getIntegerAttribute(RELATIVE_TO);
        }

        cxo = xo.getChild(CATEGORIES);
        String categories = "";
        String states = "";
        if (cxo != null) {
            categories = cxo.getStringAttribute(CATEGORY_STRING);
            states = cxo.getStringAttribute(CATEGORY_STATES);
        }

        return new CategorySiteModel(substitutionModel, muParam, rateParam, categories, states, relativeTo);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteModel that has a gamma distributed rates across sites";
    }

    public Class getReturnType() {
        return CategorySiteModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                    new ElementRule(SubstitutionModel.class)
            }),
            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(RATE_PARAMETER, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(RELATIVE_TO, true),
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(CATEGORIES, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(CATEGORY_STRING, true),
                    AttributeRule.newStringRule(CATEGORY_STATES, true)
            }, true),
    };

}
