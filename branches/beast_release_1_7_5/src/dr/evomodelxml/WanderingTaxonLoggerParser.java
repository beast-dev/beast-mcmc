package dr.evomodelxml;

import dr.xml.*;
import dr.evomodel.tree.WanderingTaxonLogger;
import dr.evolution.util.Taxon;

/**
 * @author Marc A. Suchard
 */
public class WanderingTaxonLoggerParser extends AbstractXMLObjectParser {

    public static final String NAME = "name";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(NAME, xo.getId());

        WanderingTaxonLogger.Relative relative = WanderingTaxonLogger.Relative.SISTER;
        if (xo.getAttribute(WanderingTaxonLogger.RELATIVE,"sister").equalsIgnoreCase("parent")) {
            relative = WanderingTaxonLogger.Relative.PARENT;
        }

        Taxon taxon = (Taxon) xo.getChild(Taxon.class);

        return new WanderingTaxonLogger(name,taxon,relative);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newStringRule(NAME,true),
                AttributeRule.newStringRule(WanderingTaxonLogger.RELATIVE,true),
                new ElementRule(Taxon.class),
        };
    }

    public String getParserDescription() {
        return null;
    }

    public Class getReturnType() {
        return WanderingTaxonLogger.class;
    }

    public String getParserName() {
        return WanderingTaxonLogger.WANDERER;
    }
}
