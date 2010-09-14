package dr.evomodelxml.speciation;

import dr.evolution.util.Taxon;
import dr.evomodel.speciation.SpeciesBindings;
import dr.xml.*;

/**
 */
public class SpeciesBindingsSPinfoParser extends AbstractXMLObjectParser {
    public static final String SP = "sp";

    public String getParserName() {
        return SP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Taxon[] taxa = new Taxon[xo.getChildCount()];
        for (int nt = 0; nt < taxa.length; ++nt) {
            taxa[nt] = (Taxon) xo.getChild(nt);
        }
        return new SpeciesBindings.SPinfo(xo.getId(), taxa);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Taxon.class, 1, Integer.MAX_VALUE)
        };
    }

    public String getParserDescription() {
        return "Taxon in a species tree";
    }

    public Class getReturnType() {
        return SpeciesBindings.SPinfo.class;
    }
}
