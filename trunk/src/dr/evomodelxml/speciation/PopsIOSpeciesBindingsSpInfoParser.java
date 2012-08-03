package dr.evomodelxml.speciation;

import dr.evolution.util.Taxon;
import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.xml.*;

/*
 * User: Graham Jones
 * Date: 11/05/12
 */

public class PopsIOSpeciesBindingsSpInfoParser  extends AbstractXMLObjectParser {
    public static final String PIOSP = "pioSp";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Taxon[] taxa = new Taxon[xo.getChildCount()];
        for (int nt = 0; nt < taxa.length; ++nt) {
            taxa[nt] = (Taxon) xo.getChild(nt);
        }
        return new PopsIOSpeciesBindings.SpInfo(xo.getId(), taxa);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Taxon.class, 1, Integer.MAX_VALUE)
        };
    }

    @Override
    public String getParserDescription() {
        return "A species made of taxa (sequences)";
    }

    @Override
    public Class getReturnType() {
        return PopsIOSpeciesBindings.SpInfo.class;
    }

    public String getParserName() {
        return PIOSP;
    }
}
