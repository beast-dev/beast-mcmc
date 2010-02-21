package dr.evomodelxml.speciation;

import dr.evomodel.speciation.MultiSpeciesCoalescent;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.xml.*;

/**
 */
public class MultiSpeciesCoalescentParser extends AbstractXMLObjectParser {
    public static final String SPECIES_COALESCENT = "speciesCoalescent";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final SpeciesBindings sb = (SpeciesBindings) xo.getChild(SpeciesBindings.class);
        final SpeciesTreeModel tree = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);
        return new MultiSpeciesCoalescent(sb, tree);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciesBindings.class),
                new ElementRule(SpeciesTreeModel.class),
        };
    }

    public String getParserDescription() {
        return "Compute coalecent log-liklihood of a set of gene trees embedded inside one species tree.";
    }

    public Class getReturnType() {
        return MultiSpeciesCoalescent.class;
    }

    public String getParserName() {
        return SPECIES_COALESCENT;
    }
}
