
package dr.evomodelxml.speciation;

import dr.evomodel.speciation.MulMSCoalescent;
import dr.evomodel.speciation.MulSpeciesBindings;
import dr.evomodel.speciation.MulSpeciesTreeModel;
import dr.xml.*;

/**
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */
public class MulMSCoalescentParser extends AbstractXMLObjectParser {
    public static final String MUL_MS_COALESCENT = "mulMSCoalescent";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final MulSpeciesBindings sb = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);
        final MulSpeciesTreeModel tree = (MulSpeciesTreeModel) xo.getChild(MulSpeciesTreeModel.class);
        return new MulMSCoalescent(sb, tree);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MulSpeciesBindings.class),
                new ElementRule(MulSpeciesTreeModel.class),
        };
    }

    public String getParserDescription() {
        return "Compute coalecent log-liklihood of a set of gene trees embedded inside one species tree.";
    }

    public Class getReturnType() {
        return MulMSCoalescent.class;
    }

    public String getParserName() {
        return MUL_MS_COALESCENT;
    }
}
