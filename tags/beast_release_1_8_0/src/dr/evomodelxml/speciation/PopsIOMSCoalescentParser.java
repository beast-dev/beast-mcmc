package dr.evomodelxml.speciation;

import dr.evomodel.speciation.PopsIOMSCoalescent;
import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.xml.*;

/**
 * User: Graham Jones
 * Date: 10/05/12
 */
public class PopsIOMSCoalescentParser extends AbstractXMLObjectParser {

    public static final String POPSIO_MSCOALESCENT = "PopsIOMSCoalescent";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        System.out.println("PopsIOMSCoalescentParser");
        final PopsIOSpeciesBindings piosb = (PopsIOSpeciesBindings) xo.getChild(PopsIOSpeciesBindings.class);
        final PopsIOSpeciesTreeModel piostm =
                (PopsIOSpeciesTreeModel) xo.getChild(PopsIOSpeciesTreeModel.class);
        return new PopsIOMSCoalescent(piosb, piostm);
    }



    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(PopsIOSpeciesBindings.class),
                new ElementRule(PopsIOSpeciesTreeModel.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Likelihood of a set of gene trees embedded in a species tree.";
    }

    @Override
    public Class getReturnType() {
        return PopsIOMSCoalescent.class;
    }

    public String getParserName() {
        return POPSIO_MSCOALESCENT;
    }
}
