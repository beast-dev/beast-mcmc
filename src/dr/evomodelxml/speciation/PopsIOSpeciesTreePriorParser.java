package dr.evomodelxml.speciation;

import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.evomodel.speciation.PopsIOSpeciesTreePrior;
import dr.evomodel.speciation.SpeciationModel;
import dr.xml.*;

/**
 * User: Graham Jones
 * Date: 10/05/12
 */
public class PopsIOSpeciesTreePriorParser extends AbstractXMLObjectParser {
    public static final String POPSIO_SPECIES_TREE_PRIOR = "PopsIOSpeciesTreePrior";
    public static final String MODEL = "model";
    public static final String PIO_TREE = "pioTree";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        System.out.println("PopsIOSpeciesTreePriorParser");
        final XMLObject mxo = xo.getChild(MODEL);
        final SpeciationModel sppm = (SpeciationModel) mxo.getChild(SpeciationModel.class);
        final XMLObject mulsptxo = xo.getChild(PIO_TREE);
        final PopsIOSpeciesTreeModel piostm = (PopsIOSpeciesTreeModel) mulsptxo.getChild(PopsIOSpeciesTreeModel.class);
        return new PopsIOSpeciesTreePrior(sppm, piostm);
    }


    private  XMLSyntaxRule[] speciationModelSyntax() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciationModel.class)
        };

    }

    private  XMLSyntaxRule[] piostmSyntax() {
        return new XMLSyntaxRule[]{
                new ElementRule(PopsIOSpeciesTreeModel.class)
        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MODEL, speciationModelSyntax()),
                new ElementRule(PIO_TREE, piostmSyntax()),

        };
    }

    @Override
    public String getParserDescription() {
        return "Prior for a species tree.";
    }

    @Override
    public Class getReturnType() {
        return PopsIOSpeciesTreePrior.class;
    }

    public String getParserName() {
        return POPSIO_SPECIES_TREE_PRIOR;
    }
}
