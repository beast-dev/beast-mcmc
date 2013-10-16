package dr.evomodelxml.speciation;

import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * User: Graham Jones
 * Date: 10/05/12
 */


public class PopsIOSpeciesBindingsParser extends AbstractXMLObjectParser {
    public static final String PIO_SPECIES_BINDINGS = "pioSpeciesBindings";
    public static final String GENE_TREES = "geneTrees";
    public static final String GTREE = "gtree";
    public static final String POPFACTOR = "popfactor";
    public static final String MIN_GENENODE_HEIGHT = "minGeneNodeHeight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        System.out.println("PopsIOSpeciesBindingsParser");
        List<PopsIOSpeciesBindings.SpInfo> piosb = new ArrayList<PopsIOSpeciesBindings.SpInfo>();
        for (int k = 0; k < xo.getChildCount(); ++k) {
            final Object child = xo.getChild(k);
            if (child instanceof PopsIOSpeciesBindings.SpInfo) {
                piosb.add((PopsIOSpeciesBindings.SpInfo) child);
            }
        }
        final double mingenenodeheight = xo.getDoubleAttribute(MIN_GENENODE_HEIGHT);
        final XMLObject xogt = xo.getChild(GENE_TREES);
        final int nTrees = xogt.getChildCount();
        final TreeModel[] trees = new TreeModel[nTrees];
        double[] popFactors = new double[nTrees];

        for (int nt = 0; nt < trees.length; ++nt) {
            Object child = xogt.getChild(nt);
            if (!(child instanceof TreeModel)) {
                assert child instanceof XMLObject;
                popFactors[nt] = ((XMLObject) child).getDoubleAttribute(POPFACTOR);
                child = ((XMLObject) child).getChild(TreeModel.class);

            } else {
                popFactors[nt] = -1;
            }
            trees[nt] = (TreeModel) child;
        }

        try {
            return new PopsIOSpeciesBindings(piosb.toArray(new PopsIOSpeciesBindings.SpInfo[piosb.size()]),
                    trees, mingenenodeheight, popFactors);
        } catch (Error e) {
            throw new XMLParseException(e.getMessage());
        }
    }


    private XMLSyntaxRule[] treeWithPopFactorsSyntax() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(POPFACTOR),
                new ElementRule(TreeModel.class)
        };
    }

    private ElementRule treeWithPopFactors() {
        return new ElementRule(GTREE, treeWithPopFactorsSyntax(), 0, Integer.MAX_VALUE);
    }

    private XMLSyntaxRule[] geneTreesSyntax() {
        return new XMLSyntaxRule[]{
                 new ElementRule(TreeModel.class, 0, Integer.MAX_VALUE),
                 treeWithPopFactors()
        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MIN_GENENODE_HEIGHT),
                new ElementRule(PopsIOSpeciesBindings.SpInfo.class, 2, Integer.MAX_VALUE),
                new ElementRule(GENE_TREES, geneTreesSyntax())
        };
    }

    @Override
    public String getParserDescription() {
        return "Binds taxa in gene trees with species information.";
    }

    @Override
    public Class getReturnType() {
        return PopsIOSpeciesBindings.class;
    }

    public String getParserName() {
        return PIO_SPECIES_BINDINGS;
    }
}
