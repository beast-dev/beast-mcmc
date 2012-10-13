package dr.evomodelxml.coalescent;

import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class NewCoalescentSimulatorParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_SIMULATOR = "coalescentSimulator";
    public static final String HEIGHT = "height";

    public String getParserName() {
        return COALESCENT_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoalescentSimulator simulator = new CoalescentSimulator();

        DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);
        List<TaxonList> taxonLists = new ArrayList<TaxonList>();
        List<Tree> subtrees = new ArrayList<Tree>();

        double height = xo.getAttribute(HEIGHT, Double.NaN);

        // should have one child that is node
        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);

            // AER - swapped the order of these round because Trees are TaxonLists...
            if (child instanceof Tree) {
                subtrees.add((Tree) child);
            } else if (child instanceof TaxonList) {
                taxonLists.add((TaxonList) child);
            }
          }

        if (taxonLists.size() == 0) {
            if (subtrees.size() == 1) {
                return subtrees.get(0);
            }
            throw new XMLParseException("Expected at least one taxonList or two subtrees in "
                    + getParserName() + " element.");
        }

        Taxa remainingTaxa = new Taxa();
        for (int i = 0; i < taxonLists.size(); i++) {
            remainingTaxa.addTaxa(taxonLists.get(i));
        }

        for (int i = 0; i < subtrees.size(); i++) {
            remainingTaxa.removeTaxa(subtrees.get(i));
        }

        try {
            Tree[] trees = new Tree[subtrees.size() + remainingTaxa.getTaxonCount()];
            // add the preset trees
            for (int i = 0; i < subtrees.size(); i++) {
                trees[i] = subtrees.get(i);
            }

            // add all the remaining taxa in as single tip trees...
            for (int i = 0; i < remainingTaxa.getTaxonCount(); i++) {
                Taxa tip = new Taxa();
                tip.addTaxon(remainingTaxa.getTaxon(i));
                trees[i + subtrees.size()] = simulator.simulateTree(tip, demoModel);
            }

            return simulator.simulateTree(trees, demoModel, height, trees.length != 1);

        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a simulated tree under the given demographic model. The element can " +
                "be nested to simulate with monophyletic clades. The tree will be rescaled to the given height.";
    }

    public Class getReturnType() {
        return CoalescentSimulator.class; //Object.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(HEIGHT, true, ""),
            new ElementRule(Tree.class, 0, Integer.MAX_VALUE),
            new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
            new ElementRule(DemographicModel.class, 0, Integer.MAX_VALUE),
    };
}
