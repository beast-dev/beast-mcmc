package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 */
public class SpeciationLikelihoodParser extends AbstractXMLObjectParser {

    public static final String SPECIATION_LIKELIHOOD = "speciationLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "speciesTree";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public String getParserName() {
        return SPECIATION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) {

        XMLObject cxo = xo.getChild(MODEL);
        SpeciationModel specModel = (SpeciationModel) cxo.getChild(SpeciationModel.class);

        cxo = xo.getChild(TREE);
        Tree tree = (Tree) cxo.getChild(Tree.class);

        Set<Taxon> excludeTaxa = null;

        if (xo.hasChildNamed(INCLUDE)) {
            excludeTaxa = new HashSet<Taxon>();
            for (int i = 0; i < tree.getTaxonCount(); i++) {
                excludeTaxa.add(tree.getTaxon(i));
            }

            cxo = xo.getChild(INCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                TaxonList taxonList = (TaxonList) cxo.getChild(i);
                for (int j = 0; j < taxonList.getTaxonCount(); j++) {
                    excludeTaxa.remove(taxonList.getTaxon(j));
                }
            }
        }

        if (xo.hasChildNamed(EXCLUDE)) {
            excludeTaxa = new HashSet<Taxon>();
            cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                TaxonList taxonList = (TaxonList) cxo.getChild(i);
                for (int j = 0; j < taxonList.getTaxonCount(); j++) {
                    excludeTaxa.add(taxonList.getTaxon(j));
                }
            }
        }
        if (excludeTaxa != null) {
            Logger.getLogger("dr.evomodel").info("Speciation model excluding " + excludeTaxa.size() + " taxa from prior - " +
                    (tree.getTaxonCount() - excludeTaxa.size()) + " taxa remaining.");
        }

        return new SpeciationLikelihood(tree, specModel, excludeTaxa, null);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the speciation.";
    }

    public Class getReturnType() {
        return SpeciationLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(SpeciationModel.class)
            }),
            new ElementRule(TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),

            new ElementRule(INCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be included from calculate the likelihood (the remaining taxa are excluded)", true),

            new ElementRule(EXCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be excluded from calculate the likelihood (which is calculated on the remaining subtree)", true)
    };

}
