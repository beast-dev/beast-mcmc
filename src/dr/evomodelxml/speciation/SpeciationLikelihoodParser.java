package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.inference.model.Statistic;
import dr.math.distributions.Distribution;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
   // public static final String COEFFS = "coefficients";

    public static final String CALIBRATION = "calibration";

    public String getParserName() {
        return SPECIATION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

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

        final XMLObject cal = xo.getChild(CALIBRATION);
        if( cal != null ) {
            if( excludeTaxa != null ) {
                throw new XMLParseException("Sorry, not implemented: internal calibration prior + excluded taxa");
            }

            if( ! specModel.supportsInternalCalibration() ) {
              throw new XMLParseException("Sorry, not implemented: internal calibration prior for this model.");
            }

            //final double[] coef = cal.hasAttribute(COEFFS) ? cal.getDoubleArrayAttribute(COEFFS) : null;

            List<Distribution> dists = new ArrayList<Distribution>();
            List<Taxa> taxa = new ArrayList<Taxa>();
            for(int k = 0; k < cal.getChildCount(); ++k) {
                final Object ck = cal.getChild(k);
                if ( Distribution.class.isInstance(ck) ) {
                    dists.add((Distribution) ck);
                }
                if ( Taxa.class.isInstance(ck) ) {
                    taxa.add((Taxa) ck);
                }
            }

            if( dists.size() != taxa.size() ) {
                throw new XMLParseException("Mismatch in number of distributions and taxa specs");
            }

            final Statistic s = (Statistic) cal.getChild(Statistic.class);
            if( dists.size() > 1 && s == null ) {
                throw new XMLParseException("Sorry, not implemented: multiple internal calibrations - please provide the " +
                        "correction explicitly.");
            }

            return new SpeciationLikelihood(tree, specModel, null, dists, taxa, s);
//                    (Distribution) cal.getChild(Distribution.class),
//                    (Taxa) cal.getChild(Taxa.class),
//                    coef);
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


    private final XMLSyntaxRule[] calibration = {
//            AttributeRule.newDoubleArrayRule(COEFFS,true, "use log(lam) -lam * c[0] + sum_k=1..n (c[k+1] * e**(-k*lam*x)) " +
//                    "as a calibration correction instead of default - used when additional constarints are put on the topology."),
            new ElementRule(Statistic.class, true),
            new ElementRule(Distribution.class, 1, 100),
            new ElementRule(Taxa.class,1, 100)
    };

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
            }, "One or more subsets of taxa which should be excluded from calculate the likelihood (which is calculated on the remaining subtree)", true),

            new ElementRule(CALIBRATION, calibration, true),
    };

}
