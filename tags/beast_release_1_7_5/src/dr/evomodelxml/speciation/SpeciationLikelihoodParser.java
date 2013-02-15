package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.speciation.CalibrationPoints;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.inference.distribution.DistributionLikelihood;
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

    public static final String CALIBRATION = "calibration";
    public static final String CORRECTION = "correction";
    public static final String POINT = "point";

    private final String EXACT = CalibrationPoints.CorrectionType.EXACT.toString();
    private final String APPROX = CalibrationPoints.CorrectionType.APPROXIMATED.toString();
    private final String PEXACT = CalibrationPoints.CorrectionType.PEXACT.toString();
    private final String NONE = CalibrationPoints.CorrectionType.NONE.toString();

    public static final String PARENT = dr.evomodelxml.tree.TMRCAStatisticParser.PARENT;

    public String getParserName() {
        return SPECIATION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MODEL);
        final SpeciationModel specModel = (SpeciationModel) cxo.getChild(SpeciationModel.class);

        cxo = xo.getChild(TREE);
        final Tree tree = (Tree) cxo.getChild(Tree.class);

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

            List<Distribution> dists = new ArrayList<Distribution>();
            List<Taxa> taxa = new ArrayList<Taxa>();
            List<Boolean> forParent = new ArrayList<Boolean>();
            Statistic userPDF = null; // (Statistic) cal.getChild(Statistic.class);

            for(int k = 0; k < cal.getChildCount(); ++k) {
                final Object ck = cal.getChild(k);
                if ( DistributionLikelihood.class.isInstance(ck) ) {
                    dists.add( ((DistributionLikelihood) ck).getDistribution() );
                } else if ( Distribution.class.isInstance(ck) ) {
                    dists.add((Distribution) ck);
                } else if ( Taxa.class.isInstance(ck) ) {
                    final Taxa tx = (Taxa) ck;
                    taxa.add(tx);
                    forParent.add( tx.getTaxonCount() == 1 );
                } else if ( Statistic.class.isInstance(ck) ) {
                    if( userPDF != null ) {
                        throw new XMLParseException("more than one userPDF correction???");
                    }
                    userPDF = (Statistic) cal.getChild(Statistic.class);
                }
                else {
                    XMLObject cko = (XMLObject) ck;
                    assert cko.getChildCount() == 2;

                    for(int i = 0; i < 2; ++i) {
                        final Object chi = cko.getChild(i);
                        if ( DistributionLikelihood.class.isInstance(chi) ) {
                            dists.add( ((DistributionLikelihood) chi).getDistribution() );
                        } else if ( Distribution.class.isInstance(chi) ) {
                            dists.add((Distribution) chi);
                        } else if ( Taxa.class.isInstance(chi) ) {
                            taxa.add((Taxa) chi);
                            boolean fp = ((Taxa) chi).getTaxonCount() == 1;
                            if( cko.hasAttribute(PARENT) ) {
                                boolean ufp = cko.getBooleanAttribute(PARENT);
                                if( fp && ! ufp ) {
                                   throw new XMLParseException("forParent==false for a single taxon?? (must be true)");
                                }
                                fp = ufp;
                            }
                            forParent.add(fp);
                        } else {
                            assert false;
                        }
                    }
                }
            }

            if( dists.size() != taxa.size() ) {
                throw new XMLParseException("Mismatch in number of distributions and taxa specs");
            }

            try {
                final String correction = cal.getAttribute(CORRECTION, EXACT);

                final CalibrationPoints.CorrectionType type = correction.equals(EXACT) ? CalibrationPoints.CorrectionType.EXACT :
                        (correction.equals(APPROX) ? CalibrationPoints.CorrectionType.APPROXIMATED :
                                (correction.equals(NONE) ? CalibrationPoints.CorrectionType.NONE :
                                        (correction.equals(PEXACT) ? CalibrationPoints.CorrectionType.PEXACT :  null)));

                if( cal.hasAttribute(CORRECTION) && type == null ) {
                   throw new XMLParseException("correction type == " + correction + "???");
                }

                final CalibrationPoints calib =
                        new CalibrationPoints(tree, specModel.isYule(), dists, taxa, forParent, userPDF, type);
                final SpeciationLikelihood speciationLikelihood = new SpeciationLikelihood(tree, specModel, null, calib);
                return speciationLikelihood;
            } catch( IllegalArgumentException e ) {
                throw new XMLParseException( e.getMessage() );
            }
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

    private final XMLSyntaxRule[] calibrationPoint = {
            AttributeRule.newBooleanRule(PARENT, true),
            new XORRule(
                    new ElementRule(Distribution.class),
                    new ElementRule(DistributionLikelihood.class)),
            new ElementRule(Taxa.class)
    };

    private final XMLSyntaxRule[] calibration = {
//            AttributeRule.newDoubleArrayRule(COEFFS,true, "use log(lam) -lam * c[0] + sum_k=1..n (c[k+1] * e**(-k*lam*x)) " +
//                    "as a calibration correction instead of default - used when additional constarints are put on the topology."),
            AttributeRule.newStringRule(CORRECTION, true),
            new ElementRule(Statistic.class, true),
            new XORRule(
                    new ElementRule(Distribution.class, 1, 100),
                    new ElementRule(DistributionLikelihood.class, 1, 100)),
            new ElementRule(Taxa.class, 1, 100),
            new ElementRule("point", calibrationPoint, 0, 100)
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
