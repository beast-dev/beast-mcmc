package dr.evomodelxml.coalescent;

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
public class CoalescentSimulatorParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_TREE = "coalescentTree";
    public static final String COALESCENT_SIMULATOR = "coalescentSimulator";
    public static final String RESCALE_HEIGHT = "rescaleHeight";
    public static final String ROOT_HEIGHT = TreeModelParser.ROOT_HEIGHT;
    public static final String CONSTRAINED_TAXA = "constrainedTaxa";
    public static final String TMRCA_CONSTRAINT = "tmrca";
    public static final String IS_MONOPHYLETIC = "monophyletic";

    public String getParserName() {
        return COALESCENT_TREE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoalescentSimulator simulator = new CoalescentSimulator();

        DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);
        List<TaxonList> taxonLists = new ArrayList<TaxonList>();
        List<Tree> subtrees = new ArrayList<Tree>();

        double rootHeight = xo.getAttribute(ROOT_HEIGHT, -1.0);

        if (xo.hasAttribute(RESCALE_HEIGHT)) {
            rootHeight = xo.getDoubleAttribute(RESCALE_HEIGHT);
        }

        // should have one child that is node
        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);

            // AER - swapped the order of these round because Trees are TaxonLists...
            if (child instanceof Tree) {
                subtrees.add((Tree) child);
            } else if (child instanceof TaxonList) {
                taxonLists.add((TaxonList) child);
            } else if (xo.getChildName(i).equals(CONSTRAINED_TAXA)) {
                XMLObject constrainedTaxa = (XMLObject) child;

                // all taxa
                final TaxonList taxa = (TaxonList) constrainedTaxa.getChild(TaxonList.class);

                List<CoalescentSimulator.TaxaConstraint> constraints = new ArrayList<CoalescentSimulator.TaxaConstraint>();
                final String setsNotCompatibleMessage = "taxa sets not compatible";

                // pick up all constraints. order in partial order, where taxa_1 @in taxa_2 implies
                // taxa_1 is before taxa_2.


                for (int nc = 0; nc < constrainedTaxa.getChildCount(); ++nc) {

                    final Object object = constrainedTaxa.getChild(nc);
                    if (object instanceof XMLObject) {
                        final XMLObject constraint = (XMLObject) object;

                        if (constraint.getName().equals(TMRCA_CONSTRAINT)) {
                            TaxonList taxaSubSet = (TaxonList) constraint.getChild(TaxonList.class);
                            ParametricDistributionModel dist =
                                    (ParametricDistributionModel) constraint.getChild(ParametricDistributionModel.class);
                            boolean isMono = constraint.getAttribute(IS_MONOPHYLETIC, true);

                            final CoalescentSimulator.TaxaConstraint taxaConstraint = new CoalescentSimulator.TaxaConstraint(taxaSubSet, dist, isMono);
                            int insertPoint;
                            for (insertPoint = 0; insertPoint < constraints.size(); ++insertPoint) {
                                // if new <= constraints[insertPoint] insert before insertPoint

                                final CoalescentSimulator.TaxaConstraint iConstraint = constraints.get(insertPoint);
                                if (iConstraint.isMonophyletic) {
                                    if (!taxaConstraint.isMonophyletic) {
                                        continue;
                                    }

                                    final TaxonList taxonsip = iConstraint.taxons;
                                    final int nIn = simulator.sizeOfIntersection(taxonsip, taxaSubSet);
                                    if (nIn == taxaSubSet.getTaxonCount()) {
                                        break;
                                    }
                                    if (nIn > 0 && nIn != taxonsip.getTaxonCount()) {
                                        throw new XMLParseException(setsNotCompatibleMessage);
                                    }
                                } else {
                                    // reached non mono area
                                    if (!taxaConstraint.isMonophyletic) {
                                        if (iConstraint.upper >= taxaConstraint.upper) {
                                            break;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }
                            constraints.add(insertPoint, taxaConstraint);
                        }
                    }
                }
                final int nConstraints = constraints.size();

                if (nConstraints == 0) {
                    if (taxa != null) {
                        taxonLists.add(taxa);
                    }
                } else {
                    for (int nc = 0; nc < nConstraints; ++nc) {
                        CoalescentSimulator.TaxaConstraint cnc = constraints.get(nc);
                        if (!cnc.isMonophyletic) {
                            for (int nc1 = nc - 1; nc1 >= 0; --nc1) {
                                CoalescentSimulator.TaxaConstraint cnc1 = constraints.get(nc1);
                                int x = simulator.sizeOfIntersection(cnc.taxons, cnc1.taxons);
                                if (x > 0) {
                                    Taxa combinedTaxa = new Taxa(cnc.taxons);
                                    combinedTaxa.addTaxa(cnc1.taxons);
                                    cnc = new CoalescentSimulator.TaxaConstraint(combinedTaxa, cnc.lower, cnc.upper, cnc.isMonophyletic);
                                    constraints.set(nc, cnc);
                                }
                            }
                        }
                    }
                    // determine upper bound for each set.
                    double[] upper = new double[nConstraints];
                    for (int nc = nConstraints - 1; nc >= 0; --nc) {
                        final CoalescentSimulator.TaxaConstraint cnc = constraints.get(nc);
                        if (cnc.realLimits()) {
                            upper[nc] = cnc.upper;
                        } else {
                            upper[nc] = Double.POSITIVE_INFINITY;
                        }
                    }

                    for (int nc = nConstraints - 1; nc >= 0; --nc) {
                        final CoalescentSimulator.TaxaConstraint cnc = constraints.get(nc);
                        if (upper[nc] < Double.POSITIVE_INFINITY) {
                            for (int nc1 = nc - 1; nc1 >= 0; --nc1) {
                                final CoalescentSimulator.TaxaConstraint cnc1 = constraints.get(nc1);
                                if (simulator.contained(cnc1.taxons, cnc.taxons)) {
                                    upper[nc1] = Math.min(upper[nc1], upper[nc]);
                                    if (cnc1.realLimits() && cnc1.lower > upper[nc1]) {
                                        throw new XMLParseException(setsNotCompatibleMessage);
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    // collect subtrees here
                    List<Tree> st = new ArrayList<Tree>();
                    for (int nc = 0; nc < constraints.size(); ++nc) {
                        final CoalescentSimulator.TaxaConstraint nxt = constraints.get(nc);
                        // collect all previously built subtrees which are a subset of taxa set to be added
                        List<Tree> subs = new ArrayList<Tree>();
                        Taxa newTaxons = new Taxa(nxt.taxons);
                        for (int k = 0; k < st.size(); ++k) {
                            final Tree stk = st.get(k);
                            int x = simulator.sizeOfIntersection(stk, nxt.taxons);
                            if (x == st.get(k).getTaxonCount()) {
                                final Tree tree = st.remove(k);
                                --k;
                                subs.add(tree);
                                newTaxons.removeTaxa(tree);
                            }
                        }

                        SimpleTree tree = simulator.simulateTree(newTaxons, demoModel);
                        final double lower = nxt.realLimits() ? nxt.lower : 0;
                        if (upper[nc] < Double.MAX_VALUE) {
                            simulator.attemptToScaleTree(tree, (lower + upper[nc]) / 2);
                        }
                        if (subs.size() > 0) {
                            if (tree.getTaxonCount() > 0) subs.add(tree);
                            double h = -1;
                            if (upper[nc] < Double.MAX_VALUE) {
                                for (Tree t : subs) {
                                    // protect against 1 taxa tree with height 0
                                    final double rh = t.getNodeHeight(t.getRoot());
                                    h = Math.max(h, rh > 0 ? rh : (lower + upper[nc]) / 2);
                                }
                                h = (h + upper[nc]) / 2;
                            }
                            tree = simulator.simulateTree(subs.toArray(new Tree[subs.size()]), demoModel, h, true);
                        }
                        st.add(tree);

                    }

                    // add a taxon list for remaining taxa
                    if (taxa != null) {
                        final Taxa list = new Taxa();
                        for (int j = 0; j < taxa.getTaxonCount(); ++j) {
                            Taxon taxonj = taxa.getTaxon(j);
                            for (Tree aSt : st) {
                                if (aSt.getTaxonIndex(taxonj) >= 0) {
                                    taxonj = null;
                                    break;
                                }
                            }
                            if (taxonj != null) {
                                list.addTaxon(taxonj);
                            }
                        }
                        if (list.getTaxonCount() > 0) {
                            taxonLists.add(list);
                        }
                    }
                    if (st.size() > 1) {
                        final Tree t = simulator.simulateTree(st.toArray(new Tree[st.size()]), demoModel, -1, false);
                        subtrees.add(t);
                    } else {
                        subtrees.add(st.get(0));
                    }
                }
            }
        }

        if (taxonLists.size() == 0) {
            if (subtrees.size() == 1) {
                return subtrees.get(0);
            }
            throw new XMLParseException("Expected at least one taxonList or two subtrees in "
                    + getParserName() + " element.");
        }

        try {
            Tree[] trees = new Tree[taxonLists.size() + subtrees.size()];
            // simulate each taxonList separately
            for (int i = 0; i < taxonLists.size(); i++) {
                trees[i] = simulator.simulateTree(taxonLists.get(i), demoModel);
            }
            // add the preset trees
            for (int i = 0; i < subtrees.size(); i++) {
                trees[i + taxonLists.size()] = subtrees.get(i);
            }

            return simulator.simulateTree(trees, demoModel, rootHeight, trees.length != 1);

        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a simulated tree under the given demographic model.";
    }

    public Class getReturnType() {
        return CoalescentSimulator.class; //Object.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(RESCALE_HEIGHT, true, "Attempt to rescale the tree to the given root height"),
            AttributeRule.newDoubleRule(ROOT_HEIGHT, true, ""),
            new ElementRule(Tree.class, 0, Integer.MAX_VALUE),
            new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
            new ElementRule(CONSTRAINED_TAXA, new XMLSyntaxRule[]{
                    new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
                    new ElementRule(TMRCA_CONSTRAINT, new XMLSyntaxRule[]{
                          new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
                    }),
            }, true),
            new ElementRule(DemographicModel.class),
    };
}
