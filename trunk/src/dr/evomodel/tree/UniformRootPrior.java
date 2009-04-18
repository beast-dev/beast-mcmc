/*
 * UniformRootPrior.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.inference.model.*;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;

/**
 * Two priors for the tree that are relatively non-informative on the internal node heights given the root height.
 * The first further assumes that the root height is truncated uniform, see Nicholls, G. & R.D. Gray (2004) for details.
 * The second allows any marginal specification over the root height given that it is larger than the oldest
 * sampling time (Bloomquist and Suchard, unpublished).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Erik Bloomquist
 * @author Marc Suchard
 *
 * @version $Id: UniformRootPrior.java,v 1.10 2005/05/24 20:25:58 rambaut Exp $
 */
public class UniformRootPrior extends AbstractModelLikelihood {

    // PUBLIC STUFF

    public static final String UNIFORM_ROOT_PRIOR = "uniformRootPrior";
    public static final String MAX_ROOT_HEIGHT = "maxRootHeight";

    private int k = 0;
    private double logFactorialK;

    private double maxRootHeight;

    private boolean isNicholls;

    Set<Double> tipDates = new TreeSet<Double>();
    List<Double> reversedTipDateList = new ArrayList<Double>();
    Map<Double, Integer> intervals = new TreeMap<Double, Integer>();

    public UniformRootPrior(Tree tree) {
        this(UNIFORM_ROOT_PRIOR, tree);
    }

    private UniformRootPrior(String name, Tree tree) {

        super(name);

        this.tree = tree;
        this.isNicholls = false;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double h = tree.getNodeHeight(tree.getExternalNode(i));
            tipDates.add(h);
        }

        if (tipDates.size() == 1) {
            // the tips are contemporaneous so these are constant...
            k = tree.getInternalNodeCount() - 1;
            Logger.getLogger("dr.evomodel").info("Uniform Root Prior, Intervals = " + (k + 1));

            logFactorialK = logFactorial(k);
        } else {
            reversedTipDateList.addAll(tipDates);
            Collections.reverse(reversedTipDateList);
        }

    }

    public UniformRootPrior(Tree tree, double maxRootHeight) {
        this(UNIFORM_ROOT_PRIOR, tree, maxRootHeight);
    }

    private UniformRootPrior(String name, Tree tree, double maxRootHeight) {

        super(name);

        this.tree = tree;
        this.maxRootHeight = maxRootHeight;
        isNicholls = true;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }

    UniformRootPrior(String name) {
        super(name);
    }

    // **************************************************************
    // Extendable methods
    // **************************************************************

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected final void handleModelChangedEvent(Model model, Object object, int index) {

        likelihoodKnown = false;
    }

    // **************************************************************
    // ParameterListener IMPLEMENTATION
    // **************************************************************

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected final void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected final void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public double getLogLikelihood() {

        return calculateLogLikelihood();

        //if (!likelihoodKnown) {
        //	logLikelihood = calculateLogLikelihood();
        //	likelihoodKnown = true;
        //}
        //return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
    }

    public double calculateLogLikelihood() {

        double rootHeight = tree.getNodeHeight(tree.getRoot());

        if (isNicholls) {
            int nodeCount = tree.getExternalNodeCount();

            if (rootHeight < 0 || rootHeight > (0.999 * maxRootHeight)) return Double.NEGATIVE_INFINITY;

            // from Nicholls, G. & R.D. Gray (2004)
            return rootHeight * (2 - nodeCount) - Math.log(maxRootHeight - rootHeight);

        } else {
            // the Bloomquist & Suchard variant
            // Let the sampling times and rootHeight specify the boundaries between a fixed number of intervals.
            // Internal node heights are equally likely to fall in any of these intervals and uniformly distributed
            // in an interval before sorting (i.e. the intercoalescent times in an interval form a scaled Dirchelet(1,1,\ldots,1)
            // This is a conditional density on the rootHeight, so it is possible to specify a marginal distribution
            // on the rootHeight given it is greater than the oldest sampling time.
            double logLike;

            if (k > 0) {
                // the tips are contemporaneous
                logLike = logFactorialK - (double) k * Math.log(rootHeight);

            } else {
//                int k1 = -1; // the root will always be > maxRootHeight
//                for (int i = 0; i < tree.getInternalNodeCount(); i++) {
//                    double h = tree.getNodeHeight(tree.getInternalNode(i));
//                    if (h > maxTipHeight) {
//                        k1++;
//                    }
//                }
//
//                logLike = logFactorial(k1) - (double) k1 * Math.log(rootHeight - maxTipHeight);
                intervals.clear();
                for (Double date : tipDates) {
                    intervals.put(date, 0);
                }

                traverse(tree, tree.getRoot());

                logLike = 0.0;
                int k = 0;
                for (Double date : reversedTipDateList) {
                    double s = rootHeight - date;
                    k += intervals.get(date);

                    logLike += logFactorial(k) - (double) k * Math.log(s);
                }
            }

            assert !Double.isInfinite(logLike) && !Double.isNaN(logLike);
            return logLike;
        }
    }


    private double logFactorial(int n) {
        if (n == 0) {
            return 0;
        }

        double rValue = 0;

        for (int i = n; i > 0; i--) {
            rValue += Math.log(i);
        }
        return rValue;
    }


    private Double traverse(Tree tree, NodeRef node) {
        Double date;
        if (tree.isExternal(node)) {
            date = tree.getNodeHeight(node);
            if (!intervals.keySet().contains(date)) {
                throw new RuntimeException("Tip date not found");
            }

        } else {
            Double date1 = traverse(tree, tree.getChild(node, 0));
            Double date2 = traverse(tree, tree.getChild(node, 1));
            date = (date1 > date2 ? date1 : date2);
            if (!tree.isRoot(node)) {
                intervals.put(date, intervals.get(date) + 1);
            }
        }

        return date;
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public final dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return UNIFORM_ROOT_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            if (xo.hasAttribute(MAX_ROOT_HEIGHT)) {
                // the Nicholls & Gray variant
                double maxRootHeight = xo.getDoubleAttribute(MAX_ROOT_HEIGHT);
                return new UniformRootPrior(treeModel, maxRootHeight);
            } else {
                // the Bloomquist & Suchard variant
                return new UniformRootPrior(treeModel);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the demographic function.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MAX_ROOT_HEIGHT, true),
                new ElementRule(TreeModel.class)
        };
    };

    /**
     * The tree.
     */
    Tree tree = null;

    double logLikelihood;
    private double storedLogLikelihood;
    boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}