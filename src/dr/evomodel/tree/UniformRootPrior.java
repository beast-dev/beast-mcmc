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

    private static final double tolerance = 1E-6;

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

            // Prune out intervals smaller in length than tolerance
            double intervalStart = tree.getNodeHeight(tree.getRoot());
            List<Double> pruneDates = new ArrayList<Double>();

            for (Double intervalEnd : reversedTipDateList) {
                if (intervalStart - intervalEnd < tolerance) {
                    pruneDates.add(intervalStart);
                }
                intervalStart = intervalEnd;
            }

            for (Double date : pruneDates)
                reversedTipDateList.remove(date);         
        }

        intervalCountsKnown = false;

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


                if (!treePolynomialKnown) {
                    treePolynomial = recursivelyComputePolynomial(tree,tree.getRoot());
                    System.err.println("Final polynomial: "+treePolynomial);
                    treePolynomialKnown = true; // TODO Hit flag when tree updated.
                }

//                System.err.println("treePolynomial = "+treePolynomial);
//                System.exit(-1);
                logLike = -Math.log(treePolynomial.evaluate(rootHeight));

//                // TODO Recalculate only when a tip sampling time changes
//                if (!intervalCountsKnown) {
//
//                    intervals.clear();
//                    for (Double date : reversedTipDateList) {
//                        intervals.put(date, 0);
//                    }
////
////                    traverse(tree, tree.getRoot());
//
//                    for(int i=0; i<tree.getExternalNodeCount(); i++) {
//                        double tipDate = tree.getNodeHeight(tree.getExternalNode(i));
//                        if (!intervals.keySet().contains(tipDate)) {
//                            // Find closest within tolerance
//                            boolean found = false;
//                            for(Double prunedDate : reversedTipDateList) {
//                                if (Math.abs(tipDate - prunedDate) < tolerance) {
//                                    tipDate = prunedDate;
//                                    found = true;
//                                    break;
//                                }
//                            }
//                            if (!found)
//                                throw new RuntimeException("Tip date not found");
//                        }
//                        intervals.put(tipDate, intervals.get(tipDate)+1);
//                    }
//                    intervalCountsKnown = true;
//                }
//
//                logLike = 0.0;
//                for (Double date : reversedTipDateList) {
//                    double s = rootHeight - date;
//                    int k = intervals.get(date);  // There was a bug here, was +=
//
////                    System.err.println("date = "+date+" ("+k+")");
//
//                    k--; // ONLY FOR ABOVE INTERVALS, NOT TRAVERSE
//                    if (k > 0)
//                        logLike += logFactorial(k) - (double) k * Math.log(s);
//
//                }



//                System.exit(-1);

//                Polynomial a = new LabeledPolynomial(new double[] {1,2,3});
//                Polynomial b = new LabeledPolynomial(new double[] {1,0,0,1});
//
//                System.err.println("a = "+a);
//                System.err.println("b = "+b);
//                System.err.println("a * b = "+a.multiply(b));
//                System.err.println("integral a = "+a.integrateWRTX(9));
//                System.exit(-1);


                // Does not work
//                tmpLogLikelihood = 0;
//                recursivelyComputeConditionalUniform(tree, tree.getRoot(), 0.0);
//                logLike = tmpLogLikelihood;

                // New Idea

//                intervals.clear();
//                for (Double date : reversedTipDateList) {
//                    intervals.put(date, 0);
//                }
//                for(int i=0; i<tree.getInternalNodeCount(); i++) {
//                    NodeRef node = tree.getInternalNode(i);
//                    if( !tree.isRoot(node)) {
//                        Double date = inInterval(tree.getNodeHeight(node));
//                        intervals.put(date, intervals.get(date) + 1);
//                    }
//                }
//                logLike = 0.0;
//                double startInterval = rootHeight;
//                for(Double endInterval : reversedTipDateList) {
//                    double length = startInterval - endInterval;
//                    int k = intervals.get(endInterval);
//                    if (k > 0)
//                        logLike += logFactorial(k) - k*Math.log(length) +
//                                   k*(Math.log(length) - Math.log(rootHeight)) -
//                                   logFactorial(k);
//                    startInterval = endInterval;
//                }
            }

            assert !Double.isInfinite(logLike) && !Double.isNaN(logLike);
            return logLike;
        }
    }


    private LabeledPolynomial recursivelyComputePolynomial(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            return new LabeledPolynomial(new double[] {1.0}, tree.getNodeHeight(node));
        }
        LabeledPolynomial childPolynomial1 = recursivelyComputePolynomial(tree, tree.getChild(node,0));
        LabeledPolynomial childPolynomial2 = recursivelyComputePolynomial(tree, tree.getChild(node,1));
        LabeledPolynomial polynomial = childPolynomial1.multiply(childPolynomial2);
        if (!tree.isRoot(node)) {
            polynomial = polynomial.integrateWithLowerBound(polynomial.label);
        }
        return polynomial;
    }

    class Polynomial {

        private static final String X = " x^";

        Polynomial(double[] coefficient) {
            this.coefficient = coefficient;
        }

        Polynomial(Polynomial polynomial) {
            this.coefficient = polynomial.coefficient;
        }

        public int getDegree() {
            return coefficient.length - 1;
        }

        public Polynomial multiply(Polynomial b) {
            double[] newCoefficient = new double[getDegree() + b.getDegree()+1];
            for(int n=0; n<=getDegree(); n++) {
                for(int m=0; m<=b.getDegree(); m++) {
                    newCoefficient[n+m] += coefficient[n] * b.coefficient[m];
                }
            }
            return new Polynomial(newCoefficient);
        }

        public Polynomial integrate() {
            double[] newCoefficient = new double[getDegree()+2];
            for(int n=0; n<=getDegree(); n++) {
                newCoefficient[n+1] = coefficient[n] / (n+1);
            }
            return new Polynomial(newCoefficient);
        }

        public double evaluate(double x) {
            double result = 0;
            double xn = 1;
            for(int n=0; n<=getDegree(); n++) {
                result += xn * coefficient[n];
                xn *= x;
            }
            return result;
        }

        public Polynomial integrateWithLowerBound(double bound) {
            Polynomial integrand = integrate();
            integrand.coefficient[0] = -integrand.evaluate(bound);
            return integrand;
        }

        public String toString() {
            StringBuffer bf = new StringBuffer();
            for(int n=getDegree(); n>=0; n--) {
                bf.append(coefficient[n]);
                bf.append(X);
                bf.append(n);
                if (n > 0)
                    bf.append(" + ");
            }
            return bf.toString();
        }

        double[] coefficient;
    }

    class LabeledPolynomial extends Polynomial {

        LabeledPolynomial(double[] coefficients, double label) {
            super(coefficients);
            this.label = label;
        }

        LabeledPolynomial(double[] coefficients) {
            this(coefficients,0.0);
        }

        LabeledPolynomial(Polynomial polynomial, double label) {
            super(polynomial);
            this.label = label;
        }

        public LabeledPolynomial multiply(LabeledPolynomial b) {
            double maxLabel = Math.max(label,b.label);
            return new LabeledPolynomial(super.multiply(b),maxLabel);
        }

        public LabeledPolynomial integrateWithLowerBound(double bound) {
            return new LabeledPolynomial(super.integrateWithLowerBound(bound),label);
        }

        public String toString() {
            return super.toString() + " {"+label+"}";
        }

        double label;
    }

    private Double inInterval(double date) {
        double startInterval = tree.getNodeHeight(tree.getRoot());
        for(Double endInterval : reversedTipDateList) {
            if (date > endInterval)
                return endInterval;
        }
        throw new RuntimeException("Date in no interval???");
    }

    private double recursivelyComputeConditionalUniform(Tree tree, NodeRef node, double parentHeight) {

        if(tree.isExternal(node))
            return tree.getNodeHeight(node);

        double thisNodeHeight = tree.getNodeHeight(node);
        double childHeight1 = recursivelyComputeConditionalUniform(tree, tree.getChild(node,0), thisNodeHeight);
        double childHeight2 = recursivelyComputeConditionalUniform(tree, tree.getChild(node,1), thisNodeHeight);


        if (!tree.isRoot(node)) {
            double maxChildHeight = (childHeight1 > childHeight2 ? childHeight1 : childHeight2);
            tmpLogLikelihood -= Math.log(parentHeight - maxChildHeight);
        }

        return thisNodeHeight;
    }


    private double logFactorial(int n) {
        if (n == 0 || n == 1) {
            return 0;
        }

        double rValue = 0;

        for (int i = n; i > 1; i--) {
            rValue += Math.log(i);
        }
        return rValue;
    }

//    private void addCountToInterval(double date) {
//        for(Double prunedDate : reversedTipDateList) {
//            if( date > prunedDate) {
//                intervals.put(prunedDate, intervals.get(prunedDate) + 1);
//                break;
//            }
//        }
//    }

    private Double traverse(Tree tree, NodeRef node) {
        Double date;
        if (tree.isExternal(node)) {
            date = tree.getNodeHeight(node);
            if (!intervals.keySet().contains(date)) {
                // Find closest within tolerance
                boolean found = false;
                for(Double prunedDate : reversedTipDateList) {
                    if (Math.abs(date - prunedDate) < tolerance) {
                        date = prunedDate;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    throw new RuntimeException("Tip date not found");
            }

        } else {
            Double date1 = traverse(tree, tree.getChild(node, 0));
            Double date2 = traverse(tree, tree.getChild(node, 1));
            date = (date1 > date2 ? date1 : date2);
            if (!tree.isRoot(node)) {
                // Only increases counts for internal nodes
                intervals.put(date, intervals.get(date) + 1);
            }
        }

        return date;
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
    private boolean intervalCountsKnown = false;
    private boolean treePolynomialKnown = false;
    private double tmpLogLikelihood;
    private LabeledPolynomial treePolynomial;
}