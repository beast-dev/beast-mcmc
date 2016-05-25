/*
 * UniformNodeHeightPrior.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodelxml.tree.UniformNodeHeightPriorParser;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.math.Polynomial;

import java.util.*;
import java.util.logging.Logger;

//import org.jscience.mathematics.number.Rational;

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
 * @version $Id: UniformRootPrior.java,v 1.10 2005/05/24 20:25:58 rambaut Exp $
 */
public class UniformNodeHeightPrior extends AbstractModelLikelihood {

    // PUBLIC STUFF
    public static final int MAX_ANALYTIC_TIPS = 60; // TODO Determine this value!
    public static final int DEFAULT_MC_SAMPLE = 100000;

    private static final double tolerance = 1E-6;

    private int k = 0;
    private double logFactorialK;

    private double maxRootHeight;

    private boolean isNicholls;
    private boolean useAnalytic;
    private boolean useMarginal;
    private boolean leadingTerm;
    private int mcSampleSize;

    Set<Double> tipDates = new TreeSet<Double>();
    List<Double> reversedTipDateList = new ArrayList<Double>();
    Map<Double, Integer> intervals = new TreeMap<Double, Integer>();

    public UniformNodeHeightPrior(Tree tree, boolean useAnalytic, boolean marginal, boolean leadingTerm) {
        this(UniformNodeHeightPriorParser.UNIFORM_NODE_HEIGHT_PRIOR, tree, useAnalytic, DEFAULT_MC_SAMPLE, marginal, leadingTerm);
    }

    public UniformNodeHeightPrior(Tree tree, boolean useAnalytic, int mcSampleSize) {
        this(UniformNodeHeightPriorParser.UNIFORM_NODE_HEIGHT_PRIOR,tree,useAnalytic,mcSampleSize, false, false);
    }

    private UniformNodeHeightPrior(String name, Tree tree, boolean useAnalytic, int mcSampleSize,
                                   boolean marginal, boolean leadingTerm) {

        super(name);

        this.tree = tree;
        this.isNicholls = false;
        this.useAnalytic = useAnalytic;
        this.useMarginal = marginal;
        this.mcSampleSize = mcSampleSize;
        this.leadingTerm = leadingTerm;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double h = tree.getNodeHeight(tree.getExternalNode(i));
            tipDates.add(h);
        }

        if (tipDates.size() == 1 || leadingTerm) {
            // the tips are contemporaneous so these are constant...
            k = tree.getInternalNodeCount() - 1;
            Logger.getLogger("dr.evomodel").info("Uniform Node Height Prior, Intervals = " + (k + 1));

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

            if (!useAnalytic) {
                logLikelihoods = new double[mcSampleSize];
                drawNodeHeights = new double[tree.getNodeCount()][mcSampleSize];
                minNodeHeights = new double[tree.getNodeCount()];
            }
        }

        // Leading coefficient on tree polynomial is X = (# internal nodes)!
        // To keep X > 10E-40, should use log-space polynomials for more than ~30 tips
        if (tree.getExternalNodeCount() < 30) {
            polynomialType = Polynomial.Type.DOUBLE; // Much faster
        } else if (tree.getExternalNodeCount() < 45){
            polynomialType = Polynomial.Type.LOG_DOUBLE;
        } else {
//            polynomialType = Polynomial.Type.APDOUBLE;
            polynomialType = Polynomial.Type.LOG_DOUBLE;
        }
        Logger.getLogger("dr.evomodel").info("Using "+polynomialType+" polynomials!");

    }

    public UniformNodeHeightPrior(Tree tree, double maxRootHeight) {
        this(UniformNodeHeightPriorParser.UNIFORM_NODE_HEIGHT_PRIOR, tree, maxRootHeight);
    }

    private UniformNodeHeightPrior(String name, Tree tree, double maxRootHeight) {

        super(name);

        this.tree = tree;
        this.maxRootHeight = maxRootHeight;
        isNicholls = true;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }

    UniformNodeHeightPrior(String name) {
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
        treePolynomialKnown = false;
        return;

        // Only set treePolynomialKnown = false when a topology change occurs
        // Only set likelihoodKnown = false when a topology change occurs or the rootHeight is changed

//        if (model == tree) {
//            if (object instanceof TreeModel.TreeChangedEvent) {
//                TreeModel.TreeChangedEvent event = (TreeModel.TreeChangedEvent) object;
//                if (event.isHeightChanged()) {
//                    if (event.getNode() == tree.getRoot()) {
//                        likelihoodKnown = false;
//                        return;
//                    } // else
//                    return;
//                }
//                if (event.isNodeParameterChanged())
//                    return;
//                // All others are probably tree structure changes
//                likelihoodKnown = false;
//                treePolynomialKnown = false;
//                return;
//            }
//            // TODO Why are not all node height changes invoking TreeChangedEvents?
//            if (object instanceof Parameter.Default) {
//                Parameter parameter = (Parameter) object;
//                if (tree.getNodeHeight(tree.getRoot()) == parameter.getParameterValue(index)) {
//                    likelihoodKnown = false;
//                    treePolynomialKnown = false;
//                    return;
//                }
//                return;
//            }
//        }
//
//        throw new RuntimeException("Unexpected event!");
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
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
//        storedTreePolynomialKnown = treePolynomialKnown;
//        if (treePolynomial != null)
//            storedTreePolynomial = treePolynomial.copy(); // TODO Swap pointers
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected final void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
//        treePolynomialKnown = storedTreePolynomialKnown;
//        treePolynomial = storedTreePolynomial;
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

//        return calculateLogLikelihood();

        if (!likelihoodKnown) {
        	logLikelihood = calculateLogLikelihood();
        	likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
        treePolynomialKnown = false;
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

            if (k > 0) {    // Also valid for leading-term approximation
                // the tips are contemporaneous
                logLike = logFactorialK - (double) k * Math.log(rootHeight);

//                double cutoff = 62;
//                int count = 0;
//                for (int i = 0; i < tree.getNodeCount(); i++) {
//                    if (tree.getNodeHeight(tree.getNode(i)) >  cutoff) {
//                        count++;
////                        if (tree.isExternal(tree.getNode(i))) {
////                            System.err.println("Problem");
////                            System.exit(-1);
////                        }
//                    }
//                }
//                count -= 1; // ignore root
////                System.err.println("c = " + count);
//                logLike = logFactorial(count) - (double) count * Math.log(rootHeight - cutoff);
//                logLike = logFactorial(k - count) - (double) (k - count) * Math.log(cutoff);
                
            } else {
                // TODO Rewrite description above to discuss this new prior
                if (useAnalytic) {

//                    long startTime1 = System.nanoTime();
                    if (useMarginal) {

                        if (!treePolynomialKnown) {
//                            treePolynomial = recursivelyComputePolynomial(tree, tree.getRoot(), polynomialType).getPolynomial();
                            treePolynomials = constructRootPolyonmials(tree,polynomialType); // Each polynomial is of lower degree
                            treePolynomialKnown = true;
                        }

//                        logLike = -treePolynomial.logEvaluate(rootHeight);
                        logLike = -treePolynomials[0].logEvaluate(rootHeight) - treePolynomials[1].logEvaluate(rootHeight);

                        if (Double.isNaN(logLike)) {
                            // Try using Horner's method
//                            logLike = -treePolynomial.logEvaluateHorner(rootHeight); // TODO this could be causing the problem!!
                            logLike = -treePolynomials[0].logEvaluateHorner(rootHeight) - treePolynomials[1].logEvaluateHorner(rootHeight);
                            if (Double.isNaN(logLike)) {
                                logLike = Double.NEGATIVE_INFINITY;
                            }
                        }
                    } else {

                        tmpLogLikelihood = 0;
                        recursivelyComputeDensity(tree, tree.getRoot(), 0);
                        logLike = tmpLogLikelihood;

                    }

//                    long stopTime1 = System.nanoTime();


                } else {

//                    long startTime2 = System.nanoTime();

                   // Copy over current root height
                    final double[] drawRootHeight = drawNodeHeights[tree.getRoot().getNumber()];
                    Arrays.fill(drawRootHeight,rootHeight); // TODO Only update when rootHeight changes

                    // Determine min heights for each node in tree
                    recursivelyFindNodeMinHeights(tree,tree.getRoot()); // TODO Only update when topology changes

                    // Simulate from prior
                    Arrays.fill(logLikelihoods,0.0);
                    recursivelyComputeMCIntegral(tree, tree.getRoot(), tree.getRoot().getNumber()); // TODO Only update when topology or rootHeight changes

                    // Take average
                    logLike = -LogTricks.logSum(logLikelihoods) + Math.log(mcSampleSize);

//                    long stopTime2 = System.nanoTime();

                }
            }

            assert !Double.isInfinite(logLike) && !Double.isNaN(logLike);
            return logLike;
        }
    }

//    Map<Double,Integer> boxCounts;
//
//    private double recursivelyComputeMarcDensity(Tree tree, NodeRef node, double rootHeight) {
//        if (tree.isExternal(node))
//            return tree.getNodeHeight(node);
//
////        double thisHeight = tree.getNodeHeight(node);
////        double thisHeight = rootHeight;
//        double heightChild1 = recursivelyComputeMarcDensity(tree, tree.getChild(node, 0), rootHeight);
//        double heightChild2 = recursivelyComputeMarcDensity(tree, tree.getChild(node, 1), rootHeight);
//        double minHeight = (heightChild1 > heightChild2) ? heightChild1 : heightChild2;
//
//        if (!tree.isRoot(node)) {
//            double diff = rootHeight - minHeight;
//            if (diff <= 0)
//                tmpLogLikelihood = Double.NEGATIVE_INFINITY;
//            else
//                tmpLogLikelihood -= Math.log(diff);
//
//            Integer count = boxCounts.get(minHeight);
//            if (count == null) {
//                boxCounts.put(minHeight,1);
////                System.err.println("new height: "+minHeight);
//            } else {
//                boxCounts.put(minHeight,count+1);
////                System.err.println("old height: "+minHeight);
//            }
//            // TODO Could do the logFactorial right here
//        } else {
//            // Do nothing
//        }
//        return minHeight;
//    }


    private double recursivelyComputeDensity(Tree tree, NodeRef node, double parentHeight) {
        if (tree.isExternal(node))
            return tree.getNodeHeight(node);

        double thisHeight = tree.getNodeHeight(node);
        double heightChild1 = recursivelyComputeDensity(tree, tree.getChild(node, 0), thisHeight);
        double heightChild2 = recursivelyComputeDensity(tree, tree.getChild(node, 1), thisHeight);
        double minHeight = (heightChild1 > heightChild2) ? heightChild1 : heightChild2;

        if (!tree.isRoot(node)) {
            double diff = parentHeight - minHeight;
            if (diff <= 0)
                tmpLogLikelihood = Double.NEGATIVE_INFINITY;
            else
                tmpLogLikelihood -= Math.log(diff);
//                tmpLogLikelihood -= Math.log(parentHeight-minHeight);
        } else {
            // Do nothing
        }
        return minHeight;
    }

    private double recursivelyFindNodeMinHeights(Tree tree, NodeRef node) {

        double minHeight;

        if (tree.isExternal(node))
            minHeight = tree.getNodeHeight(node);
        else {
            double minHeightChild0 = recursivelyFindNodeMinHeights(tree, tree.getChild(node,0));
            double minHeightChild1 = recursivelyFindNodeMinHeights(tree, tree.getChild(node,1));
            minHeight = (minHeightChild0 > minHeightChild1) ? minHeightChild0 : minHeightChild1;
        }

        minNodeHeights[node.getNumber()] = minHeight;
        return minHeight;        
    }

    private void recursivelyComputeMCIntegral(Tree tree, NodeRef node, int parentNodeNumber) {

        if (tree.isExternal(node))
            return;

        final int nodeNumber = node.getNumber();

        if (!tree.isRoot(node)) {

            final double[] drawParentHeight = drawNodeHeights[parentNodeNumber];
            final double[] drawThisNodeHeight = drawNodeHeights[nodeNumber];
            final double minHeight = minNodeHeights[nodeNumber];

            final boolean twoChild = (tree.isExternal(tree.getChild(node,0)) && tree.isExternal(tree.getChild(node,1)));

            for(int i=0; i<mcSampleSize; i++) {

                final double diff = drawParentHeight[i] - minHeight;
                if (diff <= 0) {
                    logLikelihoods[i] = Double.NEGATIVE_INFINITY;
                    break;
                }

                if (!twoChild)
                    drawThisNodeHeight[i] = MathUtils.nextDouble() * diff + minHeight;

                logLikelihoods[i] += Math.log(diff);
            }
        }

        recursivelyComputeMCIntegral(tree, tree.getChild(node,0), nodeNumber);
        recursivelyComputeMCIntegral(tree, tree.getChild(node,1), nodeNumber);

    }

    private static final double INV_PRECISION = 10;

    private static double round(double x) {
        return Math.round(x * INV_PRECISION) / INV_PRECISION;
    }

    private Polynomial[] constructRootPolyonmials(Tree tree, Polynomial.Type type) {
        NodeRef root = tree.getRoot();
        return new Polynomial[] {
                recursivelyComputePolynomial(tree,tree.getChild(root,0),type).getPolynomial(),
                recursivelyComputePolynomial(tree,tree.getChild(root,1),type).getPolynomial()
        };
    }

    private TipLabeledPolynomial recursivelyComputePolynomial(Tree tree, NodeRef node, Polynomial.Type type) {

        if (tree.isExternal(node)) {
            double[] value = new double[]{1.0};
            double height = round(tree.getNodeHeight(node)); // Should help in numerical stability        
            return new TipLabeledPolynomial(value, height, type, true);
        }

        TipLabeledPolynomial childPolynomial1 = recursivelyComputePolynomial(tree, tree.getChild(node, 0), type);
        TipLabeledPolynomial childPolynomial2 = recursivelyComputePolynomial(tree, tree.getChild(node, 1), type);
        // TODO The partialPolynomial below *should* be cached in an efficient reuse scheme (at least for arbitrary precision)
        TipLabeledPolynomial polynomial = childPolynomial1.multiply(childPolynomial2);
        // See AbstractTreeLikelihood for an example of how to flag cached polynomials for re-evaluation
        if (!tree.isRoot(node)) {
            polynomial = polynomial.integrateWithLowerBound(polynomial.label);
        }

        return polynomial;
    }

//    private void test() {
//
//        double[] value = new double[]{2, 0, 2};
//        Polynomial a = new Polynomial.Double(value);
//        Polynomial a2 = a.multiply(a);
//        System.err.println("a  :" + a);
//        System.err.println("a*a: " + a2);
//        System.err.println("eval :" + a2.evaluate(2));
//        Polynomial intA = a.integrate();
//        System.err.println("intA: " + intA);
//        Polynomial intA2 = a.integrateWithLowerBound(2.0);
//        System.err.println("intA2: " + intA2);
//        System.err.println("");
//
//        Polynomial b = new Polynomial.APDouble(value);
//        System.err.println("b : " + b);
//        Polynomial b2 = b.multiply(b);
//        System.err.println("b2 : " + b2);
//        System.err.println("eval : " + b2.evaluate(2));
//        Polynomial intB = b.integrate();
//        System.err.println("intB: " + intB);
//        Polynomial intB2 = b.integrateWithLowerBound(2.0);
//        System.err.println("intB2: " + intB2);
//        System.err.println("");
//
//        Polynomial c = new Polynomial.LogDouble(value);
//        System.err.println("c : " + c);
//        Polynomial c2 = c.multiply(c);
//        System.err.println("c2 : " + c2);
//        System.err.println("eval : " + c2.evaluate(2));
//        Polynomial intC = c.integrate();
//        System.err.println("intC: " + intC);
//        Polynomial intC2 = c.integrateWithLowerBound(2.0);
//        System.err.println("intC2: " + intC2);
//        System.exit(-1);
//    }

    class TipLabeledPolynomial extends Polynomial.Abstract {

        TipLabeledPolynomial(double[] coefficients, double label, Polynomial.Type type, boolean isTip) {
            switch (type) {
                case DOUBLE:
                    polynomial = new Polynomial.Double(coefficients);
                    break;
                case LOG_DOUBLE:
                    polynomial = new Polynomial.LogDouble(coefficients);
                    break;
                case BIG_DOUBLE:
                    polynomial = new Polynomial.BigDouble(coefficients);
                    break;
//                case APDOUBLE:      polynomial = new Polynomial.APDouble(coefficients);
//                    break;
//                case RATIONAL:      polynomial = new Polynomial.RationalDouble(coefficients);
//                    break;
//                case MARCRATIONAL:  polynomial = new Polynomial.MarcRational(coefficients);
//                                    break;

                default:
                    throw new RuntimeException("Unknown polynomial type");
            }
            this.label = label;
            this.isTip = isTip;
        }

        TipLabeledPolynomial(Polynomial polynomial, double label, boolean isTip) {
            this.polynomial = polynomial;
            this.label = label;
            this.isTip = isTip;
        }

        public TipLabeledPolynomial copy() {
            Polynomial copyPolynomial = polynomial.copy();
            return new TipLabeledPolynomial(copyPolynomial, this.label, this.isTip);
        }

        public Polynomial getPolynomial() {
            return polynomial;
        }

        public TipLabeledPolynomial multiply(TipLabeledPolynomial b) {
            double maxLabel = Math.max(label, b.label);
            return new TipLabeledPolynomial(polynomial.multiply(b), maxLabel, false);
        }

        public int getDegree() {
            return polynomial.getDegree();
        }

        public Polynomial multiply(Polynomial b) {
            return polynomial.multiply(b);
        }

        public Polynomial integrate() {
            return polynomial.integrate();
        }

        public void expand(double x) {
            polynomial.expand(x);
        }

        public double evaluate(double x) {
            return polynomial.evaluate(x);
        }

        public double logEvaluate(double x) {
            return polynomial.logEvaluate(x);
        }

        public double logEvaluateHorner(double x) {
            return polynomial.logEvaluateHorner(x);
        }

        public void setCoefficient(int n, double x) {
            polynomial.setCoefficient(n, x);
        }

        public TipLabeledPolynomial integrateWithLowerBound(double bound) {
            return new TipLabeledPolynomial(polynomial.integrateWithLowerBound(bound), label, isTip);
        }

        public double getCoefficient(int n) {
            return polynomial.getCoefficient(n);
        }

        public String toString() {
            return polynomial.toString() + " {" + label + "}";
        }

        public String getCoefficientString(int n) {
            return polynomial.getCoefficientString(n);
        }

        private double label;
        private Polynomial polynomial;
        private boolean isTip;
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

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /**
     * The tree.
     */
    Tree tree = null;

    double logLikelihood;
    private double storedLogLikelihood;
    boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private boolean treePolynomialKnown = false;
    private boolean storedTreePolynomialKnown = false;
    private Polynomial treePolynomial;
    private Polynomial[] treePolynomials;
    private Polynomial storedTreePolynomial;

    private double tmpLogLikelihood;
    //    private Iterator<Polynomial.Type> typeIterator = EnumSet.allOf(Polynomial.Type.class).iterator();
    //    private Polynomial.Type polynomialType = typeIterator.next();
    private Polynomial.Type polynomialType;

    private double[] logLikelihoods;
    private double[][] drawNodeHeights;
    private double[] minNodeHeights;

}