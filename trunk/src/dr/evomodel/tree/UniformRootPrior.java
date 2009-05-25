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

        // Leading coefficient on tree polynomial is X = (# internal nodes)!
        // To keep X > 10E-40, should use log-space polynomials for more than ~30 tips
        if (tree.getExternalNodeCount() < 30) {
            polynomialType = Polynomial.Type.DOUBLE; // Much faster
        } else {
            polynomialType = Polynomial.Type.LOG_DOUBLE;
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

        return calculateLogLikelihood();

//        if (!likelihoodKnown) {
//        	logLikelihood = calculateLogLikelihood();
//        	likelihoodKnown = true;
//        }
//        return logLikelihood;
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

            if (k > 0) {
                // the tips are contemporaneous
//                logLike = logFactorialK - (double) k * Math.log(rootHeight);
                // Try new prior.... should behave the same.
                tmpLogLikelihood = 0;
                recursivelyComputeDensity(tree,tree.getRoot(),0);
                logLike = tmpLogLikelihood;

            } else {
                // TODO Rewrite description above to discuss this new prior
//                if (!treePolynomialKnown) {
//                    polynomialType = Polynomial.Type.MARCRATIONAL;
//                    treePolynomial  = recursivelyComputePolynomial(tree,tree.getRoot(),polynomialType).getPolynomial();
//                    System.err.println("poly1: "+treePolynomial);
//                    System.err.println("eval = "+treePolynomial.logEvaluate(rootHeight));
//                    System.err.println("last: "+treePolynomial.getCoefficientString(0));
//                   Polynomial test = recursivelyComputePolynomial(tree,tree.getRoot(),Polynomial.Type.RATIONAL).getPolynomial();
//                    System.err.println("poly2: "+test);
//                    System.err.println("eval = "+test.logEvaluate(rootHeight));
//                    System.err.println("last: "+test.getCoefficientString(0));
//
//                    System.exit(-1);
//                    treePolynomialKnown = true;
//                }
//
//                logLike  = -treePolynomial .logEvaluate(rootHeight);
//
//                if (Double.isNaN(logLike)) {
//                    // Try using Horner's method
//                    logLike = -treePolynomial.logEvaluateHorner(rootHeight);
//                    if (Double.isNaN(logLike)) {
//                        logLike = Double.NEGATIVE_INFINITY;
//                    }
//                }

                // Try new prior!
                tmpLogLikelihood = 0;
                recursivelyComputeDensity(tree,tree.getRoot(),0);
                logLike = tmpLogLikelihood;
            }

            assert !Double.isInfinite(logLike) && !Double.isNaN(logLike);
            return logLike;
        }
    }


    private double recursivelyComputeDensity(Tree tree, NodeRef node, double parentHeight) {
        if (tree.isExternal(node))
            return tree.getNodeHeight(node);

        double thisHeight = tree.getNodeHeight(node);
        double heightChild1 = recursivelyComputeDensity(tree,tree.getChild(node,0),thisHeight);
        double heightChild2 = recursivelyComputeDensity(tree,tree.getChild(node,1),thisHeight);
        double minHeight = (heightChild1 > heightChild2) ? heightChild1 : heightChild2;

        if (!tree.isRoot(node)) {
            double diff = parentHeight - minHeight;
            if (diff <= 0)
            	tmpLogLikelihood = Double.NEGATIVE_INFINITY;
            else
	            tmpLogLikelihood -= Math.log(diff);
            //tmpLogLikelihood -= Math.log(parentHeight-minHeight);
        } else {
            // Do nothing
        }
        return minHeight;
    }

//    private TipLabeledPolynomial recursivelyComputePolynomial(Tree tree, NodeRef node, Polynomial.Type type) {
//
//        if (tree.isExternal(node)) {
//            double[] value = new double[] {1.0};
//            return new TipLabeledPolynomial(value, tree.getNodeHeight(node), type, true);
//        }
//
//        TipLabeledPolynomial childPolynomial1 = recursivelyComputePolynomial(tree, tree.getChild(node,0),type);
//        TipLabeledPolynomial childPolynomial2 = recursivelyComputePolynomial(tree, tree.getChild(node,1),type);
//        // TODO The partialPolynomial below *should* be cached in an efficient reuse scheme (at least for arbitrary precision)
//        TipLabeledPolynomial polynomial = childPolynomial1.multiply(childPolynomial2);
//        // See AbstractTreeLikelihood for an example of how to flag cached polynomials for re-evaluation
//        System.err.println("B> "+polynomial);
//        if (!tree.isRoot(node)) {
//            polynomial = polynomial.integrateWithLowerBound(polynomial.label);
//            System.err.println("<A "+polynomial);
//        } else {
//            System.err.println("<= ROOT");
//        }
//
//        return polynomial;
//    }

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

//    class TipLabeledPolynomial extends Polynomial.Abstract {
//
//        TipLabeledPolynomial(double[] coefficients, double label, Polynomial.Type type, boolean isTip) {
//            switch (type) {
//                case DOUBLE:        polynomial = new Polynomial.Double(coefficients);
//                                    break;
//                case LOG_DOUBLE:    polynomial = new Polynomial.LogDouble(coefficients);
//                                    break;
//                case BIG_DOUBLE:    polynomial = new Polynomial.BigDouble(coefficients);
//                                    break;
//                case APDOUBLE:      polynomial = new Polynomial.APDouble(coefficients);
//                                    break;
//                case RATIONAL:      polynomial = new Polynomial.RationalDouble(coefficients);
//                                    break;
//                case MARCRATIONAL:  polynomial = new Polynomial.MarcRational(coefficients);
//                                    break;
//
//                default: throw new RuntimeException("Unknown polynomial type");
//            }
//            this.label = label;
//            this.isTip = isTip;
//        }
//
//        TipLabeledPolynomial(Polynomial polynomial, double label, boolean isTip) {
//            this.polynomial = polynomial;
//            this.label = label;
//            this.isTip = isTip;
//        }
//
//        public TipLabeledPolynomial copy() {
//            Polynomial copyPolynomial = polynomial.copy();
//            return new TipLabeledPolynomial(copyPolynomial,this.label,this.isTip);
//        }
//
//        public Polynomial getPolynomial() { return polynomial; }
//
//        public TipLabeledPolynomial multiply(TipLabeledPolynomial b) {
//            double maxLabel = Math.max(label,b.label);
//            return new TipLabeledPolynomial(polynomial.multiply(b),maxLabel,false);
//        }
//
//        public int getDegree() {
//            return polynomial.getDegree();
//        }
//
//        public Polynomial multiply(Polynomial b) {
//            return polynomial.multiply(b);
//        }
//
//        public Polynomial integrate() {
//            return polynomial.integrate();
//        }
//
//        public void expand(double x) {
//            polynomial.expand(x);
//        }
//
//        public double evaluate(double x) {
//            return polynomial.evaluate(x);
//        }
//
//        public double logEvaluate(double x) {
//            return polynomial.logEvaluate(x);
//        }
//
//        public double logEvaluateHorner(double x) {
//            return polynomial.logEvaluateHorner(x);
//        }
//
//        public void setCoefficient(int n, double x) {
//            polynomial.setCoefficient(n,x);
//        }
//
//        public TipLabeledPolynomial integrateWithLowerBound(double bound) {
//            return new TipLabeledPolynomial(polynomial.integrateWithLowerBound(bound),label,isTip);
//        }
//
//        public double getCoefficient(int n) {
//            return polynomial.getCoefficient(n);
//        }
//
//        public String toString() {
//            return polynomial.toString() + " {"+label+"}";
//        }
//
//        public String getCoefficientString(int n) {
//            return polynomial.getCoefficientString(n);
//        }
//
//        private double label;
//        private Polynomial polynomial;
//        private boolean isTip;
//    }

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
    private boolean treePolynomialKnown = false;
    private boolean storedTreePolynomialKnown = false;
    private Polynomial treePolynomial;
    private Polynomial storedTreePolynomial;

    private double tmpLogLikelihood;
//    private Iterator<Polynomial.Type> typeIterator = EnumSet.allOf(Polynomial.Type.class).iterator();
//    private Polynomial.Type polynomialType = typeIterator.next();
    private Polynomial.Type polynomialType;

}