package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class IntegratedMultivariateTraitLikelihood extends AbstractMultivariateTraitLikelihood {

    public static final String TRAIT_LIKELIHOOD = "integratedMultivariateTraitLikelihood";

    public static final double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 TreeModel treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel, Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 MultivariateNormalDistribution rootPrior) {
        super(traitName, treeModel, diffusionModel, traitParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate);

        dimTrait = diffusionModel.getPrecisionmatrix().length;
        dim = treeModel.getMultivariateNodeTrait(treeModel.getExternalNode(0), traitName).length;


        if (dim % dimTrait != 0)
            throw new RuntimeException("dim is not divisible by dimTrait");

        dimData = dim / dimTrait;

        meanCache = new double[dim * treeModel.getNodeCount()];
        upperPrecisionCache = new double[dim * treeModel.getNodeCount()];
        lowerPrecisionCache = new double[dim * treeModel.getNodeCount()];
        logRemainderDensityCache = new double[dim * treeModel.getNodeCount()];

        // Set tip data values
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            int index = node.getNumber();
            double[] traitValue = treeModel.getMultivariateNodeTrait(node, traitName);
            System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        }

        setRootPrior(rootPrior);

        StringBuffer sb = new StringBuffer();
        sb.append("\tDiffusion dimension: ").append(dimTrait).append("\n");
        sb.append("\tNumber of observations: ").append(dimData).append("\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());


//        System.err.println("starting mean cache: "+new Vector(meanCache));
//        System.err.println("hello?");
//        System.exit(-1);
    }

    private void setRootPrior(MultivariateNormalDistribution rootPrior) {
        rootPriorMean = rootPrior.getMean();
        rootPriorPrecision = rootPrior.getScaleMatrix();

//        System.err.println("root prior mean: " + new Vector(rootPriorMean));
//        System.err.println("root prior prec: \n" + new Matrix(rootPriorPrecision));

        try {
            rootPriorPrecisionDeterminant = new Matrix(rootPriorPrecision).determinant();
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
    }

    protected String extraInfo() {
        return "\tSample internal node traits: false\n";
    }

    public double getLogDataLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public double calculateLogLikelihood() {

        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double logDetTreePrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());
        postOrderTraverse(treeModel, treeModel.getRoot(), treePrecision, logDetTreePrecision);

        // Integrate root trait out against rootPrior
        final int rootIndex = treeModel.getRoot().getNumber();

        System.err.println("mean: " + new Vector(meanCache));
        System.err.println("upre: " + new Vector(upperPrecisionCache));
        System.err.println("lpre: " + new Vector(lowerPrecisionCache));

        double logLikelihood = 0;

        for (int datum = 0; datum < dimData; datum++) {

            double[] rootMean = new double[dimTrait]; // TODO Allocate once
            System.arraycopy(meanCache, rootIndex + datum*dimTrait, rootMean, 0, dimTrait);
            double rootPrecision = lowerPrecisionCache[rootIndex];

            System.err.println("root mean: " + new Vector(rootMean));
            System.err.println("root prec: " + rootPrecision);

            double[][] diffusionPrecision = diffusionModel.getPrecisionmatrix();
            System.err.println("diffusion prec: " + new Matrix(diffusionPrecision));

            // B = root precision
            // z = root mean
            // A = likelihood precision
            // y = likelihood mean

            double[] Bz = new double[dimTrait]; // TODO Allocate once
            double[] Ay = new double[dimTrait]; // TODO Allocate once
            double[][] AplusB = new double[dimTrait][dimTrait]; // TODO Allocate once

            // z'Bz -- TODO only need to calculate once
            double zBz = 0;
            for (int i = 0; i < dimTrait; i++) {
                Bz[i] = 0;
                for (int j = 0; j < dimTrait; j++)
                    Bz[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
                zBz += rootPriorMean[i] * Bz[i];
            }

            // y'Ay
            double yAy = 0;
            for (int i = 0; i < dimTrait; i++) {
                Ay[i] = 0;
                for (int j = 0; j < dimTrait; j++)
                    Ay[i] += diffusionPrecision[i][j] * rootMean[j] * rootPrecision;
                yAy += rootMean[i] * Ay[i];
            }

            // (Ay + Bz)' (A+B)^{-1} (Ay + Bz)
            for (int i = 0; i < dimTrait; i++) {
                Ay[i] += Bz[i];   // Ay is filled with sum
                for (int j = 0; j < dimTrait; j++) {
                    AplusB[i][j] = diffusionPrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
                }
            }

            Matrix mat = new Matrix(AplusB);
            double detAplusB = 0;
            double detDiffusion = diffusionModel.getDeterminantPrecisionMatrix();

            try {
                detAplusB = mat.determinant();

            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

            double square = 0;
            double[][] invAplusB = mat.inverse().toComponents();
            for (int i = 0; i < dimTrait; i++) {
                for (int j = 0; j < dimTrait; j++)
                    square += Ay[i] * invAplusB[i][j] * Ay[j];
            }

//            System.err.println("zBz = " + zBz);
//            System.err.println("yAy = " + yAy);
//            System.err.println("-(Ay+Bz)(A+B)^{1}(Ay+Bz) = " + square);

            logLikelihood += -LOG_SQRT_2_PI * dimTrait
                    + 0.5 * (Math.log(detDiffusion) + dimTrait * Math.log(rootPrecision) +
                    Math.log(rootPriorPrecisionDeterminant) - Math.log(detAplusB))
                    - 0.5 * (yAy + zBz - square);

        }

        double sumLogRemainders = 0;
        for(int i=0; i<logRemainderDensityCache.length; i++) // TODO Skip remainders for tips and root
            sumLogRemainders += logRemainderDensityCache[i];

        logLikelihood += sumLogRemainders;

        return logLikelihood;
    }

    void postOrderTraverse(TreeModel treeModel, NodeRef node, double[][] precisionMatrix, double logDetPrecisionMatrix) {

        if (treeModel.isExternal(node)) {
            // Fill in precision scalar, traitValues already filled in
            upperPrecisionCache[node.getNumber()] = 1.0 / getRescaledBranchLength(node);
            lowerPrecisionCache[node.getNumber()] = Double.POSITIVE_INFINITY;
            return;
        }

        final NodeRef childNode0 = treeModel.getChild(node, 0);
        final NodeRef childNode1 = treeModel.getChild(node, 1);

        postOrderTraverse(treeModel, childNode0, precisionMatrix, logDetPrecisionMatrix);
        postOrderTraverse(treeModel, childNode1, precisionMatrix, logDetPrecisionMatrix);

        final int childNumber0 = childNode0.getNumber();
        final int childNumber1 = childNode1.getNumber();
        final int thisNumber = node.getNumber();
        final int meanOffset0 = dim * childNumber0;
        final int meanOffset1 = dim * childNumber1;
        final int meanThisOffset = dim * thisNumber;

        final double precision0 = upperPrecisionCache[childNumber0];
        final double precision1 = upperPrecisionCache[childNumber1];
        final double totalPrecision = precision0 + precision1;

        lowerPrecisionCache[thisNumber] = totalPrecision;

        // Multiple child0 and child1 densities
        for (int i = 0; i < dim; i++)
            meanCache[meanThisOffset + i] = (meanCache[meanOffset0 + i] * precision0 +
                    meanCache[meanOffset1 + i] * precision1)
                    / totalPrecision;

        if (!treeModel.isRoot(node)) {

            // Integrate out trait value at this node
            double thisPrecision = 1.0 / getRescaledBranchLength(node);
            upperPrecisionCache[thisNumber] = totalPrecision * thisPrecision / (totalPrecision + thisPrecision);
        }

        // Compute logRemainderDensity

        logRemainderDensityCache[thisNumber] = 0;

        for(int k=0; k<dimData; k++) {

            double childSS0 = 0;
            double childSS1 = 0;
            double crossSS  = 0;

            final double remainderPrecision = precision0 * precision1 / (precision0 + precision1);

            for(int i=0; i<dimTrait; i++) {

                final double wChild0i = meanCache[meanOffset0+i] * precision0;
                final double wChild1i = meanCache[meanOffset1+i] * precision1;

                for(int j=0; j<dimTrait; j++) {

                    final double child0j = meanCache[meanOffset0+j];
                    final double child1j = meanCache[meanOffset1+j];

                    childSS0 += wChild0i * precisionMatrix[i][j] * child0j;
                    childSS1 += wChild1i * precisionMatrix[i][j] * child1j;

                    crossSS += (wChild0i + wChild1i) * precisionMatrix[i][j] * meanCache[meanThisOffset + j];
                }

                logRemainderDensityCache[thisNumber] +=
                        - dimTrait * LOG_SQRT_2_PI
                        + 0.5 * dimTrait * (Math.log(remainderPrecision)+logDetPrecisionMatrix)
                        - 0.5 * (childSS0 + childSS1 - crossSS);
            }

        }

    }

    class MeanPrecision {
        MeanPrecision(double[] mean, double precisionScalar) {
            this.mean = mean;
            this.precisionScalar = precisionScalar;
        }

        double[] mean;
        double precisionScalar;
    }

    protected double[] traitForNode(TreeModel tree, NodeRef node, String traitName) {
        // TODO Must to a pre-order traversal to draw values, see AncestralStateTreeLikelihood
        return new double[0];
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRAIT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            CompoundParameter traitParameter = (CompoundParameter) xo.getElementFirstChild(TRAIT_PARAMETER);
//            CompoundParameter traitParameter = null;


            boolean cacheBranches = xo.getAttribute(CACHE_BRANCHES, false);
//            boolean integrate = xo.getAttribute(INTEGRATE,false);
//            boolean integrate = true;

            BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

            List<Integer> missingIndices = null;
            String traitName = DEFAULT_TRAIT_NAME;

            if (xo.hasAttribute(TRAIT_NAME)) {

                traitName = xo.getStringAttribute(TRAIT_NAME);


                // Fill in attributeValues
                int taxonCount = treeModel.getTaxonCount();
                for (int i = 0; i < taxonCount; i++) {
                    String taxonName = treeModel.getTaxonId(i);
                    String paramName = taxonName + "." + traitName;
                    Parameter traitParam = getTraitParameterByName(traitParameter, paramName);
                    if (traitParam == null)
                        throw new RuntimeException("Missing trait parameters at tree tips");
                    String object = (String) treeModel.getTaxonAttribute(i, traitName);
                    if (object == null)
                        throw new RuntimeException("Trait \"" + traitName + "\" not found for taxa \"" + taxonName + "\"");
                    else {
                        StringTokenizer st = new StringTokenizer(object);
                        int count = st.countTokens();
                        if (count != traitParam.getDimension())
                            throw new RuntimeException("Trait length must match trait parameter dimension");
                        for (int j = 0; j < count; j++) {
                            String oneValue = st.nextToken();
                            double value = Double.NaN;
                            if (oneValue.compareTo("NA") == 0) {
                                // Missing values not yet handled.
                            } else {
                                try {
                                    value = new Double(oneValue);
                                } catch (NumberFormatException e) {
                                    throw new RuntimeException(e.getMessage());
                                }
                            }
                            traitParam.setParameterValue(j, value);
                        }
                    }
                }

                // Find missing values
                double[] allValues = traitParameter.getParameterValues();
                missingIndices = new ArrayList<Integer>();
                for (int i = 0; i < allValues.length; i++) {
                    if ((new Double(allValues[i])).isNaN()) {
                        traitParameter.setParameterValue(i, 0);
                        missingIndices.add(i);
                    }
                }

                if (xo.hasChildNamed(MISSING)) {
                    XMLObject cxo = (XMLObject) xo.getChild(MISSING);
                    Parameter missingParameter = new Parameter.Default(allValues.length, 0.0);
                    for (int i : missingIndices) {
                        missingParameter.setParameterValue(i, 1.0);
                    }
                    missingParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, allValues.length));
/*					CompoundParameter missingParameter = new CompoundParameter(MISSING);
					System.err.println("TRAIT: "+traitParameter.toString());
					System.err.println("CNT:   "+traitParameter.getNumberOfParameters());
					for(int i : missingIndices) {
						Parameter thisParameter = traitParameter.getIndicatorParameter(i);
						missingParameter.addParameter(thisParameter);
					}*/
                    ParameterParser.replaceParameter(cxo, missingParameter);
                }

            }

//            Parameter traits = null;
//            Parameter check = null;
            Model samplingDensity = null;

//            if (xo.hasChildNamed(SAMPLING_DENSITY)) {
//                XMLObject cxo = (XMLObject) xo.getChild(SAMPLING_DENSITY);
//                samplingDensity = (Model) cxo.getChild(Model.class);
//            }
//            if (xo.hasChildNamed(RANDOMIZE)) {
//                XMLObject cxo = (XMLObject) xo.getChild(RANDOMIZE);
//                traits = (Parameter) cxo.getChild(Parameter.class);
//            }

//            if (xo.hasChildNamed(CHECK)) {
//                XMLObject cxo = (XMLObject) xo.getChild(CHECK);
//                check = (Parameter) cxo.getChild(Parameter.class);
//            }

            boolean useTreeLength = xo.getAttribute(USE_TREE_LENGTH, false);

            boolean scaleByTime = xo.getAttribute(SCALE_BY_TIME, false);

            boolean reportAsMultivariate = false;
            if (xo.hasAttribute(REPORT_MULTIVARIATE) && xo.getBooleanAttribute(REPORT_MULTIVARIATE))
                reportAsMultivariate = true;

            MultivariateDistributionLikelihood rootPrior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
            if (!(rootPrior.getDistribution() instanceof MultivariateDistribution))
                throw new XMLParseException("Only multivariate normal priors allowed for Gibbs sampling the root trait");

            MultivariateNormalDistribution rootDistribution =
                    (MultivariateNormalDistribution) rootPrior.getDistribution();
//            if (integrate)
            return new IntegratedMultivariateTraitLikelihood(traitName, treeModel, diffusionModel,
                    traitParameter, missingIndices, cacheBranches,
                    scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate, rootDistribution);

//            AbstractMultivariateTraitLikelihood like =
//                    new SampledMultivariateTraitLikelihood(traitName, treeModel, diffusionModel,
//                            traitParameter, missingIndices, cacheBranches,
//                            scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate);

//            if (traits != null) {
//                like.randomize(traits);
//            }
//
//            if (check != null) {
//                like.check(check);
//            }
//
//            return like;
        }


        private Parameter getTraitParameterByName(CompoundParameter traits, String name) {

            for (int i = 0; i < traits.getNumberOfParameters(); i++) {
                Parameter found = traits.getParameter(i);
                if (found.getStatisticName().compareTo(name) == 0)
                    return found;
            }
            return null;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of a continuous trait evolving on a tree by a " +
                    "given diffusion model.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
//                new ElementRule(TRAIT_PARAMETER, new XMLSyntaxRule[]{
//                        new ElementRule(Parameter.class)
//                }),
//                AttributeRule.newBooleanRule(INTEGRATE,true),
                new ElementRule(MultivariateDistributionLikelihood.class),
                new ElementRule(MultivariateDiffusionModel.class),
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class, true),
                AttributeRule.newDoubleArrayRule("cut", true),
                AttributeRule.newBooleanRule(REPORT_MULTIVARIATE, true),
                AttributeRule.newBooleanRule(USE_TREE_LENGTH, true),
                AttributeRule.newBooleanRule(SCALE_BY_TIME, true),
//                new ElementRule(Parameter.class, true),
//                new ElementRule(RANDOMIZE, new XMLSyntaxRule[]{
//                        new ElementRule(Parameter.class)
//                }, true),
//                new ElementRule(CHECK, new XMLSyntaxRule[]{
//                        new ElementRule(Parameter.class)
//                }, true)
        };


        public Class getReturnType() {
            return AbstractMultivariateTraitLikelihood.class;
        }
    };

    private double[] meanCache;
    private double[] upperPrecisionCache;
    private double[] lowerPrecisionCache;
    private double[] logRemainderDensityCache;

    private int dimData;
    private int dimTrait;
    private int dim;

    private double[] rootPriorMean;
    private double[][] rootPriorPrecision;
    private double rootPriorPrecisionDeterminant;

}
