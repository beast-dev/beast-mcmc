package dr.evomodel.continuous;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */


public class MultivariateTraitLikelihood extends AbstractModelLikelihood implements NodeAttributeProvider {

    public static final String TRAIT_LIKELIHOOD = "multivariateTraitLikelihood";
    public static final String TRAIT_NAME = "traitName";
    public static final String ROOT_PRIOR = "rootPrior";
    public static final String MODEL = "diffusionModel";
    public static final String TREE = "tree";
    public static final String TRAIT_PARAMETER = "traitParameter";
    public static final String SET_TRAIT = "setOutcomes";
    public static final String MISSING = "missingIndicator";
    public static final String CACHE_BRANCHES = "cacheBranches";
    public static final String REPORT_MULTIVARIATE = "reportAsMultivariate";
    public static final String DEFAULT_TRAIT_NAME = "trait";
    public static final String RANDOMIZE = "randomize";
    public static final String CHECK = "check";
    public static final String USE_TREE_LENGTH = "useTreeLength";
    public static final String SCALE_BY_TIME = "scaleByTime";
    public static final String SUBSTITUTIONS = "substitutions";
    public static final String SAMPLING_DENSITY = "samplingDensity";

    public MultivariateTraitLikelihood(String traitName,
                                       TreeModel treeModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       CompoundParameter traitParameter,
                                       List<Integer> missingIndices,
                                       boolean cacheBranches,
                                       boolean scaleByTime,
                                       boolean useTreeLength,
                                       BranchRateModel rateModel,
                                       Model samplingDensity,
                                       boolean reportAsMultivariate) {

        super(TRAIT_LIKELIHOOD);

        this.traitName = traitName;
        this.treeModel = treeModel;
        this.rateModel = rateModel;
        this.diffusionModel = diffusionModel;
        this.traitParameter = traitParameter;
        this.missingIndices = missingIndices;
        addModel(treeModel);
        addModel(diffusionModel);

        if (rateModel != null) {
            hasRateModel = true;
            addModel(rateModel);
        }

        if (samplingDensity != null) {
            Model samplingDensity1 = samplingDensity;
            addModel(samplingDensity);
        }
        addParameter(traitParameter);

        this.reportAsMultivariate = reportAsMultivariate;

//		if (cacheBranches)
//			cachedLikelihoods = new HashMap<NodeRef, Double>();


        this.scaleByTime = scaleByTime;
        this.useTreeLength = useTreeLength;

        StringBuffer sb = new StringBuffer("Creating multivariate diffusion model:\n");
        sb.append("\tTrait: " + traitName + "\n");
        sb.append("\tDiffusion process: " + diffusionModel.getId() + "\n");
        sb.append("\tHeterogenity model: " + (hasRateModel ? rateModel.getId() : "homogeneous") + "\n");
//        if (!hasRateModel) {
        sb.append("\tTree normalization: " + (scaleByTime ? (useTreeLength ? "length" : "height") : "off") + "\n");
        if (scaleByTime) {
            recalculateTreeLength();
            if (useTreeLength) {
                sb.append("\tInitial tree length: " + treeLength + "\n");
            } else {
                sb.append("\tInitial tree height: " + treeLength + "\n");
            }
//            }
        }
        sb.append("\tPlease cite Suchard, Lemey and Rambaut (in preparation) if you publish results using this model.");

        Logger.getLogger("dr.evomodel").info(sb.toString());

        recalculateTreeLength();

    }

    public String getTraitName() {
        return traitName;
    }

    public double getRescaledBranchLength(NodeRef node) {

        double length = treeModel.getBranchLength(node);

        if (hasRateModel)
            length *= rateModel.getBranchRate(treeModel, node);

        if (scaleByTime)
            return length / treeLength;

        return length;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        likelihoodKnown = false;
        if (model == treeModel)
            recalculateTreeLength();
    }


    public void recalculateTreeLength() {

        if (!scaleByTime)
            return;

        if (useTreeLength) {
            treeLength = 0;
            for (int i = 0; i < treeModel.getNodeCount(); i++) {
                NodeRef node = treeModel.getNode(i);
                if (!treeModel.isRoot(node))
                    treeLength += treeModel.getBranchLength(node); // Bug was here
            }
        } else { // Normalizing by tree height.
            treeLength = treeModel.getNodeHeight(treeModel.getRoot());
        }

    }

    // **************************************************************
    // ParameterListener IMPLEMENTATION
    // **************************************************************

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {

        likelihoodKnown = false;

    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
        storedTreeLength = treeLength;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        treeLength = storedTreeLength;
    }

    protected void acceptState() {
    } // nothing to do

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public MultivariateDiffusionModel getDiffusionModel() {
        return diffusionModel;
    }

//	public boolean getInSubstitutionTime() {
//		return inSubstitutionTime;
//	}

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public Model getModel() {
        return this;
    }

    public String toString() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";

    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        double logLikelihood = traitLogLikelihood(null, treeModel.getRoot());
        if (logLikelihood > maxLogLikelihood) {
            maxLogLikelihood = logLikelihood;
        }
        return logLikelihood;
    }

    public double getMaxLogLikelihood() {
        return maxLogLikelihood;
    }

    private double traitLogLikelihood(double[] parentTrait, NodeRef node) {

        double logL = 0.0;
        double[] childTrait = treeModel.getMultivariateNodeTrait(node, traitName);

        if (parentTrait != null) {
            double time = getRescaledBranchLength(node);

//			if (inSubstitutionTime) {
//				time *= treeModel.getNodeRate(node);
//			}

//			if (cachedLikelihoods != null && cachedLikelihoods.containsKey(node)) {
//				logL = cachedLikelihoods.get(node);
//			} else {
            logL = diffusionModel.getLogLikelihood(parentTrait, childTrait, time);
            if (new Double(logL).isNaN()) {
                System.err.println("MultivariateTraitLikelihood: likelihood is undefined");
                System.err.println("time = " + time);
                System.err.println("parent trait value = " + new Vector(parentTrait));
                System.err.println("child trait value = " + new Vector(childTrait));
                System.err.println("precision matrix = " + new Matrix(diffusionModel.getPrecisionmatrix()));
                if (diffusionModel.getPrecisionParameter() instanceof CompoundSymmetricMatrix) {
                    CompoundSymmetricMatrix csMatrix = (CompoundSymmetricMatrix) diffusionModel.getPrecisionParameter();
                    System.err.println("diagonals = " + new Vector(csMatrix.getDiagonals()));
                    System.err.println("off diagonal = " + csMatrix.getOffDiagonal());
                }
            }
//				if (cachedLikelihoods != null) {
//					cachedLikelihoods.put(node, logL);
//				}
//			}
        }
        int childCount = treeModel.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            logL += traitLogLikelihood(childTrait, treeModel.getChild(node, i));
        }

        if (new Double(logL).isNaN()) {
            System.err.println("logL = " + logL);
            System.err.println(new Matrix(diffusionModel.getPrecisionmatrix()));
            System.exit(-1);
        }

        return logL;
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    private String[] attributeLabel = null;

    public String[] getNodeAttributeLabel() {
        if (attributeLabel == null) {
            double[] trait = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName);
            if (trait.length == 1 || reportAsMultivariate)
                attributeLabel = new String[]{traitName};
            else {
                attributeLabel = new String[trait.length];
                for (int i = 1; i <= trait.length; i++)
                    attributeLabel[i - 1] = traitName + i;
            }
        }
        return attributeLabel;
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        double trait[] = treeModel.getMultivariateNodeTrait(node, traitName);
        String[] value;
        if (!reportAsMultivariate || trait.length == 1) {
            value = new String[trait.length];
            for (int i = 0; i < trait.length; i++)
                value[i] = Double.toString(trait[i]);
        } else {
            StringBuffer sb = new StringBuffer("{");
            for (int i = 0; i < trait.length - 1; i++)
                sb.append(Double.toString(trait[i])).append(",");
            sb.append(Double.toString(trait[trait.length - 1])).append("}");
            value = new String[]{sb.toString()};
        }
        return value;
    }

    private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }


    public void randomize(Parameter trait) {
        diffusionModel.randomize(trait);
    }

    public void check(Parameter trait) throws XMLParseException {
        diffusionModel.check(trait);
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
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

            boolean cacheBranches = xo.getAttribute(CACHE_BRANCHES, false);

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
                    replaceParameter(cxo, missingParameter);
                }


            }

            Parameter traits = null;
            Parameter check = null;
            Model samplingDensity = null;

            if (xo.hasChildNamed(SAMPLING_DENSITY)) {
                XMLObject cxo = (XMLObject) xo.getChild(SAMPLING_DENSITY);
                samplingDensity = (Model) cxo.getChild(Model.class);
            }
            if (xo.hasChildNamed(RANDOMIZE)) {
                XMLObject cxo = (XMLObject) xo.getChild(RANDOMIZE);
                traits = (Parameter) cxo.getChild(Parameter.class);
            }

            if (xo.hasChildNamed(CHECK)) {
                XMLObject cxo = (XMLObject) xo.getChild(CHECK);
                check = (Parameter) cxo.getChild(Parameter.class);
            }

            boolean useTreeLength = xo.getAttribute(USE_TREE_LENGTH, false);

            boolean scaleByTime = xo.getAttribute(SCALE_BY_TIME, false);

            boolean reportAsMultivariate = false;
            if (xo.hasAttribute(REPORT_MULTIVARIATE) && xo.getBooleanAttribute(REPORT_MULTIVARIATE))
                reportAsMultivariate = true;

            MultivariateTraitLikelihood like =
                    new MultivariateTraitLikelihood(traitName, treeModel, diffusionModel,
                            traitParameter, missingIndices, cacheBranches,
                            scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate);

            if (traits != null) {
                like.randomize(traits);
            }

            if (check != null) {
                like.check(check);
            }

            return like;
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
                new ElementRule(TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(MultivariateDiffusionModel.class),
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class, true),
                AttributeRule.newDoubleArrayRule("cut", true),
                AttributeRule.newBooleanRule(REPORT_MULTIVARIATE, true),
                AttributeRule.newBooleanRule(USE_TREE_LENGTH, true),
                AttributeRule.newBooleanRule(SCALE_BY_TIME, true),
                new ElementRule(Parameter.class, true),
                new ElementRule(RANDOMIZE, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true),
                new ElementRule(CHECK, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true)
        };


        public Class getReturnType() {
            return MultivariateTraitLikelihood.class;
        }
    };

    TreeModel treeModel = null;
    MultivariateDiffusionModel diffusionModel = null;
    String traitName = null;
    CompoundParameter traitParameter;
    List<Integer> missingIndices;

    ArrayList dataList = new ArrayList();

    private double logLikelihood;
    private double maxLogLikelihood = Double.NEGATIVE_INFINITY;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private BranchRateModel rateModel = null;
    private boolean hasRateModel = false;

    private double treeLength;
    private double storedTreeLength;

    private final boolean reportAsMultivariate;

    private final boolean scaleByTime;
    private final boolean useTreeLength;
}

