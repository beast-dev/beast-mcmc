package dr.evomodel.continuous;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.xml.*;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.*;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */


public abstract class AbstractMultivariateTraitLikelihood extends AbstractModelLikelihood implements NodeAttributeProvider {

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
    public static final String RANDOMIZE_LOWER = "lower";
    public static final String RANDOMIZE_UPPER = "upper";
    public static final String CHECK = "check";
    public static final String USE_TREE_LENGTH = "useTreeLength";
    public static final String SCALE_BY_TIME = "scaleByTime";
    public static final String SUBSTITUTIONS = "substitutions";
    public static final String SAMPLING_DENSITY = "samplingDensity";
    public static final String INTEGRATE = "integrateInternalTraits";
    public static final String JITTER = "jitter";
    public static final String WINDOW = "window";
    public static final String DUPLICATES = "duplicatesOnly";

    public AbstractMultivariateTraitLikelihood(String traitName,
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
            addModel(samplingDensity);
        }

        if (traitParameter != null)
            addVariable(traitParameter);

        this.reportAsMultivariate = reportAsMultivariate;

        this.cacheBranches = cacheBranches;
        if (cacheBranches) {
            cachedLogLikelihoods = new double[treeModel.getNodeCount()];
            storedCachedLogLikelihood = new double[treeModel.getNodeCount()];
            validLogLikelihoods = new boolean[treeModel.getNodeCount()];
            storedValidLogLikelihoods = new boolean[treeModel.getNodeCount()];            
        }
		
        this.scaleByTime = scaleByTime;
        this.useTreeLength = useTreeLength;

        StringBuffer sb = new StringBuffer("Creating multivariate diffusion model:\n");
        sb.append("\tTrait: ").append(traitName).append("\n");
        sb.append("\tDiffusion process: ").append(diffusionModel.getId()).append("\n");
        sb.append("\tHeterogenity model: ").append(rateModel != null ? rateModel.getId() : "homogeneous").append("\n");
        sb.append("\tTree normalization: ").append(scaleByTime ? (useTreeLength ? "length" : "height") : "off").append("\n");
        if (scaleByTime) {
            recalculateTreeLength();
            if (useTreeLength) {
                sb.append("\tInitial tree length: ").append(treeLength).append("\n");
            } else {
                sb.append("\tInitial tree height: ").append(treeLength).append("\n");
            }
        }
        sb.append(extraInfo());
        sb.append("\tPlease cite Suchard, Lemey and Rambaut (in preparation) if you publish results using this model.");

        Logger.getLogger("dr.evomodel").info(sb.toString());

        recalculateTreeLength();

    }

    protected abstract String extraInfo();

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

        if (!cacheBranches) {
            likelihoodKnown = false;
            if (model == treeModel)
                recalculateTreeLength();
            return;
        }

        if (model == diffusionModel) {
            updateAllNodes();
        }

        // fireTreeEvents sends two events here when a node trait is changed,
        // ignoring object instance Parameter case

        else if (model == treeModel) {
            if (object instanceof TreeModel.TreeChangedEvent) {
                TreeModel.TreeChangedEvent event = (TreeModel.TreeChangedEvent) object;
                if (event.isHeightChanged()) {
                    recalculateTreeLength();
                    if (useTreeLength || (scaleByTime && treeModel.isRoot(event.getNode())))
                        updateAllNodes();
                    else {
                        updateNodeAndChildren(event.getNode());
                    }
                } else if (event.isNodeParameterChanged()) {
                    updateNodeAndChildren(event.getNode());
                } else if (event.isNodeChanged()) {
                    recalculateTreeLength();
                    if (useTreeLength || (scaleByTime && treeModel.isRoot(event.getNode())))
                        updateAllNodes();
                    else {
                        updateNodeAndChildren(event.getNode());
                    }
                } else {
                    throw new RuntimeException("Unexpected TreeModel TreeChangedEvent occuring in AbstractMultivariateTraitLikelihood");
                }
            } else if (object instanceof Parameter) {
                // Ignoring                
            } else {
                throw new RuntimeException("Unexpected TreeModel event occuring in AbstractMultivariateTraitLikelihood");
            }
        } else if (model == rateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                if (((Parameter)object).getDimension() == 2*(treeModel.getNodeCount()-1))
                    updateNode(treeModel.getNode(index)); // This is a branch specific update
                else
                    updateAllNodes(); // Probably an epoch model
            }
        } else {
            throw new RuntimeException("Unknown componentChangedEvent");
        }
    }

    private void updateAllNodes() {
        for(int i=0; i<treeModel.getNodeCount(); i++)
            validLogLikelihoods[i] = false;
        likelihoodKnown = false;
    }

    private void updateNode(NodeRef node) {
        validLogLikelihoods[node.getNumber()] = false;
        likelihoodKnown = false;
    }

    private void updateNodeAndChildren(NodeRef node) {
        validLogLikelihoods[node.getNumber()] = false;
        for(int i=0; i<treeModel.getChildCount(node); i++)
            validLogLikelihoods[treeModel.getChild(node,i).getNumber()] = false;
        likelihoodKnown = false;
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
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        // All parameter changes are handled first by the treeModel
        if (!cacheBranches)
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

        if (cacheBranches) {
            System.arraycopy(cachedLogLikelihoods,0,storedCachedLogLikelihood,0,treeModel.getNodeCount());
            System.arraycopy(validLogLikelihoods,0,storedValidLogLikelihoods,0,treeModel.getNodeCount());
        }
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        treeLength = storedTreeLength;

        if (cacheBranches) {
            double[] tmp = storedCachedLogLikelihood;
            storedCachedLogLikelihood = cachedLogLikelihoods;
            cachedLogLikelihoods = tmp;
            boolean[] tmp2 = storedValidLogLikelihoods;
            storedValidLogLikelihoods = validLogLikelihoods;
            validLogLikelihoods = tmp2;
        }
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

    public abstract double getLogDataLikelihood();

    public void makeDirty() {
        likelihoodKnown = false;
        if (cacheBranches)
            updateAllNodes();
    }

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId()+".joint"),
                new NumberColumn(getId()+".data") {
                    public double getDoubleValue() {
                        return getLogDataLikelihood();
                    }
                }
        };
    }

    public abstract double calculateLogLikelihood();

    public double getMaxLogLikelihood() {
        return maxLogLikelihood;
    }


    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    private String[] attributeLabel = null;

    public String[] getNodeAttributeLabel() {
        if (attributeLabel == null) {
            double[] trait = getRootNodeTrait();
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

    protected double[] getRootNodeTrait() {
        return treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName);
    }

    protected abstract double[] traitForNode(TreeModel tree, NodeRef node, String traitName);

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        double trait[] = traitForNode(treeModel, node, traitName);
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

    public void randomize(Parameter trait, double[] lower, double[] upper) {
        // Draws each dimension in each trait from U[lower, upper)
        for(int i = 0; i < trait.getDimension(); i++) {
            final int whichLower = i % lower.length;
            final int whichUpper = i % upper.length;
            final double newValue = MathUtils.uniform(lower[whichLower],upper[whichUpper]);
            trait.setParameterValue(i, newValue);    
        }
        //diffusionModel.randomize(trait);
    }

    public void jitter(Parameter trait, int dim, double[] window, boolean duplicates, boolean verbose) {
        int numTraits = trait.getDimension() / dim;
        boolean[] update = new boolean[numTraits];
        if (!duplicates) {
            Arrays.fill(update, true);
        } else {
            DoubleArray[] traitArray = new DoubleArray[numTraits];
            for (int i = 0; i < numTraits; i++) {
                double[] x = new double[dim];
                for (int j = 0; j < dim; j++) {
                    x[j] = trait.getParameterValue(i * dim + j);
                }
                traitArray[i] = new DoubleArray(x,i);
            }
            Arrays.sort(traitArray);
            // Mark duplicates
            for (int i = 1; i < numTraits; i++) {
                if (traitArray[i].compareTo(traitArray[i-1]) == 0) {
                    update[traitArray[i-1].getIndex()] = true;
                    update[traitArray[i].getIndex()] = true;
                }
            }
        }
        for (int i = 0; i < numTraits; i++) {
            if (update[i]) {
                StringBuffer sb1 = null;
                StringBuffer sb2 = null;
                if (verbose) {
                    sb1 = new StringBuffer();
                    sb2 = new StringBuffer();
                }
                for (int j = 0; j < dim; j++) {
                    final double oldValue = trait.getParameterValue(i * dim + j);
                    final double newValue = window[j % window.length] * (MathUtils.nextDouble() - 0.5) +
                            oldValue;
                    trait.setParameterValue(i * dim + j, newValue);
                    if (verbose) {
                        sb1.append(" ").append(oldValue);
                        sb2.append(" ").append(newValue);
                    }
                }
                if (verbose) {
                    Logger.getLogger("dr.evomodel.continuous").info(
                            "  Replacing trait #"+(i+1)+"  Old:"+sb1.toString()+" New: "+sb2.toString()
                    );
                }
            }
        }        
    }

    class DoubleArray implements Comparable {

        double[] value;
        int index;

        DoubleArray(double[] value, int index) {
            this.value = value;
            this.index = index;
        }

        public double[] getValues() {
            return value;
        }

        public int getIndex() {
            return index;
        }

        public int compareTo(Object o) {
            double[] x = ((DoubleArray) o).getValues();
            for(int i = 0; i < value.length; i++) {
                if (value[i] > x[i]) {
                    return 1;
                } else if (value[i] < x[i]) {
                    return -1;
                }
            }
            return 0;
        }
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
            boolean integrate = xo.getAttribute(INTEGRATE,false);

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
                    XMLObject cxo = xo.getChild(MISSING);
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
						missingParameter.addVariable(thisParameter);
					}*/
                    ParameterParser.replaceParameter(cxo, missingParameter);
                }


            }

            Model samplingDensity = null;

            if (xo.hasChildNamed(SAMPLING_DENSITY)) {
                XMLObject cxo = xo.getChild(SAMPLING_DENSITY);
                samplingDensity = (Model) cxo.getChild(Model.class);
            }

            boolean useTreeLength = xo.getAttribute(USE_TREE_LENGTH, false);

            boolean scaleByTime = xo.getAttribute(SCALE_BY_TIME, false);

            boolean reportAsMultivariate = false;
            if (xo.hasAttribute(REPORT_MULTIVARIATE) && xo.getBooleanAttribute(REPORT_MULTIVARIATE))
                reportAsMultivariate = true;

            if (integrate) {

                MultivariateDistributionLikelihood rootPrior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
                if (rootPrior == null || !(rootPrior.getDistribution() instanceof MultivariateDistribution))
                    throw new XMLParseException("Only multivariate normal priors allowed for Gibbs sampling the root trait");

                MultivariateNormalDistribution rootDistribution =
                    (MultivariateNormalDistribution) rootPrior.getDistribution();

                return new IntegratedMultivariateTraitLikelihood(traitName, treeModel, diffusionModel,
                            traitParameter, missingIndices, cacheBranches,
                            scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate, rootDistribution);
            }

            AbstractMultivariateTraitLikelihood like =
                    new SampledMultivariateTraitLikelihood(traitName, treeModel, diffusionModel,
                            traitParameter, missingIndices, cacheBranches,
                            scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate);

            if (xo.hasChildNamed(RANDOMIZE)) {
                XMLObject cxo = xo.getChild(RANDOMIZE);
                Parameter traits = (Parameter) cxo.getChild(Parameter.class);
                double[] randomizeLower;
                double[] randomizeUpper;
                if (cxo.hasAttribute(RANDOMIZE_LOWER)) {
                    randomizeLower = cxo.getDoubleArrayAttribute(RANDOMIZE_LOWER);
                } else {
                    randomizeLower = new double[] { -90.0 };
                }
                if (cxo.hasAttribute(RANDOMIZE_UPPER)) {
                    randomizeUpper = cxo.getDoubleArrayAttribute(RANDOMIZE_UPPER);
                } else {
                    randomizeUpper = new double[] { +90.0 };
                }
                like.randomize(traits, randomizeLower, randomizeUpper);
            }

            if (xo.hasChildNamed(JITTER)) {
                XMLObject cxo = xo.getChild(JITTER);
                Parameter traits = (Parameter) cxo.getChild(Parameter.class);
                double[] window = cxo.getDoubleArrayAttribute(WINDOW); // Must be included, no default value
                boolean duplicates = cxo.getAttribute(DUPLICATES,true); // default = true
                like.jitter(traits, diffusionModel.getPrecisionmatrix().length, window, duplicates, true);
            }

            if (xo.hasChildNamed(CHECK)) {
                XMLObject cxo = xo.getChild(CHECK);
                Parameter check = (Parameter) cxo.getChild(Parameter.class);
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
                AttributeRule.newBooleanRule(INTEGRATE,true),
                new ElementRule(MultivariateDistributionLikelihood.class,true),
                new ElementRule(MultivariateDiffusionModel.class),
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class, true),
                AttributeRule.newDoubleArrayRule("cut", true),
                AttributeRule.newBooleanRule(REPORT_MULTIVARIATE, true),
                AttributeRule.newBooleanRule(USE_TREE_LENGTH, true),
                AttributeRule.newBooleanRule(SCALE_BY_TIME, true),
                new ElementRule(Parameter.class, true),
                new ElementRule(RANDOMIZE, new XMLSyntaxRule[]{
                        AttributeRule.newDoubleRule(RANDOMIZE_LOWER,true),
                        AttributeRule.newDoubleRule(RANDOMIZE_UPPER,true),
                        new ElementRule(Parameter.class)
                }, true),
                new ElementRule(JITTER, new XMLSyntaxRule[] {
                        AttributeRule.newDoubleArrayRule(WINDOW),
                        AttributeRule.newBooleanRule(DUPLICATES, true),
                        new ElementRule(Parameter.class),

                }, true),
                new ElementRule(CHECK, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true)
        };


        public Class getReturnType() {
            return AbstractMultivariateTraitLikelihood.class;
        }
    };

    TreeModel treeModel = null;
    MultivariateDiffusionModel diffusionModel = null;
    String traitName = null;
    CompoundParameter traitParameter;
    List<Integer> missingIndices;

    ArrayList dataList = new ArrayList();

    protected double logLikelihood;
    protected double maxLogLikelihood = Double.NEGATIVE_INFINITY;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private BranchRateModel rateModel = null;
    private boolean hasRateModel = false;

    private double treeLength;
    private double storedTreeLength;

    private final boolean reportAsMultivariate;

    private final boolean scaleByTime;
    private final boolean useTreeLength;

    protected boolean cacheBranches;
    protected double[] cachedLogLikelihoods;
    protected double[] storedCachedLogLikelihood;
    protected boolean[] validLogLikelihoods;
    protected boolean[] storedValidLogLikelihoods;
}

