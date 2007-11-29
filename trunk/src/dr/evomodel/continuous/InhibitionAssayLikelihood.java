package dr.evomodel.continuous;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */


public class InhibitionAssayLikelihood extends AbstractModel implements Likelihood, NodeAttributeProvider {

    public static final String TRAIT_LIKELIHOOD = "inhibitionLikelihood";
    public static final String TRAIT_NAME = "traitName";
    public static final String ROOT_PRIOR = "rootPrior";
    public static final String MODEL = "diffusionModel";
    public static final String TREE = "tree";
    public static final String TRAIT_PARAMETER = "traitParameter";
    public static final String SET_TRAIT = "setOutcomes";
    public static final String MISSING = "missingIndicator";
    public static final String CACHE_BRANCHES = "cacheBranches";
    public static final String IN_REAL_TIME = "inRealTime";
    public static final String PRECISION = "precision";

    public InhibitionAssayLikelihood(TreeModel treeModel,

//	                                   List<Integer> missingIndices,
MatrixParameter dataParameter,
Parameter precision) {

        super(TRAIT_LIKELIHOOD);
        this.treeModel = treeModel;
        this.dataParameter = dataParameter;

        addModel(treeModel);

        addParameter(dataParameter);
        addParameter(precision);
        StringBuffer sb = new StringBuffer("Creating inhibition assay model:\n");
        sb.append("\tPlease cite O'Brien and Suchard (in preparation) if you publish results using this model.");

        Logger.getLogger("dr.evomodel").info(sb.toString());

    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        likelihoodKnown = false;

    }


    private double getBranchMean(TreeModel tree, NodeRef node) {

        if (tree.isRoot(node)) {
            return tree.getNodeTrait(node, "mean");
        } else {
            double rate;
            if (isClusterChangeOnBranchAbove(tree, node)) {
                rate = tree.getNodeTrait(node, "mean");

            } else {
                rate = getBranchMean(tree, tree.getParent(node));
            }
            return rate;
        }
    }

    public final boolean isClusterChangeOnBranchAbove(TreeModel tree, NodeRef node) {
        return tree.getNodeTrait(node, "indicator") == 1;
    }


    public final boolean areNodesInSameCluster(TreeModel tree, NodeRef node1, NodeRef node2) {
        return clusterStart(tree, node1) == clusterStart(tree, node2);
    }

    public final NodeRef clusterStart(TreeModel tree, NodeRef node) {
        if (tree.isRoot(node) || isClusterChangeOnBranchAbove(tree, node))
            return node;
        return clusterStart(tree, tree.getParent(node));
    }

    // **************************************************************
    // ParameterListener IMPLEMENTATION
    // **************************************************************

    protected void handleParameterChangedEvent(Parameter parameter, int index) {

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

//	public MultivariateDiffusionModel getDiffusionModel() {
//		return diffusionModel;
//	}
//
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

        logLikelihood = 0;

//        for(int i=0; i<treeModel.getNodeCount(); i++) {
//            treeModel.setNodeTrait(treeModel.getNode(i),"indicator",0.0);
//        }
//
//        treeModel.setNodeTrait(treeModel.getInternalNode(2),"indicator",1.0);
//
//
//        treeModel.setNodeTrait(treeModel.getRoot(),"mean",2.0);
//
//
//        for(int i=0; i<treeModel.getExternalNodeCount(); i++) {
//            System.err.println("mean for tip "+i+" = "+getBranchMean(treeModel,treeModel.getExternalNode(i)));
//        }


        final int N = treeModel.getExternalNodeCount();
        double[] mean = new double[N];

        for (int i = 0; i < N; i++)
            mean[i] = getBranchMean(treeModel, treeModel.getExternalNode(i));

        final double[][] data = dataParameter.getParameterAsMatrix();

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                if (i != j) {
                    NodeRef nodeI = treeModel.getExternalNode(i);
                    NodeRef nodeJ = treeModel.getExternalNode(j);
                    if (areNodesInSameCluster(treeModel, nodeI, nodeJ)) {
                        // todo do something to logLikelihood
                    } else {
                        // todo do something else to logLikelihood   
                    }
                }
            }
        }


        return logLikelihood;
    }

    public double getMaxLogLikelihood() {
        return maxLogLikelihood;
    }

    private double traitLogLikelihood(NodeRef node) {

        double logL = 0.0;
        double[] childTrait = treeModel.getMultivariateNodeTrait(node, "trait");
        // todo need to read in trait name in XML

        if (!treeModel.isRoot(node)) {
            NodeRef parent = treeModel.getParent(node);
            double[] parentTrait = treeModel.getMultivariateNodeTrait(parent, "trait");
            // todo need to read in trait name in XML

            double time = (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node));
            time /= treeLength;

            if (inSubstitutionTime) {
                time *= treeModel.getNodeRate(node);
            }

            if (cachedLikelihoods != null && cachedLikelihoods.containsKey(node)) {
                logL = cachedLikelihoods.get(node);
            } else {
                logL = diffusionModel.getLogLikelihood(parentTrait, childTrait, time);
                if (new Double(logL).isNaN()) {
                    System.err.println("time = " + time);
                    System.err.println("parent = " + new Vector(parentTrait));
                    System.err.println("child = " + new Vector(childTrait));
                    System.err.println(new Matrix(diffusionModel.getPrecisionmatrix()));
                }
                if (cachedLikelihoods != null) {
                    cachedLikelihoods.put(node, logL);
                }
            }
        }
        int childCount = treeModel.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            logL += traitLogLikelihood(treeModel.getChild(node, i));
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
            double[] trait = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), "trait");
            attributeLabel = new String[trait.length];
            if (trait.length == 1)
                attributeLabel[0] = traitName;
            else {
                for (int i = 1; i <= trait.length; i++)
                    attributeLabel[i - 1] = new String(traitName + i);
            }
        }
        return attributeLabel;
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        double trait[] = treeModel.getMultivariateNodeTrait(node, "trait");
//		StringBuffer sb = new StringBuffer();
//		sb.append("{");
//		for(int i=0; i<trait.length-1; i++) {
//			sb.append(trait[i]);
//			sb.append(",");
//		}
//		sb.append(trait[trait.length-1]);
//		sb.append("}");
        String[] value = new String[trait.length];
        for (int i = 0; i < trait.length; i++)
            value[i] = new Double(trait[i]).toString();

//		return new String[] {sb.toString()};  //To change body of implemented methods use File | Settings | File Templates.
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

//            System.err.println("did i get here?");

//            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
//			CompoundParameter traitParameter = (CompoundParameter) xo.getSocketChild(TRAIT_PARAMETER);

            MatrixParameter data = (MatrixParameter) xo.getChild(MatrixParameter.class);

            XMLObject cxo = (XMLObject) xo.getChild(PRECISION);

            Parameter precision = (Parameter) cxo.getChild(Parameter.class);

            int numTips = treeModel.getExternalNodeCount();
            if (numTips != data.getColumnDimension() || numTips != data.getRowDimension())
                throw new XMLParseException("Dimensions of matrix '" + data.getId() + "' do not match the number of taxa in '" + treeModel.getId() + "'");

//            boolean cacheBranches = false;
//			if (xo.hasAttribute(CACHE_BRANCHES))
//				cacheBranches = xo.getBooleanAttribute(CACHE_BRANCHES);
//
//			boolean inSubstitutionTime = false;
//			if (xo.hasAttribute(IN_REAL_TIME))
//				inSubstitutionTime = !xo.getBooleanAttribute(IN_REAL_TIME);
//
            List<Integer> missingIndices = null;
//			String traitName = "trait";

//			if (xo.hasAttribute(TRAIT_NAME)) {
//
//				traitName = xo.getStringAttribute(TRAIT_NAME);
//
//				// Fill in attributeValues
//				int taxonCount = treeModel.getTaxonCount();
//				for (int i = 0; i < taxonCount; i++) {
//					String taxonName = treeModel.getTaxonId(i);
//                    System.err.println("taxon "+i+" = "+taxonName);
//                    String paramName = taxonName + ".trait";
//					Parameter traitParam = getTraitParameterByName(traitParameter, paramName);
//					if (traitParam == null)
//						throw new RuntimeException("Missing trait parameters at tree tips");
//					String object = (String) treeModel.getTaxonAttribute(i, traitName);
//					if (object == null)
//						throw new RuntimeException("Trait \"" + traitName + "\" not found for taxa \"" + taxonName + "\"");
//					else {
//						StringTokenizer st = new StringTokenizer(object);
//						int count = st.countTokens();
//						if (count != traitParam.getDimension())
//							throw new RuntimeException("Trait length must match trait parameter dimension");
//						for (int j = 0; j < count; j++) {
//							String oneValue = st.nextToken();
//							double value = Double.NaN;
//							if (oneValue.compareTo("NA") == 0) {
//								// Missing values not yet handled.
//							} else {
//								try {
//									value = (new Double(oneValue)).doubleValue();
//								} catch (NumberFormatException e) {
//									throw new RuntimeException(e.getMessage());
//								}
//							}
//							traitParam.setParameterValue(j, value);
//						}
//					}
//				}

//				// Find missing values
//				double[] allValues = traitParameter.getParameterValues();
//				missingIndices = new ArrayList<Integer>();
//				for (int i = 0; i < allValues.length; i++) {
//					if ((new Double(allValues[i])).isNaN()) {
//						traitParameter.setParameterValue(i, 0);
//						missingIndices.add(i);
//					}
//				}
//
//				if (xo.hasSocket(MISSING)) {
//					XMLObject cxo = (XMLObject) xo.getChild(MISSING);
//					Parameter missingParameter = new Parameter.Default(allValues.length, 0.0);
//					for (int i : missingIndices) {
//						missingParameter.setParameterValue(i, 1.0);
//					}
//					missingParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, allValues.length));
///*					CompoundParameter missingParameter = new CompoundParameter(MISSING);
//					System.err.println("TRAIT: "+traitParameter.toString());
//					System.err.println("CNT:   "+traitParameter.getNumberOfParameters());
//					for(int i : missingIndices) {
//						Parameter thisParameter = traitParameter.getIndicatorParameter(i);
//						missingParameter.addParameter(thisParameter);
//					}*/
//					replaceParameter(cxo, missingParameter);
//				}

//			}
            return new InhibitionAssayLikelihood(treeModel, data, precision);
        }


        private Parameter getTraitParameterByName(CompoundParameter traits, String name) {
//			Parameter found = null;
//			System.err.println("LOOKING FOR: "+name);
            for (int i = 0; i < traits.getNumberOfParameters(); i++) {
                Parameter found = traits.getParameter(i);
//				System.err.println("COMPARE TO: "+found.getStatisticName());
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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//				new StringAttributeRule(TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
//				AttributeRule.newBooleanRule(IN_REAL_TIME, true),
//				new ElementRule(MultivariateDiffusionModel.class),
                new ElementRule(TreeModel.class),
                new ElementRule(MatrixParameter.class),
//                new ElementRule("precision", Parameter)
//                new ElementRule(Parameter.class)


                new ElementRule(PRECISION,
                                    new XMLSyntaxRule[]{
                                           new ElementRule(Parameter.class)  })


        };


        public Class getReturnType() {
            return MultivariateTraitLikelihood.class;
        }

        public void replaceParameter(XMLObject xo, Parameter newParam) throws XMLParseException {

            for (int i = 0; i < xo.getChildCount(); i++) {

                if (xo.getChild(i) instanceof Parameter) {

                    XMLObject rxo = null;
                    Object obj = xo.getRawChild(i);

                    if (obj instanceof Reference) {
                        rxo = ((Reference) obj).getReferenceObject();
                    } else if (obj instanceof XMLObject) {
                        rxo = (XMLObject) obj;
                    } else {
                        throw new XMLParseException("object reference not available");
                    }

                    if (rxo.getChildCount() > 0) {
                        throw new XMLParseException("No child elements allowed in parameter element.");
                    }

                    if (rxo.hasAttribute(XMLParser.IDREF)) {
                        throw new XMLParseException("References to " + xo.getName() + " parameters are not allowed in treeModel.");
                    }

                    if (rxo.hasAttribute(ParameterParser.VALUE)) {
                        throw new XMLParseException("Parameters in " + xo.getName() + " have values set automatically.");
                    }

                    if (rxo.hasAttribute(ParameterParser.UPPER)) {
                        throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                    }

                    if (rxo.hasAttribute(ParameterParser.LOWER)) {
                        throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                    }

                    if (rxo.hasAttribute(XMLParser.ID)) {

                        newParam.setId(rxo.getStringAttribute(XMLParser.ID));
                    }

                    rxo.setNativeObject(newParam);

                    return;
                }
            }
        }
    };

    TreeModel treeModel = null;
    MatrixParameter dataParameter = null;
    Parameter precision = null;
    MultivariateDiffusionModel diffusionModel = null;
    String traitName = null;
    //	private boolean jeffreysPrior = false;
    CompoundParameter traitParameter;
    List<Integer> missingIndices;

    ArrayList dataList = new ArrayList();

    private double logLikelihood;
    private double maxLogLikelihood = Double.NEGATIVE_INFINITY;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    //	private double[] cachedLikelihoods = null;
    private HashMap<NodeRef, Double> cachedLikelihoods = null;

    private double treeLength;
    private double storedTreeLength;

    private boolean inSubstitutionTime;
}

