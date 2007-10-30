package dr.evomodel.continuous;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Marc Suchard
 */


public class MultivariateTraitLikelihood extends AbstractModel implements Likelihood, NodeAttributeProvider {

	public static final String TRAIT_LIKELIHOOD = "multivariateTraitLikelihood";
	public static final String TRAIT_NAME = "traitName";
	public static final String ROOT_PRIOR = "rootPrior";
	public static final String MODEL = "diffusionModel";
	public static final String TREE = "tree";
	public static final String TRAIT_PARAMETER = "traitParameter";
	public static final String SET_TRAIT = "setOutcomes";
	public static final String MISSING = "missingIndicator";
	public static final String CACHE_BRANCHES = "cacheBranches";

	public MultivariateTraitLikelihood(String traitName,
	                                   TreeModel treeModel,
	                                   MultivariateDiffusionModel diffusionModel,
	                                   CompoundParameter traitParameter,
	                                   List<Integer> missingIndices,
	                                   boolean cacheBranches) {

		super(TRAIT_LIKELIHOOD);

		this.traitName = traitName;
		this.treeModel = treeModel;
		this.diffusionModel = diffusionModel;
		this.traitParameter = traitParameter;
		this.missingIndices = missingIndices;
//		this.jeffreysPrior = jeffreysPrior;
		addModel(treeModel);
		addModel(diffusionModel);
		addParameter(traitParameter);

		if (cacheBranches)
			cachedLikelihoods = new HashMap<NodeRef, Double>();
//		System.err.println("MADE IT");

//		this.traitName = traitName;
	}

	// **************************************************************
	// ModelListener IMPLEMENTATION
	// **************************************************************

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// todo should cache branch-specific calculations
		likelihoodKnown = false;
	}

	// **************************************************************
	// ParameterListener IMPLEMENTATION
	// **************************************************************

	protected void handleParameterChangedEvent(Parameter parameter, int index) {

		// Should recalculate the matrix-functions (determinant) of diffusion once
		// Maybe do this in the MultivariateDiffusionModel class

		// todo should cache branch-specific calculations

		if (cachedLikelihoods != null) {
			String paramName = parameter.getDimensionName(index);
			if (paramName.startsWith("root")) {
				System.err.println("ROOT");
//				treeModel.setNodeAttribute();
			} else if (paramName.startsWith("node")) {
				System.err.println("INTERNAL " + paramName);
//				NodeRef node = treeModel.getNodeOfParameter()
			} else { // leaf
				System.err.println("LEAF " + paramName);

			}
		}
//		System.err.println("P changed = "+parameter.getDimensionName(index)+" "+index);

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
	}

	/**
	 * Restores the precalculated state: that is the intervals of the tree.
	 */
	protected void restoreState() {
		likelihoodKnown = storedLikelihoodKnown;
		logLikelihood = storedLogLikelihood;
	}

	protected void acceptState() {
	} // nothing to do

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
//		System.err.println("TRAITLIKELIHOOD CALCULATION CALLED " + (likelihoodKnown ? "yes" : "no"));
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

		double logLikelihood = traitLogLikelihood(treeModel.getRoot());
		if (logLikelihood > maxLogLikelihood) {
			maxLogLikelihood = logLikelihood;
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

			double time = (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node)) * treeModel.getNodeRate(node);
			if (cachedLikelihoods != null && cachedLikelihoods.containsKey(node)) {
				logL = cachedLikelihoods.get(node);
			} else {
				logL = diffusionModel.getLogLikelihood(parentTrait, childTrait, time);
				if (new Double(logL).isNaN()) {
					System.err.println("time = " + time);
					System.err.println(new Matrix(diffusionModel.getPrecisionmatrix()));
				}
				if (cachedLikelihoods != null)
					cachedLikelihoods.put(node, logL);
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

			MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			CompoundParameter traitParameter = (CompoundParameter) xo.getSocketChild(TRAIT_PARAMETER);

			boolean cacheBranches = false;
			if (xo.hasAttribute(CACHE_BRANCHES))
				cacheBranches = xo.getBooleanAttribute(CACHE_BRANCHES);

			List<Integer> missingIndices = null;
			String traitName = "trait";

			if (xo.hasAttribute(TRAIT_NAME)) {

				traitName = xo.getStringAttribute(TRAIT_NAME);

				// Fill in attributeValues
				int taxonCount = treeModel.getTaxonCount();
				for (int i = 0; i < taxonCount; i++) {
					String taxonName = treeModel.getTaxonId(i);
					String paramName = taxonName + ".trait";
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
									value = (new Double(oneValue)).doubleValue();
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

				if (xo.hasSocket(MISSING)) {
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
						Parameter thisParameter = traitParameter.getParameter(i);
						missingParameter.addParameter(thisParameter);
					}*/
					replaceParameter(cxo, missingParameter);
				}


			}
			return new MultivariateTraitLikelihood(traitName, treeModel, diffusionModel,
					traitParameter, missingIndices, cacheBranches);
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
				new StringAttributeRule(TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
//			AttributeRule.newBooleanRule(JEFFERYS_PRIOR),
				new ElementRule(MultivariateDiffusionModel.class),
				new ElementRule(TreeModel.class)


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
}

