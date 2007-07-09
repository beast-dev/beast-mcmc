package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 3, 2007
 * Time: 4:07:48 PM
 * To change this template use File | Settings | File Templates.
 */


public class MultivariateTraitLikelihood extends AbstractModel implements Likelihood {

	public static final String TRAIT_LIKELIHOOD = "multivariateTraitLikelihood";
	public static final String TRAIT_NAME = "traitName";
	public static final String ROOT_PRIOR = "rootPrior";
	public static final String MODEL = "diffusionModel";
	public static final String TREE = "tree";
	public static final String TRAIT_PARAMETER = "traitParameter";
	public static final String SET_TRAIT = "setOutcomes";
	public static final String MISSING = "missingIndicator";

	public MultivariateTraitLikelihood(TreeModel treeModel,
	                                   MultivariateDiffusionModel diffusionModel,
	                                   CompoundParameter traitParameter,
	                                   List<Integer> missingIndices) {

		super(TRAIT_LIKELIHOOD);

		this.treeModel = treeModel;
		this.diffusionModel = diffusionModel;
		this.traitParameter = traitParameter;
		this.missingIndices = missingIndices;
//		this.jeffreysPrior = jeffreysPrior;
		addModel(treeModel);
		addModel(diffusionModel);
		addParameter(traitParameter);

//		System.err.println("MADE IT");

//		this.traitName = traitName;
	}

	// **************************************************************
	// ModelListener IMPLEMENTATION
	// **************************************************************

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		likelihoodKnown = false;
	}

	// **************************************************************
	// ParameterListener IMPLEMENTATION
	// **************************************************************

	protected void handleParameterChangedEvent(Parameter parameter, int index) {

		// Should recalculate the matrix-functions (determinant) of diffusion once
		// Maybe do this in the MultivariateDiffusionModel class

	} // No parameters to respond to

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

	/**
	 * Adopt the state of the model from source.
	 */
	protected void adoptState(Model source) {
		// all we need to do is force a recalculation of intervals
		makeDirty();
	}

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

			double time = treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node);
			logL = diffusionModel.getLogLikelihood(parentTrait, childTrait, time);
//			System.err.println("Intermediate calculation:");
//			Taxon parentTaxon = treeModel.getNodeTaxon(parent);
//			System.err.println("Parent: "+ (parentTaxon == null ? "null" : parentTaxon.getId()));
//			Taxon childTaxon = treeModel.getNodeTaxon(node);
//			System.err.println("Child: "+ (childTaxon == null ? "null" : childTaxon.getId()));
//			System.err.println("logL = "+logL);
		} // else {

		// Add prior on root trait, should define in XML

//			double[][] rootPrecision = { {1,0} , {0,1} };
//			double[]   rootMean      = { 0,0 };
//			logL = MultivariateNormalDistribution.logPdf(childTrait,rootMean,rootPrecision,0,1);

//		}
		int childCount = treeModel.getChildCount(node);
		for (int i = 0; i < childCount; i++) {
			logL += traitLogLikelihood(treeModel.getChild(node, i));
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

			List<Integer> missingIndices = null;

			if (xo.hasAttribute(TRAIT_NAME)) {

				String traitName = xo.getStringAttribute(TRAIT_NAME);

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
			return new MultivariateTraitLikelihood(treeModel, diffusionModel,
					traitParameter, missingIndices);
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
//			new StringAttributeRule(TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
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

}

