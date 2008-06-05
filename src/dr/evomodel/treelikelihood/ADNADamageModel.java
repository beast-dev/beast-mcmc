package dr.evomodel.treelikelihood;

import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class ADNADamageModel extends TipPartialsModel {

	public static final String ADNA_DAMAGE_MODEL = "aDNADamageModel";
	public static final String BASE_DAMAGE_RATE = "baseDamageRate";
	public static final String AGE_RATE_FACTOR = "ageRateFactor";
	public static final String EXCLUDE = "exclude";
	public static final String INCLUDE = "include";

	public ADNADamageModel(TreeModel treeModel, TaxonList includeTaxa, TaxonList excludeTaxa, Parameter baseDamageRateParameter, Parameter ageRateFactorParameter) {
		super(ADNA_DAMAGE_MODEL, treeModel, includeTaxa, excludeTaxa);

		this.baseDamageRateParameter = baseDamageRateParameter;
		this.ageRateFactorParameter = ageRateFactorParameter;
	}

	public double[] getTipPartials(int nodeIndex) {
		if (updatePartials) {
			double base = baseDamageRateParameter.getParameterValue(0);
			double factor = ageRateFactorParameter.getParameterValue(0);

			for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
				int[] states = this.states[i];
				double[] partials = this.partials[i];

				if (!excluded[i]) {
					double age = treeModel.getNodeHeight(treeModel.getExternalNode(i));

					double pUndamaged = (1.0 - base) * Math.exp(-factor * age);

					int k = 0;
					for (int j = 0; j < patternCount; j++) {
						switch (states[j]) {
							case 0: // is an A
								partials[k] = pUndamaged;
								partials[k + 1] = 0.0;
								partials[k + 2] = 1.0 - pUndamaged;
								partials[k + 3] = 0.0;
								break;
							case 1: // is an C
								partials[k] = 0.0;
								partials[k + 1] = pUndamaged;
								partials[k + 2] = 0.0;
								partials[k + 3] = 1.0 - pUndamaged;
								break;
							case 2: // is an G
								partials[k] = 1.0 - pUndamaged;
								partials[k + 1] = 0.0;
								partials[k + 2] = pUndamaged;
								partials[k + 3] = 0.0;
								break;
							case 3: // is an T
								partials[k] = 0.0;
								partials[k + 1] = 1.0 - pUndamaged;
								partials[k + 2] = 0.0;
								partials[k + 3] = pUndamaged;
								break;
							default: // is an ambiguity
								partials[k] = 1.0;
								partials[k + 1] = 1.0;
								partials[k + 2] = 1.0;
								partials[k + 3] = 1.0;
						}
						k += stateCount;
					}
				} else {
					int k = 0;
					for (int j = 0; j < patternCount; j++) {
						switch (states[j]) {
							case 0: // is an A
								partials[k] = 1.0;
								partials[k + 1] = 0.0;
								partials[k + 2] = 0.0;
								partials[k + 3] = 0.0;
								break;
							case 1: // is an C
								partials[k] = 0.0;
								partials[k + 1] = 1.0;
								partials[k + 2] = 0.0;
								partials[k + 3] = 0.0;
								break;
							case 2: // is an G
								partials[k] = 0.0;
								partials[k + 1] = 0.0;
								partials[k + 2] = 1.0;
								partials[k + 3] = 0.0;
								break;
							case 3: // is an T
								partials[k] = 0.0;
								partials[k + 1] = 0.0;
								partials[k + 2] = 0.0;
								partials[k + 3] = 1.0;
								break;
							default: // is an ambiguity
								partials[k] = 1.0;
								partials[k + 1] = 1.0;
								partials[k + 2] = 1.0;
								partials[k + 3] = 1.0;
						}
						k += stateCount;
					}
				}
			}
			updatePartials = false;
		}

		return partials[nodeIndex];
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return ADNA_DAMAGE_MODEL; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);

			Parameter baseDamageRateParameter = (Parameter)xo.getElementFirstChild(BASE_DAMAGE_RATE);
			Parameter ageRateFactorParameter = (Parameter)xo.getElementFirstChild(AGE_RATE_FACTOR);

			TaxonList includeTaxa = null;
			TaxonList excludeTaxa = null;

			if (xo.hasChildNamed(INCLUDE)) {
				includeTaxa = (TaxonList)xo.getElementFirstChild(INCLUDE);
			}

			if (xo.hasChildNamed(EXCLUDE)) {
				excludeTaxa = (TaxonList)xo.getElementFirstChild(EXCLUDE);
			}

			ADNADamageModel aDNADamageModel =  new ADNADamageModel(treeModel, includeTaxa, excludeTaxa, baseDamageRateParameter, ageRateFactorParameter);

			System.out.println("Using aDNA damage model.");

			return aDNADamageModel;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return
					"This element returns a model that allows for post-mortem DNA damage.";
		}

		public Class getReturnType() { return ADNADamageModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new ElementRule(TreeModel.class),
				new ElementRule(BASE_DAMAGE_RATE, Parameter.class, "The base rate of accumulation of post-mortem damage", false),
				new ElementRule(AGE_RATE_FACTOR, Parameter.class, "The factor by which rate of damage scales with age of sample", false),
				new XORRule(
						new ElementRule(INCLUDE, TaxonList.class, "A set of taxa to which to apply the damage model to"),
						new ElementRule(EXCLUDE, TaxonList.class, "A set of taxa to which to not apply the damage model to")
						, true)
		};
	};

	private final Parameter baseDamageRateParameter;
	private final Parameter ageRateFactorParameter;
}