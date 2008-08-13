package dr.evomodel.hgt;

import dr.evomodel.arg.ARGModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A likelihood function for the coalescent with recombination. Takes a tree and a demographic model.
 *
 * @author Marc Suchard
 * @version $Id: CoalescentWithRecombinationLikelihood.java,v 1.1.1.1.2.2 2006/11/06 01:38:30 msuchard Exp $
 */

public class HorizontalGeneTransferLikelihood extends AbstractModel implements Likelihood {

	private static final String HGT_LIKELIHOOD = "hgtLikelihood";

	private static final String GENE_GRAPH = "hgtGraph";

	private ARGModel arg;

	public HorizontalGeneTransferLikelihood(String name) {
		super(name);
	}

	public HorizontalGeneTransferLikelihood(ARGModel arg) {
		super(HGT_LIKELIHOOD);
		this.arg = arg;
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index) {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	protected void storeState() {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	protected void restoreState() {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	protected void acceptState() {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	public Model getModel() {
		return null;  //To change body of implemented methods use File | Settings | File Templates.
	}

	public double getLogLikelihood() {
		return calculateAnalyticalLogLikelihood();  //To change body of implemented methods use File | Settings | File Templates.
	}


	private final double calculateAnalyticalLogLikelihood() {
		int n = arg.getReassortmentNodeCount();
//		System.err.println("Prior for " + n + " reassortments.");
		double p = 0.001;
		double logL = n * Math.log(p);

		double height = arg.getNodeHeight(arg.getRoot());
		logL -= height * 1.0;   // Exp(1)
		// assumes a flat prior  over branch lengths given treeHeight
		//logL = Math.log(1.0/Math.pow(lambda,n-1));
//		System.err.println("COALESCENT PRIOR: end yoyoy");
		return logL;
	}


	public void makeDirty() {
		//To change body of implemented methods use File | Settings | File Templates.
	}

	// **************************************************************
	// Loggable IMPLEMENTATION
	// **************************************************************

	/**
	 * @return the log columns.
	 */
	public final dr.inference.loggers.LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[]{
				new LikelihoodColumn(getId())
		};
	}

	private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) {
			super(label);
		}

		public double getDoubleValue() {
			return getLogLikelihood();
		}
	}


	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return HGT_LIKELIHOOD;
		}

		public Object parseXMLObject(XMLObject xo) {

//            DemographicModel coalescentDemoModel = null;
			XMLObject cxo;
			cxo = (XMLObject) xo.getChild(GENE_GRAPH);
			ARGModel argModel = (ARGModel) cxo.getChild(ARGModel.class);

			return new HorizontalGeneTransferLikelihood(argModel);

//            return new CoalescentWithRecombinationLikelihood(argModel, coalescentDemoModel, recombinationDemoModel);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the likelihood of the ancestral recombination graph given the demographic function.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				//	new ElementRule(MODEL, new XMLSyntaxRule[] {
				//		new ElementRule(DemographicModel.class)
				//	}),
				new ElementRule(GENE_GRAPH, new XMLSyntaxRule[]{
						new ElementRule(ARGModel.class)
				}),
		};
	};

}


