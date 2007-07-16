package dr.evomodel.continuous;

import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.MultivariateNormalDistribution;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 3, 2007
 * Time: 4:14:13 PM
 * To change this template use File | Settings | File Templates.
 */


public class MultivariateDiffusionModel extends AbstractModel {

	public static final String DIFFUSION_PROCESS = "multivariateDiffusionModel";
	public static final String DIFFUSION_CONSTANT = "precisionMatrix";
	public static final String BIAS = "mu";

//	private MatrixParameter diffusionPrecisionMatrixParameter;

	/**
	 * Construct a diffusion model.
	 */

	public MultivariateDiffusionModel(MatrixParameter diffusionPrecisionMatrixParameter) {

		super(DIFFUSION_PROCESS);

		this.diffusionPrecisionMatrixParameter = diffusionPrecisionMatrixParameter;
//		dim = diffusionPrecisionMatrixParameter.getRowDimension();
		calculatePrecisionInfo();
		addParameter(diffusionPrecisionMatrixParameter);

	}

	public MultivariateDiffusionModel() {
		super(DIFFUSION_PROCESS);
	}


	public MatrixParameter getPrecisionMatrixParameter() {
		return diffusionPrecisionMatrixParameter;
	}

	public double[][] getPrecisionmatrix() {
		return diffusionPrecisionMatrixParameter.getParameterAsMatrix();
	}

	private void printVector(double[] v) {
		for (double d : v) {
			System.err.print(d + " ");
		}
		System.err.println("");
	}


	/**
	 * @return the log likelihood of going from start to stop in the given time
	 */
	public double getLogLikelihood(double[] start, double[] stop, double time) {

		double logDet = Math.log(determinatePrecisionMatrix / time);
		return MultivariateNormalDistribution.logPdf(stop, start,
				diffusionPrecisionMatrix, logDet, time);
	}


	private void calculatePrecisionInfo() {
		diffusionPrecisionMatrix = diffusionPrecisionMatrixParameter.getParameterAsMatrix();
		determinatePrecisionMatrix =
				MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(
						diffusionPrecisionMatrix);
	}

	/**
	 * @return the bias of this diffusion process.
	 */
/*	private double getBias() {
		if (biasParameter == null) return 0.0;
		return biasParameter.getParameterValue(0);
	}*/

	// *****************************************************************
	// Interface Model
	// *****************************************************************
	public void handleModelChangedEvent(Model model, Object object, int index) {
		// no intermediates need to be recalculated...
	}

	public void handleParameterChangedEvent(Parameter parameter, int index) {
		calculatePrecisionInfo();
	}

	protected void storeState() {
		savedDeterminatePrecisionMatrix = determinatePrecisionMatrix;
		savedDiffusionPrecisionMatrix = diffusionPrecisionMatrix;
	}

	protected void restoreState() {
		determinatePrecisionMatrix = savedDeterminatePrecisionMatrix;
		diffusionPrecisionMatrix = savedDiffusionPrecisionMatrix;
	}

	protected void acceptState() {
	} // no additional state needs accepting

	protected void adoptState(Model source) {
	} // no additional state needs adopting

	// **************************************************************
	// XMLElement IMPLEMENTATION
	// **************************************************************

	public Element createElement(Document document) {
		throw new RuntimeException("Not implemented!");
	}

	// **************************************************************
	// XMLObjectParser
	// **************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return DIFFUSION_PROCESS;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			XMLObject cxo = (XMLObject) xo.getChild(DIFFUSION_CONSTANT);
			MatrixParameter diffusionParam = (MatrixParameter) cxo.getChild(MatrixParameter.class);
//			MatrixParameter diffusionParam = null;
			return new MultivariateDiffusionModel(diffusionParam);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "Describes a multivariate normal diffusion process.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(DIFFUSION_CONSTANT,
						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
		};

		public Class getReturnType() {
			return MultivariateDiffusionModel.class;
		}
	};

	// **************************************************************
	// Private instance variables
	// **************************************************************

	private MatrixParameter diffusionPrecisionMatrixParameter;
	private double determinatePrecisionMatrix;
	private double savedDeterminatePrecisionMatrix;
	private double[][] diffusionPrecisionMatrix;
	private double[][] savedDiffusionPrecisionMatrix;
//	private Parameter biasParameter;
//	private int dim;
//	private double logNormalization;
}

