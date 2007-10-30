package dr.inference.distribution;

import dr.evomodel.substmodel.PositiveDefiniteSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.DesignMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class MultivariateOUModel extends GeneralizedLinearModel {

	public static final String MVOU_MODEL = "multivariateOUModel";
	public static final String MVOU_TYPE = "MVOU";
	public static final String DATA = "data";
	public static final String TIME = "times";
	public static final String DESIGN = "design";

	private SubstitutionModel Q;
	private MatrixParameter gamma;
	private double[] time;
	private double[] deltaTime;
	private double[] design;
	private double[] W;
	private double[] initialPriorMean;
	private int K;
	private int numTimeSteps;

	private MultivariateNormalDistribution initialPrior;

	private boolean likelihoodKnown;
	private double logLikelihood;

	public MultivariateOUModel(SubstitutionModel substitutionModel, Parameter dependentParam,
	                           MatrixParameter gamma, double[] time, double[] design) {
		super(dependentParam);
		this.Q = substitutionModel;
		this.time = time;
		this.design = design;
		this.gamma = gamma;

		K = substitutionModel.getDataType().getStateCount();

		W = new double[K * K];

		initialPriorMean = new double[K]; // todo send this mean in constructor

		StringBuffer sb = new StringBuffer("Constructing a multivariate OU model:\n");
		sb.append("\tOutcome dimension = ");
		sb.append(K);
		Logger.getLogger("dr.inference.distribution").info(sb.toString());

		setupTimes();
		addParameter(gamma);
		addModel(substitutionModel);

	}

	private void setupTimes() {
		numTimeSteps = time.length / K - 1;
		deltaTime = new double[numTimeSteps];
		double currentTime = time[0];
		int index = 0;
		for (int i = 0; i < numTimeSteps; i++) {
			index += K;
			deltaTime[i] = time[index] - currentTime;
			currentTime = time[index];
		}
		Logger.getLogger("dr.inference.distribution").info(
				"\tTime increments: " + new Vector(deltaTime)
		);
	}

//	public final double getLogLikelihood() {
////		if (!likelihoodKnown) {
//			logLikelihood = calculateLogLikelihood();
////			likelihoodKnown = true;
////		}
//		return logLikelihood;
//	}

	public double calculateLogLikelihood(double[] x) {
		return calculateLogLikelihood();
	}

	public double calculateLogLikelihood() {

		double logLikelihood = 0;
		double[] previous = new double[K];
		double[] current = new double[K];
		double[] tmpHolder;
		double[][] G = gamma.getParameterAsMatrix();
		double[] theta = dependentParam.getParameterValues();
		double[] Xbeta = null;
		boolean hasEffects = getNumberOfEffects() > 0;

//		System.err.println("effects: "+hasEffects);

//		System.err.println("Xbeta = "+new Vector(Xbeta));

//		double currentTime = time[0];

//		double deltaTime = 0;

		// Prior on initial time-point
		int index = 0;

		if (!hasEffects) {
			for (int i = 0; i < K; i++)
				previous[i] = theta[index++];
		} else {
			Xbeta = getXBeta();
			for (int i = 0; i < K; i++) {
				previous[i] = theta[index] - Xbeta[index];
				index++;
			}
		}

		initialPrior = new MultivariateNormalDistribution(initialPriorMean, new Matrix(G).inverse().toComponents());
		logLikelihood += initialPrior.logPdf(previous);

//		System.err.println("initial point:");
//		System.err.println("\tX: "+new Vector(previous));
//		System.err.println("\tVar:\n"+new Matrix(G));
//		System.err.println("\tlogLike: "+logLikelihood);

		for (int timeStep = 0; timeStep < numTimeSteps; timeStep++) {

//			System.err.print("TimeStep #"+timeStep+" from "+currentTime+" to ");
//			currentTime += deltaTime[timeStep];
//			System.err.println(currentTime);


			Q.getTransitionProbabilities(-deltaTime[timeStep], W);
//			Q.getTransitionProbabilities(0, W);

//			System.err.println("\ttheta_0: "+new Vector(previous));

//			System.err.println("\tW:\n" + new Vector(W));
//			System.err.println("\tG:\n" + new Matrix(G));
//			System.exit(-1);

			double[] mean = new double[K];
			int u = 0;
			for (int i = 0; i < K; i++) {
				for (int j = 0; j < K; j++)
					mean[i] += W[u++] * previous[j];
			}

			double[][] WG = new double[K][K];
			for (int i = 0; i < K; i++) {
				for (int j = 0; j < K; j++) {
					for (int k = 0; k < K; k++)
						WG[i][j] += W[i * K + k] * G[k][j];       // W
				}
			}

//			System.err.println("\tWG:\n"+new Matrix(WG));

			double[][] WGWt = new double[K][K];
			for (int i = 0; i < K; i++) {
				for (int j = 0; j < K; j++) {
					for (int k = 0; k < K; k++)
						WGWt[i][j] += WG[i][k] * W[j * K + k];
				}
			}

//			System.err.println("\tWGWt:\n"+new Matrix(WGWt));

			for (int i = 0; i < K; i++) {
				for (int j = 0; j < K; j++)
					WGWt[i][j] = G[i][j] - WGWt[i][j];

			}

			// calculate density of current time step

			if (!hasEffects) {
				for (int i = 0; i < K; i++)
					current[i] = theta[index++];
			} else {
				for (int i = 0; i < K; i++) {
					current[i] = theta[index] - Xbeta[index];
					index++;
				}
			}

			MultivariateNormalDistribution density = new MultivariateNormalDistribution(
					mean, new Matrix(WGWt).inverse().toComponents());

//			System.err.println("\ttheta_1: "+new Vector(current));
//			System.err.println("\tmean: "+ new Vector(mean));
//			System.err.println("\tvariance:\n"+new Matrix(WGWt));

			double partialLogLikelihood = density.logPdf(current);
//			System.err.println("\tpartial-logLikelihood: "+partialLogLikelihood);

			if (partialLogLikelihood > 100) {
				System.err.println("got here partial");
				System.err.println("\ttheta_1: " + new Vector(current));
				System.err.println("\tmean: " + new Vector(mean));
				System.err.println("\tvariance:\n" + new Matrix(WGWt));
				System.err.println("\tW:\n" + new Vector(W));
				double[][] qMat = ((PositiveDefiniteSubstitutionModel) Q).getRates().getParameterAsMatrix();
//				System.err.println("");
				System.err.println("\tQ:\n" + new Matrix(qMat));
				try {
					System.err.println("\tQ pd? " + new Matrix(qMat).isPD());
				} catch (Exception e) {

				}
				System.err.println("\tG:\n" + new Matrix(G));

//			System.exit(-1);
			}

			logLikelihood += partialLogLikelihood;

			// move to next point
			tmpHolder = previous;
			previous = current;
			current = tmpHolder;

//			System.err.println("");

		}
//		System.err.println("MVOU logLike = "+logLikelihood);
//		System.exit(-1);

		if (logLikelihood > 100) {
			System.err.println("got here end");
			System.exit(-1);
		}

		return logLikelihood;
	}

	protected boolean confirmIndependentParameters() {
		return true;
	}

	protected boolean requiresScale() {
		return true;
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index) {
//		if( parameter == gamma )
//			 recalculate
//			initialPrior = new MultivariateNormalDistribution(initialPriorMean,
//					new Matrix(gamma.getParameterAsMatrix()).inverse().toComponents());
	}

	protected void storeState() {
	}

	protected void restoreState() {
	}

	protected void acceptState() {
	}

//	public double logPdf(Parameter x) {
//		return 0;
//	}

//	public double[][] getScaleMatrix() {
//		return new double[0][];
//	}
//
//	public double[] getMean() {
//		return new double[0];
//	}
//
//	public String getType() {
//		return MVOU_TYPE;
//	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MVOU_MODEL;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
			Parameter effectParameter = (Parameter) ((XMLObject) xo.getChild(DATA)).getChild(Parameter.class);
			Parameter timesParameter = (Parameter) ((XMLObject) xo.getChild(TIME)).getChild(Parameter.class);
			Parameter designParameter = (Parameter) ((XMLObject) xo.getChild(DESIGN)).getChild(Parameter.class);
			MatrixParameter gammaParameter = (MatrixParameter) xo.getChild(MatrixParameter.class);

			if (effectParameter.getDimension() != timesParameter.getDimension() ||
					effectParameter.getDimension() != designParameter.getDimension()) {
//				System.err.println("dim(effect) " +effectParameter.getDimension());
//				System.err.println("dim(times) "+timesParameter.getDimension());
//				System.err.println("dim(design) "+designParameter.getDimension());
				throw new XMLParseException("dim(" + effectParameter.getStatisticName() +
						") != dim(" + timesParameter.getStatisticName() + ") != dim(" + designParameter.getStatisticName() +
						") in " + xo.getName() + " element");
			}

			MultivariateOUModel glm = new MultivariateOUModel(substitutionModel, effectParameter, gammaParameter,
					timesParameter.getParameterValues(), designParameter.getParameterValues());

			addIndependentParameters(xo, glm, effectParameter);

			// todo Confirm that design vector is consistent with substitution model
			// todo Confirm that design vector is ordered 1,\ldots,K,1,\ldots,K, etc.

			return glm;

		}

		public void addIndependentParameters(XMLObject xo, GeneralizedLinearModel glm,
		                                     Parameter dependentParam) throws XMLParseException {
			int totalCount = xo.getChildCount();

			for (int i = 0; i < totalCount; i++) {
				if (xo.getChildName(i).compareTo(INDEPENDENT_VARIABLES) == 0) {
					XMLObject cxo = (XMLObject) xo.getChild(i);
					Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);
					DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
					checkDimensions(independentParam, dependentParam, designMatrix);
					glm.addIndependentParameter(independentParam, designMatrix);
				}
			}
		}

		private void checkDimensions(Parameter independentParam, Parameter dependentParam, DesignMatrix designMatrix)
				throws XMLParseException {
			if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
					(independentParam.getDimension() != designMatrix.getColumnDimension())) {
				System.err.println(dependentParam.getDimension());
				System.err.println(independentParam.getDimension());
				System.err.println(designMatrix.getRowDimension() + " rows");
				System.err.println(designMatrix.getColumnDimension() + " cols");
				throw new XMLParseException(
						"dim(" + dependentParam.getId() + ") != dim(" + designMatrix.getId() + " %*% " + independentParam.getId() + ")"
				);
			}
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(SubstitutionModel.class),
				new ElementRule(MatrixParameter.class),
				new ElementRule(DATA, new XMLSyntaxRule[]{
						new ElementRule(Parameter.class)}),
				new ElementRule(TIME, new XMLSyntaxRule[]{
						new ElementRule(Parameter.class)}),
				new ElementRule(DESIGN, new XMLSyntaxRule[]{
						new ElementRule(Parameter.class)}),
				new ElementRule(INDEPENDENT_VARIABLES,
						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, 0, 3),
		};

		public String getParserDescription() {
			return "Describes a multivariate OU process";
		}

		public Class getReturnType() {
			return MultivariateOUModel.class;
		}

	};

}
