package dr.evomodel.arg;

import dr.evomodel.arg.ARGModel.Node;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Parameter.ChangeType;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

public class ARGRatePrior extends AbstractModelLikelihood {

	public static final String ARG_RATE_PRIOR = "argRatePrior";
	public static final String SIGMA = "sigma";
	
	private final ARGModel arg;
	private final Parameter logNormalSigma;
	
	
	public ARGRatePrior(String name,ARGModel arg, Parameter sigma) {
		super(name);
		
		this.arg = arg;
		this.logNormalSigma = sigma;
		
		addModel(arg);
		addParameter(sigma);
		
	}
	
	public double[] generateValues(){
		
		double[] values = new double[arg.getNumberOfPartitions()];
		
		double sigma = logNormalSigma.getParameterValue(0);
				
		double oneOverSigma = 1.0/sigma;
		
		for(int i = 0; i < values.length; i++){
			values[i] = MathUtils.nextGamma(oneOverSigma,oneOverSigma);
		}
		
		
							
		return values;
	}

	public double getLogLikelihood() {
		return calculateLogLikelihood();
	}
	
	public double getAddHastingsRatio(double[] values){
		return -calculateLogLikelihood(values);
	}
	
	
	private double calculateLogLikelihood(double[] values){
		double logLike = 0;
				
		double sigma = logNormalSigma.getParameterValue(0);
		double oneOverSigma = 1.0/sigma;
		
		for(double d : values){
			logLike += GammaDistribution.logPdf(d, oneOverSigma, sigma);
		}
		
		return logLike;
	}
		
	private double calculateLogLikelihood(){
		double logLike = 0;
				
		
		for(int i = 0, n = arg.getNodeCount(); i < n; i++){
			Node x = (Node)arg.getNode(i);
						
			if(!x.isRoot() && x.isBifurcation()){
				
				double[] values = x.rateParameter.getParameterValues();
			
				logLike += calculateLogLikelihood(values);
			}
		}
		
		
		return logLike;
	}
	
	

	public Model getModel() {
		return this;
	}

	public void makeDirty() {
		
	}

	public LogColumn[] getColumns() {
		return new LogColumn[]{
				new NumberColumn(getId()){
					public double getDoubleValue() {
						return getLogLikelihood();
					}
				}
			};
	}

	public String getId() {
		return super.getId();
	}

	public void setId(String id) {
		super.setId(id);
	}

	protected void acceptState() {
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index,
			ChangeType type) {
	}
	
	protected void restoreState() {
		
	}
	
	protected void storeState() {
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			return ARGRatePrior.class;
		}
		
		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String id = xo.getAttribute("id","");
			
			
			Parameter sigma = (Parameter)xo.getChild(Parameter.class);
			
			
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			return new ARGRatePrior(id,arg,sigma);
		}

		public String getParserName() {
			return ARG_RATE_PRIOR;
		}
		
	};
	
	
}
