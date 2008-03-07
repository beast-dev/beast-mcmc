package dr.evomodel.coalescent;


import dr.evomodel.tree.ARGModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ARGUniformPrior extends ARGCoalescentLikelihood{

	public static final String ARG_UNIFORM_PRIOR = "argUniformPrior";
		
	private double[] argNumber;
	
	ARGUniformPrior(ARGModel arg) {
		super(ARG_UNIFORM_PRIOR,arg);
				
		addModel(arg);
				
		argNumber = new double[5];
		for(int i = 0, n = arg.getExternalNodeCount(); i < argNumber.length; i++){
			argNumber[i] = Math.log(numberARGS(n,i));
		}
		
	}
	
	public double getLogLikelihood(){
		if(likelihoodKnown){
			return logLikelihood;
		}
		
		likelihoodKnown = true;
		logLikelihood = calculateLogLikelihood();
				
		if(!currentARGValid()){
			logLikelihood = Double.NEGATIVE_INFINITY;
		}
		
		if(arg.getReassortmentNodeCount() > 3){
			logLikelihood = Double.NEGATIVE_INFINITY;
		}
		
		return logLikelihood;
	}
		
	public double calculateLogLikelihood(){
		
		double treeHeight = arg.getNodeHeight(arg.getRoot());
		int internalNodes = arg.getInternalNodeCount() - 1;
		
		
		double logLike = logFactorial(internalNodes) - (double)internalNodes*Math.log(treeHeight)
		 	- argNumber[arg.getReassortmentNodeCount()];
		
		assert !Double.isInfinite(logLike) && !Double.isNaN(logLike);

		return logLike;
	}
		
	private double logFactorial(int n){
		double rValue = 0;
		
		for(int i = n; i > 0; i--){
			rValue += Math.log(i);
		}
		return rValue; 
	}
	
	
	private int numberARGS(int taxa, int argNumber){
		int x = taxa;
		int n = 2*argNumber + taxa - 1;
		
		return shurikoRecursion(x,n);
	}
	
	private int shurikoRecursion(int x, int n){
		int a = 0;
		if(x == 0){
			a = 0;
		}else if(x == 1){
			if(n == 0){
				a = 1;
			}else{
				a = 0;
			}
		}else if(n == 0){
			if(x == 1){
				a = 1;
			}else{
				a = 0;
			}
		}else if(x == n + 1){
			a = x*(x-1)/2*shurikoRecursion(x-1,n-1);
		}else{
			a = x*shurikoRecursion(x+1,n-1) + x*(x-1)/2*shurikoRecursion(x-1,n-1);
		}
		return a;
	}
	 public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

			public String getParserDescription() {
				return "A uniform prior for an ARG model";
			}
			public Class getReturnType() {
				return ARGUniformPrior.class;
			}
			public String getParserName() {
				return ARG_UNIFORM_PRIOR;
			}
			
			public XMLSyntaxRule[] getSyntaxRules(){
				return rules;
			}
					
			private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
					new ElementRule(ARGModel.class),
							
			};
			
			public Object parseXMLObject(XMLObject xo) throws XMLParseException {
				ARGModel argModel = (ARGModel)xo.getChild(ARGModel.class);
				
				return new ARGUniformPrior(argModel);
			}

			
			 
		 };
	
	
}
