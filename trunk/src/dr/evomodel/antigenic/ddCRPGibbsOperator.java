package dr.evomodel.antigenic;

import java.util.logging.Logger;

import dr.evomodel.antigenic.NPAntigenicLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;


/**
 * A Gibbs operator for allocation of items to clusters under a distance dependent Chinese restaurant process.
 */

public class ddCRPGibbsOperator extends SimpleMCMCOperator implements GibbsOperator  {



	    public final static String DDCRP_GIBBS_OPERATOR = "ddCRPGibbsOperator";

	   
	    private final Parameter chiParameter;
	    private final double[][] depMatrix;
	    public NPAntigenicLikelihood modelLikelihood;
	    private Parameter links = null;
	    private Parameter assignments = null;
	    
	    
	    public ddCRPGibbsOperator(Parameter links, Parameter assignments,
	                                         Parameter chiParameter,
	                                         NPAntigenicLikelihood Likelihood,
	                                         double weight ) {
	     
	    	this.links = links;
	    	this.assignments = assignments;
	    	this.modelLikelihood = Likelihood;
	    	this.chiParameter = chiParameter;
	        this.depMatrix = Likelihood.getLogDepMatrix();
	        
	        for (int i=0;i<links.getDimension();i++){
	        	links.setParameterValue(i, i);
	        }   
	        
	        
	        setWeight(weight);

	       double[][] x=modelLikelihood.getData();
	       modelLikelihood.printInformtion(x[0][0]);
	        
	    }
	   
	      
	    
	    /**
	     * @return the parameter this operator acts on.
	     */
	    public Parameter getParameter() {
	        return (Parameter) links;
	    }

	    /**
	     * @return the Variable this operator acts on.
	   
	    public Variable getVariable() {
	        return clusteringParameter;
	    }
  */
	    /**
	     * change the parameter and return the hastings ratio.
	     */
	    public final double doOperation() {

	    	
	        int index = MathUtils.nextInt(links.getDimension());

	        int oldGroup = (int)assignments.getParameterValue(index);
	       
	      /*
	       * Set index customer link to index and all connected to it to a new assignment (min value empty)
	       */
	        int minEmp = minEmpty(modelLikelihood.getLogLikelihoodsVector());
	        links.setParameterValue(index, index);
	        int[] visited = connected(index, links);
	        
	        int ii = 0;
	        while (visited[ii]!=0){
	        	assignments.setParameterValue(visited[ii]-1, minEmp);
	        	ii++;
	        }
	        
	        
	        int maxFull = maxFull( modelLikelihood.getLogLikelihoodsVector());
	      
	       
	        /*
	         * Adjust likvector for group separated
	         */
	        
	        
	       modelLikelihood.setLogLikelihoodsVector(oldGroup,getLogLikGroup(oldGroup) );

	       modelLikelihood.setLogLikelihoodsVector(minEmp,getLogLikGroup(minEmp) );
	         
	       double[] liks = modelLikelihood.getLogLikelihoodsVector();
	       /*
	        * computing likelihoods of joint groups
	        */
	       
	       double[] crossedLiks = new double[maxFull];
	       
	       for (int ll=0;ll<maxFull;ll++ ){
	    	   if (ll!=minEmp){
	    		   crossedLiks[ll]=getLogLik2Group(ll,minEmp);
	    	   }
	       }
	       
	        
	        /*
	         * Add logPrior 
	         */
	        double[] logP = new double[links.getDimension()];

	        for (int jj=0; jj<links.getDimension(); jj++){
	        	logP[jj] += depMatrix[index][jj];
	        	
	        	int n = (int)assignments.getParameterValue(jj);
	        	if (n!= minEmp){
	        logP[jj]+=crossedLiks[n] - liks[n] - liks[minEmp];
	        	}
	        }
	      
	        logP[index]= Math.log(chiParameter.getParameterValue(0));
	        
	        
	        /*
	         * possibilidade de mandar p zero as probs muito pequenas
	         */
	      
	        
	        
	        /*
	         * Fazer o Gibbs sampling
	         */
	        
	        
	        this.rescale(logP); // Improve numerical stability
	        this.exp(logP); // Transform back to probability-scale

	        int k = MathUtils.randomChoicePDF(logP);
	        
	        links.setParameterValue(index, k);
	        
	        int newGroup = (int)assignments.getParameterValue(k);
	        ii = 0;
	        while (visited[ii]!=0){
	        	assignments.setParameterValue(visited[ii]-1, newGroup);
	        	ii++;
	        }
	        
	        
	        
	       /*
	        * updating conditional likelihood vector 
	        */
	        modelLikelihood.setLogLikelihoodsVector(newGroup, getLogLikGroup(newGroup));
	        if (newGroup!=minEmp){
	        	 modelLikelihood.setLogLikelihoodsVector(minEmp, 0);
	  	       
	        }
	        
	        return 0.0;
	
	    }


	    


	    
	    /*
	     * find min Empty
	     */
	    
	    public int minEmpty(double[] logLikVector){
	    int isEmpty=0;
	    int i =0;
	    	while (isEmpty==0){
	    if(logLikVector[i]==0){
	    	isEmpty=1;}
	    else { 
	    	if(i==logLikVector.length-1){isEmpty=1;}
		   i++;}
	    	}
	    return i;
	    }
	   
	
	    
	    /*
	     * find max Full
	     */
	    
	    
	    public int maxFull(double[] logLikVector){
		    int isEmpty=1;
		    int i =logLikVector.length-1;
		    	while (isEmpty==1){
		    if(logLikVector[i]!=0){
		    	isEmpty=0;}
		    else {i--;}
		    	}
		    return i+1;
		    }
	 /*
	  * find customers conected to i   
	  */
	    
	    public int[] connected(int i, Parameter clusteringParameter){
	    	int n =  clusteringParameter.getDimension();
	    	int[] visited = new int[n+1]; 
	  	    visited[0]=i+1;
	  	    int tv=1;
	  	    
	  	    for(int j=0;j<n;j++){
	  	    	if(visited[j]!=0){
	  	    		int curr = visited[j]-1;
		  	    	
	  	    		/*look forward
		  	    	*/
	  	    		
	  	    	int forward = (int) clusteringParameter.getParameterValue(curr);
	  	    	visited[tv] = forward+1;
	  	    	tv++;
	  	    	
	  	    	for(int ii=0; ii<tv-1; ii++){
	  	    	if(visited[ii]==forward+1){
	  	    		tv=tv-1;
	  	    		visited[tv]=0;
	  	    	}
	  	    	}
	  	    	
	  	    	/*look back
	  	    	*/
	  	    	for (int jj=0; jj<n;jj++){
	  	    		if((int)clusteringParameter.getParameterValue(jj)==curr){ 
	  	    			visited[tv]= jj+1;
	  	    			tv++;
	  	    			
	  	    			for(int ii=0; ii<tv-1; ii++){
	  	  	  	    	if(visited[ii]==jj+1){
	  	  	  	    		tv--;
	  	  	  	    		visited[tv]=0;
	  	  	  	    	}		
	  	    		}
	  
	  	    	}
	  	    	}
	  	    	
	  	    	}}
	  	return visited;
	  	
	    }
	    
	
	    
			  
	    private void printInformtion(Parameter par) {
	         StringBuffer sb = new StringBuffer("parameter \n");
	         	 for(int j=0; j<par.getDimension(); j++){
	        		 sb.append(par.getParameterValue(j));
	         }
	         
	         Logger.getLogger("dr.evomodel").info(sb.toString()); };
	        
		  
	    
	    
	    
		 public double getLogLikGroup(int groupNumber){   
		    double L =0.0;
			 
		    
				  int ngroup=0;
				  for (int i=0;i<assignments.getDimension(); i++){
				  if((int) assignments.getParameterValue(i) == groupNumber){
				  ngroup++;}}
				  
				  
				  if (ngroup != 0){
				  double[] group = new double[2*ngroup];
				  
				  int count = 0;
				  for (int i=0;i<assignments.getDimension(); i++){
					  if((int) assignments.getParameterValue(i) == groupNumber){
						  group[count] = modelLikelihood.getData()[i][0];
						  group[ngroup+count] = modelLikelihood.getData()[i][1];
						  count+=1;}}
					  
				  
				  double[][] var = new double[2*ngroup][2*ngroup];
				  double[] mean = new double[2*ngroup];
				  
				  double m0 = modelLikelihood.getPriorMean().getParameterValue(0);
				  double m1 = modelLikelihood.getPriorMean().getParameterValue(1);
				  double vp = modelLikelihood.getPriorVar().getParameterValue(0);
				  double vc = modelLikelihood.getClusterVar().getParameterValue(0);
					 
				  
				  for (int i=0; i<ngroup; i++){
					  mean[i]=m0;
					  mean[ngroup+i]=m1;
					  
				  
					  for (int l=0;l<ngroup;l++){
						  var[i][ngroup+l]=0;
						  var[ngroup+i][l]=0;
						  
					  if (l==i){var[i][l]= vp+ vc;
					  			var[ngroup+i][ngroup+l]= vp+ vc;}
					  
					  else { var[i][l] = vp;
					  		var[ngroup+i][ngroup+l]= vp;}
					  	}
				  }
				  
				  
				  
				  
				  double[][] precision = new SymmetricMatrix(var).inverse().toComponents(); 
				   L = new MultivariateNormalDistribution(mean, precision).logPdf(group);
				  
				  }			  
				  
				 
		    
		    return L;
		    
		 }   
		    
		    
		    
		 
		   
		 public double getLogLik2Group(int group1, int group2){   
		    double L =0.0;
			 
		    
				  int ngroup1=0;
				  for (int i=0;i<assignments.getDimension(); i++){
				  if((int) assignments.getParameterValue(i) == group1 ){
				  ngroup1++;}}
				  
				  int ngroup2=0;
				  for (int i=0;i<assignments.getDimension(); i++){
				  if((int) assignments.getParameterValue(i) == group2 ){
				  ngroup2++;}}
				  
				  int ngroup = (ngroup1+ngroup2);
				  
				  if (ngroup != 0){
				  double[] group = new double[2*ngroup];
				  
				  int count = 0;
				  for (int i=0;i<assignments.getDimension(); i++){
					  if((int) assignments.getParameterValue(i) == group1 ){
						  group[count] = modelLikelihood.getData()[i][0];
						  group[count+ngroup] = modelLikelihood.getData()[i][1];
						  count+=1;}}
				
				  for (int i=0;i<assignments.getDimension(); i++){
					  if((int) assignments.getParameterValue(i) == group2 ){
						  group[count] = modelLikelihood.getData()[i][0];
						  group[count+ngroup] = modelLikelihood.getData()[i][1];
						  count+=1;}}
				
				  
				  
				  
				  double[][] var = new double[2*ngroup][2*ngroup];
				  double[] mean = new double[2*ngroup];
				  
				  double m0 = modelLikelihood.getPriorMean().getParameterValue(0);
				  double m1 = modelLikelihood.getPriorMean().getParameterValue(1);
				  double vp = modelLikelihood.getPriorVar().getParameterValue(0);
				  double vc = modelLikelihood.getClusterVar().getParameterValue(0);
					 
				  
				  for (int i=0; i<ngroup; i++){
					  mean[i]=m0;
					  mean[i+ngroup]=m1;
					  
				  for (int l=0;l<ngroup;l++){
					  var[i][ngroup+l]=0;
					  var[ngroup+i][l]=0;
					  
				  if (l==i){var[i][l]= vp+ vc;
				  			var[ngroup+i][ngroup+l]= vp+ vc;}
				  
				  else { var[i][l] = vp;
				  		var[ngroup+i][ngroup+l]= vp;}
				  	}
			  }
				  
				  
				  double[][] precision = new SymmetricMatrix(var).inverse().toComponents(); 
				   L = new MultivariateNormalDistribution(mean, precision).logPdf(group);
				  
				  }	
			 
	    
	    return L;
	    
	 }  
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	    
	    
	    private void exp(double[] logX) {
	        for (int i = 0; i < logX.length; ++i) {
	            logX[i] = Math.exp(logX[i]);
	          //  if(logX[i]<1E-5){logX[i]=0;}
	        }
	    }

	    private void rescale(double[] logX) {
	        double max = this.max(logX);
	        for (int i = 0; i < logX.length; ++i) {
	            logX[i] -= max;
	        }
	    }

	    private double max(double[] x) {
	        double max = x[0];
	        for (double xi : x) {
	            if (xi > max) {
	                max = xi;
	            }
	        }
	        return max;
	    }

	    //MCMCOperator INTERFACE
	    public final String getOperatorName() {
	        return DDCRP_GIBBS_OPERATOR;
	    }

	    public final void optimize(double targetProb) {

	        throw new RuntimeException("This operator cannot be optimized!");
	    }

	    public boolean isOptimizing() {
	        return false;
	    }

	    public void setOptimizing(boolean opt) {
	        throw new RuntimeException("This operator cannot be optimized!");
	    }

	    public double getMinimumAcceptanceLevel() {
	        return 0.1;
	    }

	    public double getMaximumAcceptanceLevel() {
	        return 0.4;
	    }

	    public double getMinimumGoodAcceptanceLevel() {
	        return 0.20;
	    }

	    public double getMaximumGoodAcceptanceLevel() {
	        return 0.30;
	    }

	    public String getPerformanceSuggestion() {
	        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
	            return "";
	        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
	            return "";
	        } else {
	            return "";
	        }
	    }

	
	   
	

	    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
	        public final static String CHI = "chi";
	        public final static String LIKELIHOOD = "likelihood";
	        public final static String ASSIGNMENTS = "assignments";
	        public final static String LINKS = "links";
	        public final static String DEP_MATRIX = "depMatrix";

	        public String getParserName() {
	            return DDCRP_GIBBS_OPERATOR;
	        }

	        /* (non-Javadoc)
	         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
	         */
	        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

	            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

	            
	            XMLObject cxo = xo.getChild(ASSIGNMENTS);
		        Parameter assignments = (Parameter) cxo.getChild(Parameter.class);

		        cxo = xo.getChild(LINKS);
		        Parameter links = (Parameter) cxo.getChild(Parameter.class);

		        cxo = xo.getChild(CHI);
		        Parameter chiParameter = (Parameter) cxo.getChild(Parameter.class);

		        
		        
		        cxo = xo.getChild(LIKELIHOOD);
		        NPAntigenicLikelihood likelihood = (NPAntigenicLikelihood)cxo.getChild(NPAntigenicLikelihood.class);
		        
		  
		        
		        
	            return new ddCRPGibbsOperator( links, assignments,
                        chiParameter,likelihood, weight);

	        }

	        //************************************************************************
	        // AbstractXMLObjectParser implementation
	        //************************************************************************

	        public String getParserDescription() {
	            return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
	        }

	        public Class getReturnType() {
	            return ddCRPGibbsOperator.class;
	        }


	        public XMLSyntaxRule[] getSyntaxRules() {
	            return rules;
	        }

	        private final XMLSyntaxRule[] rules = {
	                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
	                new ElementRule(CHI, new XMLSyntaxRule[] {
	                        new ElementRule(Parameter.class),
	                }),
	                new ElementRule(LIKELIHOOD, new XMLSyntaxRule[] {
	                        new ElementRule(Likelihood.class),
	                }, true),
	                new ElementRule(ASSIGNMENTS,
    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    	            new ElementRule(LINKS,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    	 
	    };
	    
	    };


	
	    public int getStepCount() {
	        return 1;
	    }

	}




