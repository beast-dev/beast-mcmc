package dr.evomodel.coalescent;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.math.util.MathUtils;
import dr.evomodel.tree.ARGModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ARGUniformPrior extends ARGCoalescentLikelihood{

	public static final String ARG_UNIFORM_PRIOR = "argUniformPrior";
	public static final String INITIAL_CALCULATIONS = "initialCalculations";

	private ArrayList<Double> argNumber;
	
	ARGUniformPrior(ARGModel arg, int max, int initial) {
		super(ARG_UNIFORM_PRIOR,arg, max);
				
		addModel(arg);
				
		argNumber = new ArrayList<Double>(initial);
		
		Logger.getLogger("dr.evomodel").info("Calulating ARGs for uniform prior:");
		Logger.getLogger("dr.evomodel").info("0 reassorments");
		
		double temp = 0;
		for(int i = arg.getExternalNodeCount() ; i > 1; i--)
			temp += Math.log(i*(i-1)/2.0);
		argNumber.add(temp);
		
		
		for(int i = 1, n = arg.getExternalNodeCount(); i < initial; i++){
			argNumber.add(logNumberARGS(n,i));
			Logger.getLogger("dr.evomodel").info(i + " reassorments ");
		}
		
	}
	
	public double getLogARGNumber(int i){
		if(i >= argNumber.size()){
			argNumber.add(logNumberARGS(arg.getExternalNodeCount(), i));
		}
		return argNumber.get(i);
	}
	
	public double getLogLikelihood(){
		if(likelihoodKnown){
			return logLikelihood;
		}
		
		likelihoodKnown = true;
		logLikelihood = calculateLogLikelihood();
				
		
		
		if(arg.getReassortmentNodeCount() > maxReassortments)
			logLikelihood = Double.NEGATIVE_INFINITY;
		else
			logLikelihood = calculateLogLikelihood();
		
		if(!currentARGValid(true)){
			logLikelihood = Double.NEGATIVE_INFINITY;
		}
		
		return logLikelihood;
	}
		
	public double calculateLogLikelihood(){
		
		double treeHeight = arg.getNodeHeight(arg.getRoot());
		int internalNodes = arg.getInternalNodeCount() - 1;
		
		
		double logLike = logFactorial(internalNodes) - (double)internalNodes*Math.log(treeHeight)
		 	- getLogARGNumber(arg.getReassortmentNodeCount());
		
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
	
	public static double logNumberARGS(int start, int reassortments){
		int[] max = new int[start - 3 + reassortments*2];
		int[] x = new int[max.length]; 
	
		int i = 0;
		while(i < reassortments){
			x[i] = max[i] = 1;
			i++;
		}
		while(i < max.length){
			x[i] = max[i] = -1;
			i++;
		}
		
		
		double approx = 0;
		while(x[0] != -9 && !stopCombination(x,start)){
			if(testCombination(x, start)){
				int[] y = generateValues(x,start);
				approx += reduceThenDivide(y, generateValues(max,start));
			}
			nextCombination(x);
		}
		
		approx = Math.log(approx);
		
		int[] y = new int[max.length + 2];
		for(i = 0; i < max.length; i++)
			y[i] = max[i];
		y[y.length - 2] = y[y.length - 1] = -1;
		
		max = generateValues(y,start);
		
		for(int k = 0; k < y.length; k++)
			approx += Math.log(max[k]);
		
		return approx;
		
	}
		
	private static double reduceThenDivide(int[] top, int[] bottom){
			
		if(false){
			for(int i = 0; i < top.length; i++){
				for(int j = 0; j < bottom.length; j++){
					int gcd = MathUtils.gcd(top[i],bottom[j]);
					
					if(gcd > 1){
						top[i] = top[i]/gcd;
						bottom[j] = bottom[j]/gcd;
					}
				}
				
			}
			
		}
		
		Arrays.sort(top);
		Arrays.sort(bottom);
			
		
		double a = 1;
		for(int i = 0; i < top.length; i++)
			a *= (double)top[i] / bottom[i];
		return a;
		
	}
	
	private static int[] generateValues(int[] x, int start){
		int[] y = new int[x.length];
		
		for(int i = 0; i < x.length; i++){
			if(x[i] == 1)
				y[i] = start;
			else
				y[i] = start*(start-1)/2;
			
			start += x[i];
			
		}
		return y;
	}
	
	private static boolean testCombination(int[] x, int start){
		
		for(int i = 0; i < x.length; i++){
			start += x[i];
			if(start == 1)
				return false;
		}
		return true;
		
	}
	
	private static boolean stopCombination(int[] x, int start){
		for(int i = 0; i < x.length; i++){
			if(x[i] == -1){
				start--;
				if(start == 1){
					return true;
				}
			}else{
				break;
			}
		}
		return false;
	}
	
	private static void nextCombination(int[] x){
		if(x[x.length - 1] == -1){
			int i = x.length - 1;
			
			while(i > -1){
				if(x[i] == 1){
					x[i] = -1;
					x[i+1] = 1;
					return;
				}else
					i--;
				
			}
		}else{
			int endOnes = 0;
			int i = x.length - 1;
			while(x[i] == 1){
				endOnes++;
				i--;
			}
			int nextOne = -1;
			while(i > -1){
				if(x[i] == 1){
					nextOne = i;
					break;
				}else
					i--;
				
			}
			if(nextOne == -1){
				x[0] = -9;
				return;
			}
			
			x[nextOne] = -1;
			x[nextOne + 1] = 1;
			
			for(i = 0; i < endOnes; i++)
				x[i + nextOne + 2] = 1;
			
			i = nextOne + 2 + endOnes;
			
			while(i < x.length){
				x[i] = -1;
				i++;
			}
		}
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
				
				int max = Integer.MAX_VALUE;
				if(xo.hasAttribute(MAX_REASSORTMENTS)){
					max = xo.getIntegerAttribute(MAX_REASSORTMENTS);
				}
				
				int initial = 5;
				if(xo.hasAttribute(INITIAL_CALCULATIONS)){
					initial = xo.getIntegerAttribute(INITIAL_CALCULATIONS);
				}
				
				return new ARGUniformPrior(argModel,max,initial+1);
			}

			
			 
		 };
	
	
}
