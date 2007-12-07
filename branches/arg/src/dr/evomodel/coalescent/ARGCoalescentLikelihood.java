package dr.evomodel.coalescent;

import java.util.ArrayList;
import java.util.logging.Logger;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.ARGModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ARGCoalescentLikelihood extends CoalescentLikelihood{

	public static final int RECOMBINATION = 3;
	
	public static final String ARG_COALESCENT_MODEL = "argCoalescentLikelihood";
	public static final String RECOMBINATION_RATE = "recombinationRate";
	public static final String POPULATION_SIZE = "populationSize";
	public static final String CLOCK_RATE = "clockRate";
	public static final String ARG_MODEL = "argModel";
		
	private Parameter popSize;
	private Parameter recomRate;
	private Parameter clockRate;
	private ARGModel arg;
	
	private ArrayList<CoalescentInterval> intervals;
	private ArrayList<CoalescentInterval> storedIntervals;
	
	public ARGCoalescentLikelihood(Parameter popSize, Parameter recomRate, 
			Parameter clockRate, ARGModel arg, boolean setupIntervals) {
		//To make the super store state work correctly, we need to initalize the
		//the intervals first.
		super(ARG_COALESCENT_MODEL);
		
		this.popSize = popSize;
		this.recomRate = recomRate;
		this.clockRate = clockRate;
		this.arg = arg;
		addModel(arg);
		addParameter(popSize);
		addParameter(recomRate);
		addParameter(clockRate);
		
		intervals = new ArrayList<CoalescentInterval>(arg.getNodeCount());
		intervalsKnown = false;
		likelihoodKnown = false;
		
		if(setupIntervals){
			intervalsKnown = true;
			calculateIntervals();
		}
	}
	
	public void calculateIntervals(){
		intervals.clear();
		intervals.ensureCapacity(arg.getNodeCount());
		
		NodeRef x;
		for(int i = 0; i < arg.getInternalNodeCount(); i++){
			x = arg.getInternalNode(i);
			if(arg.isReassortment(x)){
				intervals.add(new CoalescentInterval(arg.getNodeHeight(x),RECOMBINATION));
			}else{
				intervals.add(new CoalescentInterval(arg.getNodeHeight(x),COALESCENT));
			}
		}
		for(int i = 0; i < arg.getExternalNodeCount(); i++){
			x = arg.getExternalNode(i);
			if(arg.getNodeHeight(x) > 0.0){
				intervals.add(new CoalescentInterval(arg.getNodeHeight(x),NEW_SAMPLE));
			}
		}
		dr.util.HeapSort.sort(intervals);
		
		double a;
		for(int i = 0; i < intervals.size() - 1; i++ ){
			 a = intervals.get(i).length;
			 intervals.get(i).length = a - intervals.get(i+1).length;
		}
		
		intervalsKnown = true;
		
	}
	
	 public void handleModelChangedEvent(Model model, Object object, int index) {
	        if (model == tree) {
	            intervalsKnown = false;
	        }
	        likelihoodKnown = false;
	 }
	 
	 public void handleParameterChangedEvent(Parameter parameter, int index){
		 likelihoodKnown = false;
	 }
	 
	 public void storeState(){
		 storedIntervals = new ArrayList<CoalescentInterval>(intervals.size());
		 for(CoalescentInterval interval : intervals){
			 storedIntervals.add(interval.clone());
		 }
		 makeDirty();
	}
	 
	 public void restoreState(){
		 intervals = storedIntervals;
		 storedIntervals.clear();
		 makeDirty();
	 }
	 
	 public double getLogLikelihood(){
		 if(likelihoodKnown){
			 return logLikelihood;
		 }
		 if(!intervalsKnown){
			 calculateIntervals();
		 }
		 likelihoodKnown = true;
		 logLikelihood = calculateLogLikelihood(popSize.getParameterValue(0), 
				 recomRate.getParameterValue(0),
				 clockRate.getParameterValue(0)); 
	 
		 return logLikelihood;
	 }
	 
	 private double calculateLogLikelihood(double pSize, double rRate, double clock){
		 double logLike = 0.0;
		 int numberOfTaxa = 2;
		 
		 for(CoalescentInterval interval: intervals){
			 double rate = (double)numberOfTaxa * (numberOfTaxa - 1);
			 rate = rate/pSize;
			 
			 logLike += Math.log(rate) - rate * interval.length * clock;
			 		 
					
			 if(interval.type == COALESCENT){
				 logLike += 0;
				 numberOfTaxa++;
			 }else if(interval.type == RECOMBINATION){
				 logLike += Math.log(rRate);
				 numberOfTaxa--;
			 }else{
				 throw new RuntimeException("Not implemented yet");
			 }
					 
		 }
		 		 
		 return logLike;
	 }
	 
	 private class CoalescentInterval implements Comparable<CoalescentInterval>, Cloneable{
		 public int type;
		 public double length;
		 
		 public CoalescentInterval(double length, int type){
			 this.length = length;
			 this.type = type;
		 }
		 
		 public int compareTo(CoalescentInterval a){
		 
		 if(a.length < this.length){
				return -1;
			}else if(a.length == this.length){
				Logger.getLogger("dr.evomodel.coalescent").severe("The current ARG Model " + 
						"has 2 internal nodes at the same height");
				return 0;
			}
			return 1;
		 }
		
		public String toString(){
			if(type == 0){
				return "(" + length + ", Coalescent)";
			}
			return "(" + length + ",Recombination)";
		}
		
		public CoalescentInterval clone(){
			return new CoalescentInterval(length, type);
		}
			 
	 }
	 
	 public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "A coalescent likelihood for an ARG model";
		}
		public Class getReturnType() {
			return ARGCoalescentLikelihood.class;
		}
		public String getParserName() {
			return ARG_COALESCENT_MODEL;
		}

		
		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
					
			
			XMLObject cxo = (XMLObject) xo.getChild(RECOMBINATION_RATE);
			Parameter rRate = (Parameter) cxo.getChild(Parameter.class);
			
			cxo = (XMLObject) xo.getChild(CLOCK_RATE);
			Parameter clock = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject) xo.getChild(POPULATION_SIZE);
			Parameter pSize = (Parameter)cxo.getChild(Parameter.class);
			
			cxo = (XMLObject) xo.getChild(ARG_MODEL);
			ARGModel argModel = (ARGModel)cxo.getChild(ARGModel.class);
			
			return new ARGCoalescentLikelihood(pSize,rRate,clock,argModel,false);
		}

		
		 
	 };
	 
	 
	
	

}
