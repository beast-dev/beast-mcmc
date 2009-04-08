package dr.evomodel.coalescent;


import dr.evolution.io.NewickImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeTrace;
import dr.inference.loggers.Logger;
import dr.inference.model.Parameter;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;

public class GMRFFixedGridImportanceSampler implements Logger{

    public static final String FIXED_GRID_IMPORTANCE_SAMPLER = "gmrfFixedGridImportanceSampler";
    public static final String PROPOSAL_FILE_NAME = "proposals";
    public static final String LIKELIHOOD_FILE_NAME = "likelihood";
    public static final String TREE_FILE_NAME = "trees";
    public static final String OUTPUT_FILE_NAME = "output";
    public static final String COALESCENT_INTERVALS = "coalescentIntervals";
    public static final String BURN_IN = "burnIn";
    
    private LogFileTraces[] lft;
    private TreeTrace treeTrace;
    private PrintWriter pw;
    private int burnin;
    private int coalescentIntervals;
    private double gridHeight;
   
    public GMRFFixedGridImportanceSampler(String[] logFileNames, Reader treeReader,
    		PrintWriter pw, int bI, int coalescentIntervals, double gridHeight) 
    		throws TraceException, IOException, ImportException{
        lft = new LogFileTraces[2];
        
        for(int i = 0; i < 2; i++){
        	File file = new File(logFileNames[i]);
           	lft[i] = new LogFileTraces(logFileNames[i], file);
           	lft[i].loadTraces();
           	lft[i].setBurnIn(bI);
        }
        
        treeTrace = TreeTrace.loadTreeTrace(treeReader);
        
    	this.coalescentIntervals = coalescentIntervals;
        this.burnin = bI;
        this.pw = pw;
        this.gridHeight = gridHeight;
    }
    
    
    public double getPrior(double[] pSize, double distance, double precision){
    	
    	double a = 0;
    	for(int i = 1; i < pSize.length; i++){
    		a += Math.pow(pSize[i-1] - pSize[i],2);
    	}
    	    	
    	double b = 0.5 * precision * a / distance;
    	
    	double prior = (pSize.length - 1) * 0.5 * Math.log(precision) - b;
    	
//    	prior -= (pSize.length - 1) / 2.0 * GMRFSkyrideLikelihood.LOG_TWO_TIMES_PI;
    	    	
    	return prior;
    }
    
    public double getIntegratedPrior(double[] pSize, double distance, double alpha, double beta){
    	double a = 0;
    	for(int i = 1; i < pSize.length; i++){
    		a += Math.pow(pSize[i-1] - pSize[i],2)/(distance*2);
    	}
    	
//    	double c = -(pSize.length + 1.0)/2.0 * Math.log(1 + a);
//    	
//    	if(c > -13){
//    		for(int i = 1; i < pSize.length; i++){
//        		System.out.println(pSize[i]);
//        	}
//    		System.exit(-1);
//    	}
    	
    	return -(pSize.length/2.0 + alpha - 1.0) * Math.log(beta + a);
    	
    	
    }
    
    public double getNonInformativeIntegratedPrior(double[] pSize, double distance){
    	double a = 0;
    	for(int i = 1; i < pSize.length; i++){
    		a += Math.pow(pSize[i-1] - pSize[i],2)/(distance*2);
    	}
    	
//    	double c = -(pSize.length + 1.0)/2.0 * Math.log(1 + a);
//    	
//    	if(c > -13){
//    		for(int i = 1; i < pSize.length; i++){
//        		System.out.println(pSize[i]);
//        	}
//    		System.exit(-1);
//    	}
    	
    	return -(pSize.length/2.0 - 1.0) * Math.log(a);
    	
    	
    }
    
    

    public void report() throws IOException, ImportException{
    	ArrayList<Double> mylist = new ArrayList<Double>();
    	
    	for(int i = 0; i < lft[0].getStateCount(); i++){
    		mylist.add(getImportanceWeight(i));
    	}
    	    	
    	double max = Collections.max(mylist);
    	
    	for(int i = 0 ; i < mylist.size(); i++){
    		mylist.set(i, Math.exp(-max + mylist.get(i)));
    	}
    	
    	double sum = 0;
    	for(double a : mylist)
    		sum += a;
    	    	
    	for(int i = 0; i < mylist.size(); i++){
    		mylist.set(i, mylist.get(i)/sum);
    		
    	}
    	for(int i = 1; i < mylist.size(); i++){
    		mylist.set(i, mylist.get(i-1) + mylist.get(i));
    	}
    	
    	for(int i = 0; i < 500; i++){
    		double draw = MathUtils.nextDouble();
    		
    		int index = 0;
    		
    		for(index = 0; index < mylist.size(); index++){
    			if(draw < mylist.get(index)){
    				break;
    			}
    		}
    		
//    		System.out.println(index);
    		
    		for(int k = 0; k < 10; k++){
        		System.out.print(lft[0].getStateValue(k,index) + " ");
    		}
    		System.out.println("");
    	}
    	
    }
    
    private double getImportanceWeight(int index){
    	Tree tree = treeTrace.getTree(index, burnin);
    	
    	
    	double[] fixedPopSizes = new double[lft[0].getTraceCount() - 1];
    	double[] fixedHeights = new double[fixedPopSizes.length];
    	for(int i = 0; i < fixedPopSizes.length; i++){
    		fixedPopSizes[i] = lft[0].getStateValue(i,index);
    		fixedHeights[i] = (i+1)*gridHeight/(double)fixedPopSizes.length;
    		
    	}
    	 
    	GMRFFixedGridLikelihood like = new GMRFFixedGridLikelihood(
    				tree,new Parameter.Default(fixedPopSizes),
    				new Parameter.Default(fixedHeights),coalescentIntervals);
    	
    	    	    	
    	double top = getImportanceWeightTop(like, index, fixedPopSizes);
    	double bottom = getImportanceWeightBottom(like, index, fixedPopSizes);
    	    	
    	
    	
    	return (top - bottom);
    }
    
    public double getImportanceWeightTop(GMRFFixedGridLikelihood like, int index, double[] fixedPopSizes){
    	double rescaledPrecision = lft[1].getStateValue(0, index);
	    	
    	double top = like.getLogLikelihood()
//    	    + getNonInformativeIntegratedPrior(fixedPopSizes, gridHeight/fixedPopSizes.length)
    		+ getIntegratedPrior(fixedPopSizes, gridHeight/fixedPopSizes.length,0.001,0.001)
//    		+ getPrior(fixedPopSizes,gridHeight/fixedPopSizes.length,rescaledPrecision)
		   ;
    	 	
    	
    	return top;
    }
    
    public double getImportanceWeightBottom(GMRFFixedGridLikelihood like, int index, double[] fixedPopSizes){
    	double a = lft[1].getStateValue(2, index) + lft[0].getStateValue(fixedPopSizes.length, index);
    	
    	return a;
//    	return 0;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return "Generates importance samples from a fixed grid";
        }

        public Class getReturnType() {
            return GMRFFixedGridImportanceSampler.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
        		AttributeRule.newStringRule(PROPOSAL_FILE_NAME),
        		AttributeRule.newStringRule(LIKELIHOOD_FILE_NAME),
        		AttributeRule.newStringRule(TREE_FILE_NAME),
        		AttributeRule.newStringRule(OUTPUT_FILE_NAME),
        		AttributeRule.newStringRule(COALESCENT_INTERVALS),
        		AttributeRule.newDoubleRule(
        				dr.evomodel.coalescent.AbstractFixedGridLogger.GRID_STOP_TIME),
        		AttributeRule.newIntegerRule(BURN_IN, true),
    		};
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            PrintWriter pw;
            Reader reader;
            String[] inputFileName = {
            		xo.getStringAttribute(PROPOSAL_FILE_NAME),
            		xo.getStringAttribute(LIKELIHOOD_FILE_NAME),
            		xo.getStringAttribute(TREE_FILE_NAME),
            };		
            String outputFileName = xo.getStringAttribute(OUTPUT_FILE_NAME);
            
            int coalescentIntervals = xo.getIntegerAttribute(COALESCENT_INTERVALS);
            double gridHeight = xo.getDoubleAttribute(
            		dr.evomodel.coalescent.AbstractFixedGridLogger.GRID_STOP_TIME);
            try {
            	File[] file = new File[3];
                String[] name = new String[3];
                String[] parent = new String[3];
                
                for(int i = 0; i < 3; i++){
                	file[i] = new File(inputFileName[i]);
                	name[i] = file[i].getName();
                	parent[i] = file[i].getParent();
                	
                	if (!file[i].isAbsolute()) {
                        parent[i] = System.getProperty("user.dir");
                    }
                	
                	file[i] = new File(parent[i], name[i]);
                    inputFileName[i] = file[i].getName();
                }
                                
                File outputFile = new File(outputFileName);
				String outputName = outputFile.getName();
				String outputParent = outputFile.getParent();

				if (!outputFile.isAbsolute()) {
					outputParent = System.getProperty("user.dir");
				}
				reader = new FileReader(new File(parent[2], name[2]));
				pw = new PrintWriter(new FileOutputStream(new File(outputParent, outputName)));

            } catch (java.io.IOException e) {
                throw new XMLParseException("Something is wrong");
            }

            int burnin = xo.getAttribute(BURN_IN, -1);
           
            GMRFFixedGridImportanceSampler analysis = null;
            
            try{
            	analysis = new GMRFFixedGridImportanceSampler(inputFileName, reader, 
            				pw, burnin, coalescentIntervals, gridHeight);
            	analysis.report();
            }catch(Exception e){
            	System.err.println(e.getMessage());
            }
            System.out.println("");
            System.out.flush();

            return analysis;
        }

        public String getParserName() {
            return FIXED_GRID_IMPORTANCE_SAMPLER;
        }

    };

	public void log(int state) {
		// TODO Auto-generated method stub
		
	}

	public void startLogging() {
		// TODO Auto-generated method stub
		
	}

	public void stopLogging() {
		// TODO Auto-generated method stub
		
	}
}
