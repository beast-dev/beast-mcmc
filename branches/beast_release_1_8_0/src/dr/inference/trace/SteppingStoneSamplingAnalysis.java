package dr.inference.trace;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import dr.util.Attribute;
import dr.util.FileHelpers;

import dr.xml.*;

/**
 * @author Guy Baele
 */

public class SteppingStoneSamplingAnalysis {
	
    public static final String STEPPING_STONE_SAMPLING_ANALYSIS = "steppingStoneSamplingAnalysis";
    public static final String LIKELIHOOD_COLUMN = "likelihoodColumn";
    public static final String THETA_COLUMN = "thetaColumn";
    public static final String FORMAT = "%5.5g";
    
    private final String logLikelihoodName;
    private final List<Double> logLikelihoodSample;
    private final List<Double> thetaSample;
    private boolean logBayesFactorCalculated = false;
    private double logBayesFactor;
    private List<Double> maxLogLikelihood;
    private List<Double> mlContribution;
    private List<Double> orderedTheta;
    
    public SteppingStoneSamplingAnalysis(String logLikelihoodName, List<Double> logLikelihoodSample, List<Double> thetaSample) {
    	this.logLikelihoodSample = logLikelihoodSample;
        this.logLikelihoodName = logLikelihoodName;
        this.thetaSample = thetaSample;
    }
    
    public double getLogBayesFactor() {
        if (!logBayesFactorCalculated) {
            calculateBF();
        }
        return logBayesFactor;
    }
    
    private void calculateBF() {
    	
    	Map<Double, List<Double>> map = new HashMap<Double, List<Double>>();
        orderedTheta = new ArrayList<Double>();

        //only the log-likelihoods are needed to calculate the marginal likelihood
        for (int i = 0; i < logLikelihoodSample.size(); i++) {
        	if (!map.containsKey(thetaSample.get(i))) {
        		map.put(thetaSample.get(i), new ArrayList<Double>());
                orderedTheta.add(thetaSample.get(i));
            }
            map.get(thetaSample.get(i)).add(logLikelihoodSample.get(i));
        }

        //sort into ascending order
        Collections.sort(orderedTheta);

        //a list with the maxima of the log-likelihood values is constructed
        maxLogLikelihood = new ArrayList<Double>();
        for (double t : orderedTheta) {
        	List<Double> values = map.get(t);
            maxLogLikelihood.add(Collections.max(values));    
        }

        mlContribution = new ArrayList<Double>();
        
        logBayesFactor = 0.0;
        for (int i = 1; i < orderedTheta.size(); i++) {
        	double contribution = (orderedTheta.get(i) - orderedTheta.get(i-1)) * maxLogLikelihood.get(i-1);
        	logBayesFactor += contribution;
        	mlContribution.add(contribution);
        	//System.out.println(i + ": " + maxLogLikelihood.get(i-1));
        }
        //System.out.println(logBayesFactor);
        
        for (int i = 1; i < orderedTheta.size(); i++) {
        	double internalSum = 0.0;
        	for (int j = 0; j < map.get(orderedTheta.get(i-1)).size(); j++) {
        		internalSum += Math.exp((orderedTheta.get(i) - orderedTheta.get(i-1)) * (map.get(orderedTheta.get(i-1)).get(j) - maxLogLikelihood.get(i-1)));
        	}
        	internalSum /= map.get(orderedTheta.get(i-1)).size();
        	//System.out.print(orderedTheta.get(i) + "-" + orderedTheta.get(i-1) + ": " + Math.log(internalSum));
        	logBayesFactor += Math.log(internalSum);
        	mlContribution.set(i-1, mlContribution.get(i-1) + Math.log(internalSum));
        }
        
        logBayesFactorCalculated = true;
            
    }
    
    public String toString() {
        double bf = getLogBayesFactor();
        StringBuffer sb = new StringBuffer();
        sb.append("PathParameter\tMaxPathLikelihood\tMLContribution\n");
        for (int i = 0; i < orderedTheta.size(); ++i) {
            sb.append(String.format(FORMAT, orderedTheta.get(i)));
            sb.append("\t");
            sb.append(String.format(FORMAT, maxLogLikelihood.get(i)));
            sb.append("\t");
            if (i != (orderedTheta.size()-1)) {
            	sb.append(String.format(FORMAT, mlContribution.get(i)));
            }
            sb.append("\n");
        }

        sb.append("\nlog marginal likelihood (using stepping stone sampling) from " + logLikelihoodName + " = " + bf + "\n");
        return sb.toString();
    }
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	
    	public String getParserName() {
            return STEPPING_STONE_SAMPLING_ANALYSIS;
        }
    	
    	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    		
    		String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);
    		StringTokenizer tokenFileName = new StringTokenizer(fileName);
    		int numberOfFiles = tokenFileName.countTokens();
    		System.out.println(numberOfFiles + " file(s) found with marginal likelihood samples");
    		try {
    			
    			String likelihoodName = "";
    			List sampleLogLikelihood = null;
    			List sampleTheta = null;
    			
    			for (int j = 0; j < numberOfFiles; j++) {
    			
    				File file = new File(tokenFileName.nextToken());
    				String name = file.getName();
    				String parent = file.getParent();

    				if (!file.isAbsolute()) {
    					parent = System.getProperty("user.dir");
    				}

    				file = new File(parent, name);

    				fileName = file.getAbsolutePath();

    				XMLObject cxo = xo.getChild(LIKELIHOOD_COLUMN);
    				likelihoodName = cxo.getStringAttribute(Attribute.NAME);

    				cxo = xo.getChild(THETA_COLUMN);
    				String thetaName = cxo.getStringAttribute(Attribute.NAME);

    				LogFileTraces traces = new LogFileTraces(fileName, file);
    				traces.loadTraces();
                
    				int burnin = 0;
                
    				traces.setBurnIn(burnin);

    				int traceIndexLikelihood = -1;
    				int traceIndexTheta = -1;
    				for (int i = 0; i < traces.getTraceCount(); i++) {
    					String traceName = traces.getTraceName(i);
    					if (traceName.trim().equals(likelihoodName)) {
    						traceIndexLikelihood = i;
    					}
    					if (traceName.trim().equals(thetaName)) {
    						traceIndexTheta = i;
    					}
    				}

    				if (traceIndexLikelihood == -1) {
    					throw new XMLParseException("Column '" + likelihoodName + "' can not be found for " + getParserName() + " element.");
    				}

    				if (traceIndexTheta == -1) {
    					throw new XMLParseException("Column '" + thetaName + "' can not be found for " + getParserName() + " element.");
    				}

    				if (sampleLogLikelihood == null && sampleTheta == null) {
    					sampleLogLikelihood = traces.getValues(traceIndexLikelihood);
    					sampleTheta = traces.getValues(traceIndexTheta);
    				} else {
    					sampleLogLikelihood.addAll(traces.getValues(traceIndexLikelihood));
    					sampleTheta.addAll(traces.getValues(traceIndexTheta));
    				}
    			
    			}
    			
                SteppingStoneSamplingAnalysis analysis = new SteppingStoneSamplingAnalysis(likelihoodName, sampleLogLikelihood, sampleTheta);

                System.out.println(analysis.toString());

                return analysis;
    			
    		} catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (java.io.IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }
    		
    	}
    	
    	//************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Performs a trace analysis.";
        }

        public Class getReturnType() {
            return SteppingStoneSamplingAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately)"),
                new ElementRule(THETA_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
                new ElementRule(LIKELIHOOD_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
        };
    	
    };
	
}