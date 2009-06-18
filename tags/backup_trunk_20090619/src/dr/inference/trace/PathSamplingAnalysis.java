package dr.inference.trace;

import dr.evomodelxml.LoggerParser;
import dr.util.Attribute;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

/**
 * @author Marc A. Suchard
 * @author Alexander Alekseyenko
 */
public class PathSamplingAnalysis {

    public static final String PATH_SAMPLING_ANALYSIS = "pathSamplingAnalysis";
    public static final String LIKELIHOOD_COLUMN = "likelihoodColumn";
    public static final String THETA_COLUMN = "thetaColumn";

    PathSamplingAnalysis(double[] logLikelihoodSample, String logLikelihoodName, double[] thetaSample, String thetaName) {
        this.logLikelihoodSample = logLikelihoodSample;
        this.logLikelihoodName = logLikelihoodName;
        this.thetaSample = thetaSample;
        this.thetaName = thetaName;
    }

    public double getLogBayesFactor() {
        if (!logBayesFactorCalculated) {
            calculate();
        }
        return logBayesFactor;
    }

    private void calculate() {

//  R code from Alex Alekseyenko
//
//  psMLE = function(likelihood, pathParameter){
//      y=tapply(likelihood,pathParameter,mean)
//      L = length(y)
//      midpoints = (y[1:(L-1)] + y[2:L])/2
//      x = as.double(names(y))
//      widths = (x[2:L] - x[1:(L-1)])
//      sum(widths*midpoints)
//  }

//        List<Double> meanLogLikelihood = new ArrayList<Double>();
//        List<Double> theta = new ArrayList<Double>();

        Map<Double, List<Double>> map = new HashMap<Double,List<Double>>();
        List<Double> orderedTheta = new ArrayList<Double>();

        for(int i=0; i<logLikelihoodSample.length; i++) {
            if( !map.containsKey(thetaSample[i])) {
                map.put(thetaSample[i], new ArrayList<Double>());
                orderedTheta.add(thetaSample[i]);
            }
            map.get(thetaSample[i]).add(logLikelihoodSample[i]);
        }
        List<Double> meanLogLikelihood = new ArrayList<Double>();
         for(double t : orderedTheta) {
            double totalMean = 0;
            int lengthMean = 0;
            List<Double> values = map.get(t);
            for(double v: values) {
                totalMean += v;
                lengthMean++;
            }
            meanLogLikelihood.add(totalMean/lengthMean);
        }
//        System.err.println("classes = "+map.keySet().size());
//        System.exit(-1);

        logBayesFactor = 0;
        for(int i=0; i<meanLogLikelihood.size()-1; i++)
                logBayesFactor += (meanLogLikelihood.get(i+1) + meanLogLikelihood.get(i)) / 2.0 *
                                  (orderedTheta.get(i+1) - orderedTheta.get(i));

        logBayesFactorCalculated = true;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("log Bayes factor from ");
        sb.append(logLikelihoodName);
        sb.append(" = ");
        sb.append(String.valueOf(getLogBayesFactor()));
        return sb.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PATH_SAMPLING_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(LoggerParser.FILE_NAME);
            try {

                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

                file = new File(parent, name);

                fileName = file.getAbsolutePath();

                XMLObject cxo = (XMLObject) xo.getChild(LIKELIHOOD_COLUMN);
                String likelihoodName = cxo.getStringAttribute(Attribute.NAME);

                cxo = (XMLObject) xo.getChild(THETA_COLUMN);
                String thetaName = cxo.getStringAttribute(Attribute.NAME);

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute(MarginalLikelihoodAnalysis.BURN_IN, maxState / 10);

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
                }

//                burnin = 0;

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

                double sampleLogLikelihood[] = new double[traces.getStateCount()];
                double sampleTheta[] = new double[traces.getStateCount()];
                traces.getValues(traceIndexLikelihood, sampleLogLikelihood);
                traces.getValues(traceIndexTheta, sampleTheta);

                PathSamplingAnalysis analysis = new PathSamplingAnalysis(
                        sampleLogLikelihood, likelihoodName,
                        sampleTheta, thetaName);

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
            return PathSamplingAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(LoggerParser.FILE_NAME, "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
        };
    };

    private boolean logBayesFactorCalculated = false;
    private double logBayesFactor;
    private double[] logLikelihoodSample;
    private double[] thetaSample;
    private String logLikelihoodName;
    private String thetaName;
}
