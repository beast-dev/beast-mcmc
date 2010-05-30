package dr.inference.trace;

import dr.xml.*;
import dr.util.FileHelpers;
import dr.util.Attribute;
import dr.inferencexml.trace.MarginalLikelihoodAnalysisParser;
import dr.stats.DiscreteStatistics;

import java.io.File;
import java.io.FileNotFoundException;

import jebl.evolution.alignments.Alignment;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsPerSiteAnalysis {
    public static final String DNDS_PERSITE_ANALYSIS = "dNdSPerSiteAnalysis";
    public static final String SPERSITE_COLUMNS = "conditionalSperSite";
    public static final String UNCONDITIONAL_S_COLUMN = "unconditionalSperSite";
    public static final String NPERSITE_COLUMNS = "conditionalNperSite";
    public static final String UNCONDITIONAL_N_COLUMN = "unconditionalNperSite";
    //public static final String ALIGNMENT = "alignment";

    DnDsPerSiteAnalysis(double[][]sampleSperSite, double[]unconditionalS,
                        double[][]sampleNperSite, double[]unconditionalN){
        this.sampleSperSite = sampleSperSite;
        this.unconditionalS = unconditionalS;
        this.sampleNperSite = sampleNperSite;
        this.unconditionalN = unconditionalN;
    }

    public double[][] getdSPerSiteSample() {
        if (!dNAnddSPerSiteSampleCollected) {
            collectDnAndDsPerSite();
        }
        return dSPerSiteSample;
    }

    public double[][] getdNPerSiteSample() {
        if (!dNAnddSPerSiteSampleCollected) {
            collectDnAndDsPerSite();
        }
        return dNPerSiteSample;
    }

    public double[][] getdNdSRatioPerSiteSample() {
        if (!dNAnddSPerSiteSampleCollected) {
            collectDnAndDsPerSite();
        }
        return dNdSRatioPerSiteSample;
    }

    private void collectDnAndDsPerSite() {

        for (int r = 0; r < sampleSperSite.length; r++) {
            for (int c = 0; c < sampleSperSite[0].length; c++) {
                dSPerSiteSample[r][c] = sampleSperSite[r][c]/unconditionalS[c];
                dNPerSiteSample[r][c] = sampleNperSite[r][c]/unconditionalN[c];
                dNdSRatioPerSiteSample[r][c] = dSPerSiteSample[r][c]/dNPerSiteSample[r][c];
            }
        }
    }

    public String toString() {
        double[] dSPerSiteMean = mean(getdSPerSiteSample());
        double[] dNPerSiteMean = mean(getdNPerSiteSample());
        double[] dNdSRatioPerSiteMean = mean(getdNdSRatioPerSiteSample());

        StringBuffer sb = new StringBuffer();
        sb.append("site\tdS\tdN\tdN/dS\n");
        for(int i = 0; i < dNdSRatioPerSiteMean[i]; i++) {
            sb.append(i+1+"\t"+dSPerSiteMean+"\t"+dNPerSiteMean+"\t"+dNdSRatioPerSiteMean+"\n");    
        }

        return sb.toString();

    }

    private static double[] mean(double[][] x)    {
         double[] returnArray = null;
         for (int i = 0; i < x.length; i++) {
            returnArray[i] = DiscreteStatistics.mean(x[i]);
         }
         return returnArray;
     }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DNDS_PERSITE_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);
            try {

                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

                file = new File(parent, name);

                fileName = file.getAbsolutePath();

                XMLObject cxo = xo.getChild(SPERSITE_COLUMNS);
                String sPerSiteColumnName = cxo.getStringAttribute(Attribute.NAME);

                cxo = xo.getChild(UNCONDITIONAL_S_COLUMN);
                String unconditionalSName = cxo.getStringAttribute(Attribute.NAME);

                cxo = xo.getChild(NPERSITE_COLUMNS);
                String nPerSiteColumnName = cxo.getStringAttribute(Attribute.NAME);

                cxo = xo.getChild(UNCONDITIONAL_N_COLUMN);
                String unconditionalNName = cxo.getStringAttribute(Attribute.NAME);

                //cxo = xo.getChild(ALIGNMENT);
                //String alignmentName = cxo.getStringAttribute(Attribute.NAME);

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute(MarginalLikelihoodAnalysisParser.BURN_IN, maxState / 10);

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 5;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 20%");
                }

                traces.setBurnIn(burnin);

                int traceIndexSperSite = -1;
                int traceIndexUnconditionalS = -1;
                int traceIndexNperSite = -1;
                int traceIndexUnconditionalN = -1;

                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().contains(sPerSiteColumnName)) {
                        traceIndexSperSite = i;
                    }
                    if (traceName.trim().equals(unconditionalSName)) {
                        traceIndexUnconditionalS = i;
                    }
                    if (traceName.trim().contains(nPerSiteColumnName)) {
                        traceIndexNperSite = i;
                    }
                    if (traceName.trim().equals(unconditionalNName)) {
                        traceIndexUnconditionalN = i;
                    }
                }

                if (traceIndexSperSite == -1) {
                    throw new XMLParseException(sPerSiteColumnName+" columns can not be found for " + getParserName() + " element.");
                }
                if (traceIndexUnconditionalS == -1) {
                    throw new XMLParseException("Column '" + unconditionalSName + "' can not be found for " + getParserName() + " element.");
                }

                if (traceIndexNperSite == -1) {
                    throw new XMLParseException(nPerSiteColumnName+" columns can not be found for " + getParserName() + " element.");
                }
                if (traceIndexUnconditionalN == -1) {
                    throw new XMLParseException("Column '" + unconditionalNName + "' can not be found for " + getParserName() + " element.");
                }


                int traceEndIndexSperSite = -1;
                int numberOfSperSite = 1;
                for (int i = traceIndexSperSite; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().contains(sPerSiteColumnName)) {
                        traceEndIndexSperSite = i;
                    }
                }
                numberOfSperSite = 1 + (traceEndIndexSperSite - traceIndexSperSite);

                int traceEndIndexNperSite = -1;
                int numberOfNperSite = 1;
                for (int i = traceIndexNperSite; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().contains(nPerSiteColumnName)) {
                        traceEndIndexNperSite =  i;
                    }
                }
                numberOfNperSite = 1 + (traceEndIndexNperSite - traceIndexNperSite);

                if (numberOfSperSite != numberOfNperSite) {
                    throw new XMLParseException("different number of sites for N (" +numberOfNperSite+") and S ("+numberOfSperSite+") counts for " + getParserName() + " ??");
                }


                double sampleSperSite[][] = new double[numberOfSperSite][traces.getStateCount()];
                double sampleNperSite[][] = new double[numberOfNperSite][traces.getStateCount()];
                double unconditionalS[] = new double[traces.getStateCount()];
                double unconditionalN[] = new double[traces.getStateCount()];

                //collect all arrays
                for(int a = traceIndexSperSite; a < numberOfSperSite; a++){
                    traces.getValues(a, sampleSperSite[a]);
                }
                for(int b = traceIndexNperSite; b < numberOfNperSite; b++){
                    traces.getValues(b, sampleNperSite[b]);
                }
                traces.getValues(traceIndexUnconditionalS, unconditionalS);
                traces.getValues(traceIndexUnconditionalN, unconditionalN);

                DnDsPerSiteAnalysis analysis = new DnDsPerSiteAnalysis(
                        sampleSperSite, unconditionalS,
                        sampleNperSite, unconditionalN);

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
            return "Performs a trace dN dS analysis.";
        }

        public Class getReturnType() {
            return DnDsPerSiteAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
                new ElementRule(UNCONDITIONAL_S_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
                new ElementRule(UNCONDITIONAL_N_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
        };
    };
    private final double[][] sampleSperSite;
    private final double[][] sampleNperSite;
    private final double[] unconditionalS;
    private final double[] unconditionalN;

    private boolean dNAnddSPerSiteSampleCollected = false;
    private double[][] dSPerSiteSample;
    private double[][] dNPerSiteSample;
    private double[][] dNdSRatioPerSiteSample;

}
