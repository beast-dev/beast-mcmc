package dr.evomodel.mixturemodels;

import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;
import java.util.ArrayList;
import java.math.BigInteger;

public class HDPDrawsFromCommonBaseOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String HDP_DRAWS_OPERATOR = "hdpDrawsFromCommonBaseOperator";
    private HDPPolyaUrn hdp;
    private CompoundParameter tableCounts;
    //private CompoundParameter commonBaseDistWeights;
    private CompoundParameter stickProportions;
    private Parameter groupAssignments;
    private Parameter categoriesParameter;
    public final ArrayList<ArrayList<BigInteger>> cachedStirlingNumbers;
    private double pathWeight = 1.0;
    public static final double ACCURACY_THRESHOLD = 1E-12;

    public HDPDrawsFromCommonBaseOperator(HDPPolyaUrn hdp, double weight) {

        this.hdp = hdp;
        this.tableCounts = hdp.getTableCounts();
        //this.commonBaseDistWeights = hdp.getCommonBaseDistWeights();
        this.stickProportions = hdp.getStickProportions();
        this.categoriesParameter = hdp.getCategoriesParameter();
        this.groupAssignments = hdp.getGroupAssignments();

        cachedStirlingNumbers = new ArrayList<ArrayList<BigInteger>>();
        cachedStirlingNumbers.add(new ArrayList<BigInteger>());
        cachedStirlingNumbers.get(0).add(0, BigInteger.valueOf(1));
        cachedStirlingNumbers.get(0).add(1,BigInteger.valueOf(0));
        cachedStirlingNumbers.add(new ArrayList<>());
        cachedStirlingNumbers.get(1).add(0, BigInteger.valueOf(0));
        cachedStirlingNumbers.get(1).add(1,BigInteger.valueOf(1));

        setWeight(weight);
    }


    public double doOperation() {

        //System.err.println("log Stirling number: " + logOfBigInteger(getStirlingNumber(629,1)));

        // Let G and K denote the current number of groups and categories
        // tablesCounts is a compoundParameter consisting of G parameters
        // dimension k of parameter g indicates the the number of tables in group g corresp. to cat. k

        // counts[g][k] is the number of sites corresponding to group g and cat k
        int[][] counts = new int[hdp.maxGroupCount][hdp.maxCategoryCount];


        double commonMass = hdp.getCommonMass().getParameterValue(0);
        //double commonBaseDistMass = hdp.getCommonBaseDistMass().getParameterValue(0);

        for(int i = 0; i < categoriesParameter.getSize(); i++){
            counts[(int)groupAssignments.getParameterValue(i)][(int)categoriesParameter.getParameterValue(i)]++;
        }

        int[] isCatActive = new int[hdp.maxCategoryCount];

        double[] comBaseDistWeights = new double[hdp.maxCategoryCount+1];

        // initialize w_u = 1
        comBaseDistWeights[comBaseDistWeights.length-1] = 1;

        double weightsSum = 0;

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            for(int g = 0; g < hdp.maxGroupCount; g++){
                if(counts[g][k] > 0){
                    isCatActive[k] = 1;
                    break;
                }
            }
            // at this point, we know if cat K is active
            if(isCatActive[k] == 1) {
                // w_k = \beta_k*w_u
                comBaseDistWeights[k] = stickProportions.getParameter(k).getParameterValue(0) * comBaseDistWeights[comBaseDistWeights.length - 1];
                // w_u = (1-\beta_k)*w_u
                comBaseDistWeights[comBaseDistWeights.length - 1] = (1 - stickProportions.getParameter(k).getParameterValue(0)) * comBaseDistWeights[comBaseDistWeights.length - 1];
            }
            weightsSum = weightsSum + comBaseDistWeights[k];
        }

        weightsSum = weightsSum + comBaseDistWeights[comBaseDistWeights.length-1];

        if(Math.abs(weightsSum-1) > ACCURACY_THRESHOLD){
            throw new RuntimeException("common base dist weights must sum to 1");
        }

        // then iterate over all nonzero tableCounts, updating one count at a time according to its complete conditional

        for(int g = 0; g < hdp.maxGroupCount; g++){

            for(int c = 0; c < hdp.maxCategoryCount; c++) {

                //System.err.println("group " + g + " and category " + c + " with counts[g][c] == " +  counts[g][c]);
                //System.err.println("commonMass: " + commonMass + " commonBaseDistMass: " + commonBaseDistMass);

                // if counts[g][c] > 0, then m_{gc} can assume values 1,...,counts[g][c]
                if (counts[g][c] > 0) {

                    //double comBaseDistWeight = commonBaseDistWeights.getParameter(c).getParameterValue(0);

                    if(comBaseDistWeights[c] <= 0 || comBaseDistWeights[c] > 1){
                        throw new RuntimeException("commonBaseDistWeight has inappropriate value");
                    }

                    // tableCountProbs[m] is conditional probability that m_{gk} == m + 1
                    // m = 0,...,counts[g][k]-1
                    double[] tableCountProbs = new double[counts[g][c]];

                    // log of tableCountProbs
                    double[] logTableCountProbs = new double[counts[g][c]];

                    for (int m = 0; m < tableCountProbs.length; m++) {
                        logTableCountProbs[m] = (m + 1) * Math.log(commonMass*comBaseDistWeights[c]);

                        //System.err.println("tableCountProbs[" + m + "] after first term: " + Math.exp(logTableCountProbs[m]));

                        for (int i = 0; i < counts[g][c]; i++) {
                            logTableCountProbs[m] = logTableCountProbs[m] - pathWeight*Math.log(commonMass*comBaseDistWeights[c] + i);
                        }

                        //System.err.println("tableCountProbs[" + m + "] after second term: " + Math.exp(logTableCountProbs[m]));

                        //int tableCountsForCatExcludingGroup = tableCountsForEachCat[c] - (int) tableCounts.getParameter(g).getParameterValue(c);

                        //if(tableCountsForCatExcludingGroup == 0){
                        //    logTableCountProbs[m] = logTableCountProbs[m] + Math.log(commonBaseDistMass);
                        //}

                        //for (int i = 1; i < (m + 1); i++) {
                        //    //tableCountProbs[m] = tableCountProbs[m] * (tableCountsForCatExcludingGroup + i);
                        //    if (tableCountsForCatExcludingGroup == 0) {
                        //        if (i != 0) {
                        //            logTableCountProbs[m] = logTableCountProbs[m] + Math.log(tableCountsForCatExcludingGroup + i);
                        //        }
                        //    } else {
                        //        logTableCountProbs[m] = logTableCountProbs[m] + Math.log(tableCountsForCatExcludingGroup + i);
                        //    }
                        //}

                        //System.err.println("tableCountProbs[" + m + "] after third term: " + Math.exp(logTableCountProbs[m]));
                        //System.err.println("tableCountsForCatExcludingGroup: " + tableCountsForCatExcludingGroup);
                        //System.err.println("totalTableCount: " + totalTableCount);
                        //System.err.println("totalTableCountForGroup: " + tableCountsForEachGroup[g]);

                        //int totalTableCountExcludingGroup = totalTableCount - tableCountsForEachGroup[g];

                        //System.err.println("totalTableCountExcludingGroup: " + totalTableCountExcludingGroup);

                        //for (int i = 1; i < (m + 1); i++) {
                        //for (int i = 0; i < (m + 1); i++) {
                        //tableCountProbs[m] = tableCountProbs[m] * (1 / (totalTableCountExcludingGroup + commonBaseDistMass + i));
                        // logTableCountProbs[m] = logTableCountProbs[m] - Math.log(totalTableCountExcludingGroup + commonBaseDistMass + i);
                        //}

                        //System.err.println("tableCountProbs[" + m + "] after fourth term: " + Math.exp(logTableCountProbs[m]));

                        double logsn = logOfBigInteger(getStirlingNumber(counts[g][c], m + 1));

                        logTableCountProbs[m] = logTableCountProbs[m] + pathWeight*logsn;

                        if(Double.isNaN(logTableCountProbs[m]) || logTableCountProbs[m] == Double.POSITIVE_INFINITY){
                            System.out.println("logTableCountProbs[" + m +"]: " + logTableCountProbs[m]);
                            System.out.println("logsn: " + logsn);
                            System.out.println("counts[g][c]: " + counts[g][c]);
                            System.out.println("m + 1: " + (m+1));
                            System.out.println("getStirlingNumber(counts[g][c], m + 1): " + getStirlingNumber(counts[g][c], m + 1));
                        }

                        //System.err.println("tableCountProbs[" + m + "] final: " + Math.exp(logTableCountProbs[m]));
                    }

                    //System.err.println("tableCountProbs.length: " + tableCountProbs.length);
                    double sumTableCountProbs = 0;
                    for (int i = 0; i < tableCountProbs.length; i++) {
                        tableCountProbs[i] = Math.exp(logTableCountProbs[i]);
                        //sumTableCountProbs = sumTableCountProbs + tableCountProbs[i];
                        //System.err.println("tablesCountProbs[" + i + "]: " + tableCountProbs[i]);
                    }
                    //System.err.println("sumTableCountProbs: " + sumTableCountProbs);


                    int newTableCount = MathUtils.randomChoicePDF(tableCountProbs) + 1;

                    //int oldTableCount = (int) tableCounts.getParameter(g).getParameterValue(c);

                    if (newTableCount == 0) {
                        throw new RuntimeException("New table count is 0");
                    }

                    //set new value
                    tableCounts.getParameter(g).setParameterValue(c, newTableCount);

                    // after updating parameter, update counts for use in next iteration
                    //totalTableCount = totalTableCount - oldTableCount + newTableCount;
                    //tableCountsForEachCat[c] = tableCountsForEachCat[c] - oldTableCount + newTableCount;
                    //if (tableCountsForEachCat[c] == 0) {
                    //    System.err.println("tableCountsForEachCat[k] is 0");
                    //    System.exit(0);
                    //}
                    //tableCountsForEachGroup[g] = tableCountsForEachGroup[g] - oldTableCount + newTableCount;

                }
            }
        }

        return 0.0;
    }

    // compute Stirling number of the first kind

    public BigInteger getStirlingNumber(int n, int m) {
        if(cachedStirlingNumbers.size() < (n+1)){
            for(int i = cachedStirlingNumbers.size(); i <= n; i++){
                cachedStirlingNumbers.add(new ArrayList<BigInteger>());
                cachedStirlingNumbers.get(i).add(0, BigInteger.valueOf(0));

                // now make sure we have computed s(i,1), s(i,2), ..., s(i,m)
                for(int j = 1; j <= m; j++){
                    if(j > i){
                        cachedStirlingNumbers.get(i).add(j,BigInteger.valueOf(0));
                    }else{
                        cachedStirlingNumbers.get(i).add(j,
                                cachedStirlingNumbers.get(i-1).get(j-1).add(
                                        BigInteger.valueOf(i-1).multiply(
                                                cachedStirlingNumbers.get(i-1).get(j)))

                        );
                    }
                }
            }
        }

        if(cachedStirlingNumbers.get(n).size() < (m+1)){
            // now make sure we have computed s(n,1), s(n,2), ..., s(n,m)
            for(int j = cachedStirlingNumbers.get(n).size(); j <= m; j++){
                if(j > n){
                    cachedStirlingNumbers.get(n).add(j,BigInteger.valueOf(0));
                }else{
                    //System.err.println("j: " + j + " m: " + m + " n: " + n);
                    //System.err.println(cachedStirlingNumbers.get(n-1).size());
                    //cachedStirlingNumbers.get(n).add(j,
                    //            cachedStirlingNumbers.get(n-1).get(j-1)
                    //                    + (n-1)*cachedStirlingNumbers.get(n-1).get(j)
                    //    );
                    cachedStirlingNumbers.get(n).add(j,
                            getStirlingNumber(n-1,j-1).add(
                                    BigInteger.valueOf(n-1).multiply(
                                            getStirlingNumber(n-1,j)))
                    );
                }
            }
        }

        //System.err.println("getting stirling number with n = " + n + " and m = " + m + " : " + cachedStirlingNumbers.get(n).get(m));
        return cachedStirlingNumbers.get(n).get(m);
    }


    public static double logOfBigInteger(BigInteger val) {
        int b = val.bitLength() - 1022; // any value between 60 and 1023 works
        if (b > 0)
            val = val.shiftRight(b);
        double res = Math.log(val.doubleValue());
        return b > 0 ? res + b * Math.log(2.0) : res;
    }

    public void setPathParameter(double beta) {
        if (beta < 0 || beta > 1) {
            throw new IllegalArgumentException("Illegal path weight of " + beta);
        }
        pathWeight = beta;
    }

    public String getOperatorName() {
        return HDP_DRAWS_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return HDP_DRAWS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            HDPPolyaUrn hdp = (HDPPolyaUrn)xo.getChild(HDPPolyaUrn.class);

            return new HDPDrawsFromCommonBaseOperator(hdp, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a Gibbs operator for sampling counts of draws from the common base distribution" +
                    "for HDPs";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(HDPPolyaUrn.class, false),
        };
    };

}
