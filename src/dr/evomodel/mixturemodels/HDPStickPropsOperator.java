package dr.evomodel.mixturemodels;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

public class HDPStickPropsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String HDP_STICK_PROPS_OPERATOR = "hdpStickPropsOperator";
    private HDPPolyaUrn hdp;
    private CompoundParameter tableCounts;
    //private CompoundParameter commonBaseDistWeights;
    private CompoundParameter stickProportions;
    private Parameter groupAssignments;
    private Parameter categoriesParameter;
    private double pathWeight = 1.0;

    public HDPStickPropsOperator(HDPPolyaUrn hdp, double weight) {

        this.hdp = hdp;
        this.tableCounts = hdp.getTableCounts();
        //this.commonBaseDistWeights = hdp.getCommonBaseDistWeights();
        this.stickProportions = hdp.getStickProportions();
        this.categoriesParameter = hdp.getCategoriesParameter();
        this.groupAssignments = hdp.getGroupAssignments();

        setWeight(weight);
    }


    public double doOperation() {

        // Let G and K denote the current number of groups and categories
        // tablesCounts is a compoundParameter consisting of G parameters
        // dimension k of parameter g indicates the the number of tables in group g corresp. to cat. k

        // counts[g][k] is the number of sites corresponding to group g and cat k
        int[][] counts = new int[hdp.maxGroupCount][hdp.maxCategoryCount];

        //double commonMass = hdp.getCommonMass().getParameterValue(0);
        double comBaseDistMass = hdp.getCommonBaseDistMass().getParameterValue(0);

        for(int i = 0; i < categoriesParameter.getSize(); i++){
            counts[(int)groupAssignments.getParameterValue(i)][(int)categoriesParameter.getParameterValue(i)]++;
        }

        int[] tableCountsForEachCat = new int[hdp.maxCategoryCount];
        int[] countsForEachCat = new int[hdp.maxCategoryCount];
        int[] tablesWithGreaterCatIndex = new int[hdp.maxCategoryCount];

        for(int c = 0; c < hdp.maxCategoryCount; c++){
            for(int g = 0; g < hdp.maxGroupCount; g++){
                tableCountsForEachCat[c] = tableCountsForEachCat[c] + (int) tableCounts.getParameter(g).getParameterValue(c);
                countsForEachCat[c] = countsForEachCat[c] + counts[g][c];
            }
        }


        for(int k = 0; k < hdp.maxCategoryCount; k++){
            for(int i = k+1; i < hdp.maxCategoryCount; i++){
                tablesWithGreaterCatIndex[k] = tablesWithGreaterCatIndex[k] + tableCountsForEachCat[i];
            }
        }

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            if(countsForEachCat[k] > 0){
                double newVal = MathUtils.nextBeta(1 + tableCountsForEachCat[k], comBaseDistMass + tablesWithGreaterCatIndex[k]);
                // we are interested in simulating from the version of the Beta distribution with support (0,1)
                while(newVal == 0 || newVal == 1){
                    newVal = MathUtils.nextBeta(1 + tableCountsForEachCat[k], comBaseDistMass + tablesWithGreaterCatIndex[k]);
                }

                //System.err.println("newVal: " + newVal);
                stickProportions.getParameter(k).setParameterValue(0, newVal);
            }else{
                stickProportions.getParameter(k).setParameterValue(0,0);
            }
        }

        /*
        System.out.println("before returning: ");
        for(int k = 0; k < hdp.maxCategoryCount; k++){
            System.out.println("stick " + k + " is: " + stickProportions.getParameter(k).getParameterValue(0));
        }
        */

        return 0.0;
    }

    /*
    public double[] nextRandom(double[] concParams) {
        double[] randomGammas = new double[concParams.length];
        double sumRandomGammas = 0;
        for(int i = 0; i < concParams.length; i++){
            randomGammas[i] = MathUtils.nextGamma(concParams[i],1);
            sumRandomGammas = sumRandomGammas + randomGammas[i];
        }

        double[] returnValue = new double[concParams.length];
        for(int i = 0; i < concParams.length; i++){
            returnValue[i] = randomGammas[i]/sumRandomGammas;
        }
        return returnValue;
    }
    */

    public void setPathParameter(double beta) {
        if (beta < 0 || beta > 1) {
            throw new IllegalArgumentException("Illegal path weight of " + beta);
        }
        pathWeight = beta;
    }

    public String getOperatorName() {
        return HDP_STICK_PROPS_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return HDP_STICK_PROPS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            HDPPolyaUrn hdp = (HDPPolyaUrn)xo.getChild(HDPPolyaUrn.class);

            return new HDPStickPropsOperator(hdp, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a Gibbs operator for sampling weights corresponding to the common base distribution " +
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
