package dr.evomodel.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import java.util.Arrays;
import java.util.List;

public class HDPPolyaUrn extends GenPolyaUrnProcessPrior {

    public HDPPolyaUrn(Parameter groupAssignments,
                       Parameter categoriesParameter,
                       CompoundParameter uniquelyRealizedParameters,
                       CountableRealizationsParameter allParameters,
                       List<Parameter> massParameterList,
                       List<ParametricMultivariateDistributionModel> commonBaseDistBaseDist,
                       Parameter commonBaseDistMass,
                       int maxGroups,
                       int maxCat,
                       boolean isIHMM,
                       CompoundParameter tableCounts,
                       //CompoundParameter commonBaseDistWeights,
                       CompoundParameter stickProportions
    ) {
        super(groupAssignments, categoriesParameter, uniquelyRealizedParameters, allParameters, null,
                massParameterList, null, commonBaseDistBaseDist,
                commonBaseDistMass, maxGroups, maxCat, tableCounts, stickProportions, isIHMM, false);
    }


    // set up plausible values for initial tableCounts based upon initial category and group assignments
    public void setUpInitialTableCounts(){

        // entry counts[i][j] corresponds to number of sites corresponding to group i and cat j
        int[][] counts = getCounts();

        for(int g = 0; g < maxGroupCount; g++){
            for(int k = 0; k < maxCategoryCount; k++){
                if(counts[g][k] > 0) {
                    // set initialCount to integer drawn uniformly from 0,1,...,counts[g][k]-1
                    //int initialCount = MathUtils.nextInt(counts[g][k]);
                    // since initialCount must assume a value in {1,2,...,counts[g][k]}
                    // increment initialCount by 1
                    //initialCount = initialCount + 1;
                    // Alternative approach could be to always set initialCount to 1
                    //tableCounts.getParameter(g).setParameterValue(k, initialCount);
                    tableCounts.getParameter(g).setParameterValue(k, 1);
                }else{
                    tableCounts.getParameter(g).setParameterValue(k, 0);
                }
            }
        }
    }



    public void setUpInitialValuesForHDP(){
        // entry counts[i][j] corresponds to number of sites corresponding to group i and cat j
        int[][] counts = getCounts();
        int[] isCatOccupied = new int[maxCategoryCount];

        for(int g = 0; g < maxGroupCount; g++){
            for(int k = 0; k < maxCategoryCount; k++){
                if(counts[g][k] > 0) {
                    isCatOccupied[k] = 1;
                    tableCounts.getParameter(g).setParameterValue(k, 1);
                }else{
                    tableCounts.getParameter(g).setParameterValue(k, 0);
                }
            }
        }

        for(int c = 0; c < maxCategoryCount; c++){
            if(isCatOccupied[c] == 1){
                double simBeta = MathUtils.nextBeta(1,commonBaseDistMass.getParameterValue(0));
                // we are interested in simulating from the version of the Beta distribution with support (0,1)
                while(simBeta == 1 || simBeta == 0){
                    simBeta = MathUtils.nextBeta(1,commonBaseDistMass.getParameterValue(0));
                }

                stickProportions.getParameter(c).setParameterValue(0, simBeta);
            }else{
                stickProportions.getParameter(c).setParameterValue(0,0);
            }
        }


    }

    public void computeGroupAssignments(Parameter catParam) {
        Parameter groupAssign = new Parameter.Default(catParam.getSize());
        if(isIHMM) {
            groupAssignments.setParameterValue(0, 0.0);
            for (int i = 1; i < groupAssign.getSize(); i++) {
                groupAssignments.setParameterValue(i, 1.0 + catParam.getParameterValue(i - 1));
            }
        }else {
            // Do nothing, since groupAssignments is set and fixed a priori
        }
    }


    // may want to rewrite to accommodate HDPs aside from iHMM
    public void setGroup(Parameter gAssignments, int indexToSet, int catPrecedingSite) {
        gAssignments.setParameterValue(indexToSet, 1 + catPrecedingSite);
    }

    public int getOrder() {
        if(isIHMM) {
            return 1;
        }else {
            return 0;
        }
    }

    public boolean isFinite() {
        return false;
    }

    // log density of base distribution for common DP base distribution
    public double getLogDensity(Parameter parameter, int groupNumber) {
        double value[] = parameter.getParameterValues();

        double returnVal = 0;
        int numDistributions = commonBaseDistBaseDist.size();

        int counter = 0;

        for(int i = 0; i < numDistributions; i++){
            int distDim = commonBaseDistBaseDist.get(i).getDimension();
            double[] paramVal = new double[distDim];
            for(int k = 0; k < distDim; k++){
                paramVal[k] = value[counter];
                counter++;
            }
            returnVal = returnVal + commonBaseDistBaseDist.get(i).logPdf(paramVal);
        }

        return returnVal;
    }

    public int getGroupAssignment(int catPrecedingSite){
        System.err.println("getGroupAssignment has been called with catPrecedingSite == " + catPrecedingSite);

        if(getOrder() == 1) {
            return 1 + catPrecedingSite;
        }else{
            // NEED TO CHANGE
            throw new RuntimeException("getGroupAssignment should not be called when groups are fixed");
        }
    }


    public boolean isIHMM() {
        return isIHMM;
    }

    public double computeLogDensity() {

        double loglike = 0;

        //double betacontrib = 0;
        //double paramcontrib = 0;
        //double rest = 0;

        Arrays.fill(isCatActive,0);

        int[][] cts = getCounts();
        double comMass = massParameterList.get(0).getParameterValue(0);
        double comBaseDistMass = commonBaseDistMass.getParameterValue(0);

        //System.err.println("comMass: " + comMass);
        //System.err.println("comBaseDistMass: " + comBaseDistMass);

        int[] numSitesInGroup = new int[maxGroupCount];

        // w_k = \beta_k \prod_{j st j<k} (1-\beta_j), where the \beta are stickProportions
        double[] comBaseDistWeights = new double[maxCategoryCount+1];
        // initialize w_u to 1
        comBaseDistWeights[comBaseDistWeights.length-1] = 1;

        double weightsSum = 0;

        for(int k = 0; k < maxCategoryCount; k++){

            for(int g = 0; g < maxGroupCount; g++){
                if(cts[g][k] > 0){
                    isCatActive[k] = 1;
                }
                numSitesInGroup[g] = numSitesInGroup[g] + cts[g][k];
            }
            // isCatActive[k] is accurate for cat k by this point

            if(isCatActive[k] == 1) {
                // w_k = \beta_k*w_u
                comBaseDistWeights[k] = stickProportions.getParameter(k).getParameterValue(0) * comBaseDistWeights[comBaseDistWeights.length - 1];
                // w_u = (1-\beta_k)*w_u
                comBaseDistWeights[comBaseDistWeights.length - 1] = (1 - stickProportions.getParameter(k).getParameterValue(0)) * comBaseDistWeights[comBaseDistWeights.length - 1];

                // add logDensity of parameter corresponding to cat k
                loglike = loglike + getLogDensity(uniquelyRealizedParameters.getParameter(k),0);
                // paramcontrib = paramcontrib + getLogDensity(uniquelyRealizedParameters.getParameter(k),0);
                // add log pdf of Beta(1,comBaseDistMass) prior for stickProportion
                loglike = loglike + Math.log(comBaseDistMass) + (comBaseDistMass-1)*Math.log(1 - stickProportions.getParameter(k).getParameterValue(0));
                // System.out.println("Math.log(comBaseDistMass): " + Math.log(comBaseDistMass));
                // System.out.println("(comBaseDistMass-1)*Math.log(1 - stickProportions.getParameter(k).getParameterValue(0)): " + (comBaseDistMass-1)*Math.log(1 - stickProportions.getParameter(k).getParameterValue(0)));
                // System.out.println("stickProportions.getParameter(k).getParameterValue(0): " + stickProportions.getParameter(k).getParameterValue(0));
                // betacontrib = betacontrib + Math.log(comBaseDistMass) + (comBaseDistMass-1)*Math.log(1 - stickProportions.getParameter(k).getParameterValue(0));
                //System.out.println("betacontrib: " + betacontrib + " stickprop[" + k + "]:" + stickProportions.getParameter(k).getParameterValue(0));
            }

            weightsSum = weightsSum + comBaseDistWeights[k];
        }

        weightsSum = weightsSum + comBaseDistWeights[comBaseDistWeights.length-1];

        //System.out.println("weightsSum: " + weightsSum);

        if(Math.abs(weightsSum-1) > ACCURACY_THRESHOLD){
            throw new RuntimeException("common base dist weights must sum to 1");
        }

        for(int g = 0; g < maxGroupCount; g++){

            for(int k = 0; k < maxCategoryCount; k++){
                if(cts[g][k] > 0){
                    for (int i = 1; i <= cts[g][k]; i++) {
                        loglike = loglike + Math.log(i - 1 + comMass * comBaseDistWeights[k]);
                        // rest = rest + Math.log(i - 1 + comMass * comBaseDistWeights[k]);
                    }
                }
            }

            //if(g<4) {
            //    System.out.println("numsitesinGroup[" + g + "]: " + numSitesInGroup[g]);
            //}

            for(int s = 1; s <= numSitesInGroup[g]; s++){
                // System.out.println("g: " + g);
                loglike = loglike - Math.log(comMass + s - 1);
                // rest = rest - Math.log(comMass + s - 1);
            }
        }

        isCatActiveKnown = true;

        // System.out.println("loglike: " + loglike + " paramcontrib: " + paramcontrib + " betacontrib: " + betacontrib + " rest: " + rest);
        return loglike;
    }


    public double calculateLogLikelihood() {

        double loglike = computeLogDensity();

        return loglike;
    }

    @Override
    public Model getModel() {
        return this;
    }

}