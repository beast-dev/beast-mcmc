package dr.evomodel.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public abstract class GenPolyaUrnProcessPrior extends AbstractModelLikelihood {

    private static boolean VERBOSE = false;

    public Parameter groupAssignments;
    public Parameter categoriesParameter;
    public CompoundParameter uniquelyRealizedParameters;
    public CountableRealizationsParameter allParameters;
    public CompoundParameter tableCounts;
    //public CompoundParameter commonBaseDistWeights;
    public CompoundParameter stickProportions;
    public List<ParametricMultivariateDistributionModel> baseDistributionList;
    public List<ParametricMultivariateDistributionModel> commonBaseDistBaseDist;
    public List<Parameter> massParameterList;
    public List<Parameter> dpConcentrations;
    public Parameter commonBaseDistMass;

    public int maxCategoryCount;
    public int maxGroupCount;
    public int numStartingCategories;
    public int N;
    //public int[][] counts;
    public boolean likelihoodKnown = false;
    public double logLikelihood;
    public final List<Double> cachedLogFactorials;
    protected boolean isIHMM;
    protected boolean isHMM;
    // isCatActive[k] == 1 iff cat k is active
    protected int[] isCatActive;
    protected boolean isCatActiveKnown = false;
    public static final double ACCURACY_THRESHOLD = 1E-12;


    public GenPolyaUrnProcessPrior(Parameter groupAssign,
                                   Parameter categoriesParam,
                                   CompoundParameter uniquelyRealizedParam,
                                   CountableRealizationsParameter allParam,
                                   List<ParametricMultivariateDistributionModel> baseDistList,
                                   List<Parameter> massParamList,
                                   List<Parameter> dirichletPriorConc,
                                   List<ParametricMultivariateDistributionModel> comBaseDistBaseDist,
                                   Parameter comBaseDistMass,
                                   int maxGroups,
                                   int maxCat,
                                   CompoundParameter tCounts,
                                   //CompoundParameter comBaseDistWeights,
                                   CompoundParameter stickProps,
                                   boolean isIHMM,
                                   boolean isHMM
    ) {

        super("");

        this.isIHMM = isIHMM;

        // category assignment for each site, numbered from 0,...,K-1 where K is total number of
        // distinct categories, irrespective of groups
        this.categoriesParameter = categoriesParam;

        // Parameter j has counts m_{j0},...,m_{j,K-1}
        // m_{jk} is the number of tables in group/urn j corresponding to category k
        tableCounts = tCounts;

        // corresponds to weights w_{0},...,w_{Kmax-1},w_{u}
        //commonBaseDistWeights = comBaseDistWeights;

        // these are the beta parameters, where w_{k} = \beta_{k} \prod_{j<k} (1-beta_{j})
        stickProportions = stickProps;

        // compound parameter of all unique parameters, irrespective of groups
        // phi_0,...,phi_{K-1}
        // runs from 0,...,maxCategoryCount-1, with parameters at indices 0,...,K-1 corresponding to
        // active categories
        this.uniquelyRealizedParameters = uniquelyRealizedParam;

        this.allParameters = allParam;

        // concentration parameters for Dirichlet Prior placed on transition probabilities in the case of HMMs
        // here, it should be null, which is appropriate for Dirichlet Processes
        this.dpConcentrations = dirichletPriorConc;

        //this.categoryCount = uniquelyRealizedParameters.getDimension();
        //this.maxCategoryCount = uniquelyRealizedParameters.getDimension();
        this.maxCategoryCount = maxCat;
        this.maxGroupCount = maxGroups;

        if(maxCategoryCount != uniquelyRealizedParameters.getParameterCount()){
            throw new RuntimeException("cannot have maxCategoryCount unequal to uniquelyRealizedParameters dimension");
        }

        // base distribution for each group
        if(baseDistList == null){
            this.baseDistributionList = setUpBaseDistributionList();
        }else{
            this.baseDistributionList = baseDistList;
        }

        // base distribution for common DP base dist in HDP
        this.commonBaseDistBaseDist = comBaseDistBaseDist;

        // mass for common DP base dist in HDP
        this.commonBaseDistMass = comBaseDistMass;

        // mass parameter for each group
        // M_g, g \in G
        if(massParamList == null){
            this.massParameterList = setUpMassParameterList();
        }else {
            this.massParameterList = massParamList;
        }

        // g_i, group assignment for each of the N sites, i=0,...,N-1
        //this.groupAssignments = groupAssignments;

        // g_i, group assignment for each of the N sites, i=0,...,N-1
        // compute group assignment for each site
        // NEED TO DO: change this to accommodate situations where group assignments are
        // not computed from category assignments, and are instead fixed beforehand
        this.groupAssignments = groupAssign;

        computeGroupAssignments(categoriesParameter);

        //this.groupAssignments = computeGroupAssignments(categoriesParameter);

        //counts = new int[maxGroupCount][maxCategoryCount];
        this.N = groupAssignments.getDimension();

        cachedLogFactorials = new ArrayList<Double>();
        cachedLogFactorials.add(0, 0.0);

        // add all
        this.addVariable(groupAssignments);
        this.addVariable(categoriesParameter);

        if(tableCounts != null && stickProportions != null ) {
            setUpInitialValuesForHDP();
            //setUpInitialTableCounts();
            //System.err.println("initial table counts");
            //for(int i = 0; i < maxGroups; i++){
            //    for(int j = 0; j < maxCat; j++){
            //        System.err.println("tablecount[" + i+ "][" + j+ "]: " +
            //        tableCounts.getParameter(i).getParameterValue(j));
            //    }
            //}

            this.addVariable(tableCounts);
            this.addVariable(stickProportions);
        }



        this.addVariable(uniquelyRealizedParameters);

        if(baseDistributionList != null) {
            for (ParametricMultivariateDistributionModel baseDist : this.baseDistributionList) {
                this.addModel(baseDist);
            }
        }
        for (Parameter massParam : this.massParameterList) {
            this.addVariable(massParam);
        }

        if(commonBaseDistBaseDist != null) {
            for (ParametricMultivariateDistributionModel bDist : this.commonBaseDistBaseDist) {
                this.addModel(bDist);
            }
        }
        if(commonBaseDistMass != null) {
            this.addVariable(commonBaseDistMass);
        }

        this.likelihoodKnown = false;

        isCatActive = new int[categoriesParam.getSize()];

        updateIsCatActive();

        isCatActiveKnown = true;
    }

    // Appropriate for Dirichlet Process, override for other processes
    public List<ParametricMultivariateDistributionModel> setUpBaseDistributionList(){
        return null;
    }

    // Appropriate for Dirichlet Process, override for other processes
    public List<Parameter> setUpMassParameterList(){
        return null;
    }

    public void setUpInitialValuesForHDP(){
        // do nothing
        // override for HDPs
    }

    //public void setUpInitialTableCounts(){
    // do nothing for tableCounts
    // override for HDPs
    //}

    //public void setUpInitialCommonBaseDistWeights(){
    // do nothing
    // override for HDPs
    //}

    public void updateIsCatActive(){
        // entry counts[i][j] corresponds to number of sites corresponding to group i and cat j
        int[][] counts = getCounts();

        Arrays.fill(isCatActive,0);

        for(int g = 0; g < maxGroupCount; g++){
            for(int k = 0; k < maxCategoryCount; k++){
                if(counts[g][k] > 0) {
                    isCatActive[k] = 1;
                }
            }
        }
    }

    // Appropriate for Dirichlet Process, override for other processes
    /*
    public Parameter computeGroupAssignments(Parameter catParam){
        // Default: each site has a group assignment == 1
        Parameter groupAssign = new Parameter.Default(catParam.getSize());
        for(int i = 0; i < groupAssign.getSize(); i++){
            groupAssign.setParameterValue(i,0.0);
        }
        return groupAssign;
    }
    */

    // Appropriate for Dirichlet Process, override for other processes
    public void computeGroupAssignments(Parameter catParam){
        // Default: each site has a group assignment == 1
        for(int i = 0; i < groupAssignments.getSize(); i++){
            groupAssignments.setParameterValue(i,0.0);
        }
    }


    // Appropriate (but never used) for Dirichlet Processes, override for other processes where used
    // here, catPrecedingSite is the category at site indexToSet-1
    public void setGroup(Parameter gAssignments, int indexToSet, int catPrecedingSite){
        gAssignments.setParameterValue(indexToSet,0.0);
    }

    public double getLogFactorial(int i) {
        if ( cachedLogFactorials.size() <= i) {
            for (int j = cachedLogFactorials.size() - 1; j <= i; j++) {
                double logfactorial = cachedLogFactorials.get(j) + Math.log(j + 1);
                cachedLogFactorials.add(logfactorial);
            }
        }
        return cachedLogFactorials.get(i);
    }

    /**
     * Assumes mappings start from index 0. Entry (i,j) corresponds to number of assignments
     * to group i and category j
     * Here, j is overall category number out of categories 0,...,K-1 (not out of 0,...,K_g -1)
     * */
    public int[][] getCounts() {

        //System.err.println("getCounts() has been called");

        //int[][] counts = new int[maxGroupCount][maxCategoryCount];
        //if(!likelihoodKnown) {
        //System.err.println("maxGroupCount: " + maxGroupCount);
        //System.err.println("maxCategoryCount: " + maxCategoryCount);
        int[][] counts = new int[maxGroupCount][maxCategoryCount];
        for (int i = 0; i < N; i++) {
            int group = (int) groupAssignments.getParameterValue(i);
            int category = (int) categoriesParameter.getParameterValue(i);
            counts[group][category]++;

        }
        //  }

        //for(int i = 0; i < 4; i++){
        //    for(int j = 0; j < 4; j++){
        //        System.err.println("counts[" + i + "][" + j + "] inside getCounts(): "
        //        + counts[i][j]);
        //    }
        //}

        //System.err.println("Begin group/cat inside getCounts()");
        //for(int i = 0; i < N; i++){
        //        System.err.println("site " + i + " has cat " + categoriesParameter.getParameterValue(i)
        //        + " and group " + groupAssignments.getParameterValue(i));
        //}
        //System.err.println("End group/cat inside getCounts()");

        return counts;
    }

    // only used for HDPs
    //public int[][] getGroupTableCounts(Parameter gAssign, Parameter tAssign, int maxGCount, int maxTCount) {
    //    return null;
    //}

    public double getMassParam(int groupNumber) {
        return massParameterList.get(groupNumber).getParameterValue(0);
    }

    // Appropriate for DP, overwrite if necessary
    public int getOrder(){
        return 0;
    }

    public boolean isFinite() {
        return false;
    }

    public boolean isHMM() {
        return isHMM;
    }


    // Apppropriate for DP, overwrite if necessary
    public int getGroupAssignment(int catPrecedingSite){
        return 1;
    }

    public Parameter getGroupAssignments() {return groupAssignments;}

    public Parameter getCategoriesParameter() {return categoriesParameter;}

    public CompoundParameter getUniquelyRealizedParameters() {
        return uniquelyRealizedParameters;
    }

    public int[] getIsCatActive(){
        if(isCatActiveKnown) {
            return isCatActive;
        }else{
            updateIsCatActive();
            return isCatActive;
        }
    }

    public CountableRealizationsParameter getAllParameters() { return allParameters; }

    public List<Parameter> getMassParameterList(){return massParameterList;}

    public List<ParametricMultivariateDistributionModel> getBaseDistributionList(){return baseDistributionList;}

    public CompoundParameter getTableCounts() { return tableCounts;}

    public CompoundParameter getStickProportions() { return stickProportions;}

    //public CompoundParameter getCommonBaseDistWeights() { return commonBaseDistWeights;}

    public Parameter getCommonMass() { return massParameterList.get(0);}

    public Parameter getCommonBaseDistMass(){
        return commonBaseDistMass;
    }

    public List<ParametricMultivariateDistributionModel> getCommonBaseDistBaseDist(){
        return commonBaseDistBaseDist;
    }

    public List<ParametricMultivariateDistributionModel> getParametricBaseDist(){
        if(isIHMM){
            return commonBaseDistBaseDist;
        }else{
            return baseDistributionList;
        }
    }

    public double getLogDensity(Parameter parameter, int groupNumber) {
        double value[] = parameter.getParameterValues();

        //System.err.println("baseDistributionList.size(): " + baseDistributionList.size());
        //System.err.println("baseDistribution.get(0).getMean().length: " + baseDistributionList.get(0).getMean().length);
        //System.err.println("baseDistribution.get(1).getMean().length: " + baseDistributionList.get(1).getMean().length);

        double returnVal = 0;
        int numDistributions = baseDistributionList.size();

        int counter = 0;

        for(int i = 0; i < numDistributions; i++){
            int distDim = baseDistributionList.get(i).getDimension();
            //System.err.println("distDim: " + distDim);
            double[] paramVal = new double[distDim];
            for(int k = 0; k < distDim; k++){
                //System.err.println("k: " + k + " counter: " + counter);
                paramVal[k] = value[counter];
                counter++;
            }
            returnVal = returnVal + baseDistributionList.get(i).logPdf(paramVal);
        }

        //return baseDistributionList.get(groupNumber).logPdf(value);
        return returnVal;
    }

    /*
    public double getRealizedValuesLogDensity() {
        double total = 0.0;
        // edit this so that we don't compute counts more times than needed
        int[][] cts = getCounts();
        for (int i = 0; i < maxGroupCount; i++){
            for (int j = 0; j < maxCategoryCount; j++) {
                if (cts[i][j] > 0) {
                    Parameter param = uniquelyRealizedParameters.getParameter(j);
                    total += getLogDensity(param,i);
                }
            }
        }
        System.err.println("realizedValuesLogDensity: " + total);
        return  total;
    }
    */

    public double getSingletonProbability(Parameter param, List<ParametricMultivariateDistributionModel> baseDistList, int groupNum){
        // Singleton probability is 0.0 for DP with continuous base dist
        // override for other models
        return 0.0;
    }

    /*
    public double getCategoriesLogDensity() {

        int[][] cts = getCounts();

        //if (VERBOSE) {
        //    Utils.printArray(counts);
        //}

        // These values correspond to the K_{g_i}
        int[] groupSpecificCatCounts = new int[maxGroupCount];

        double loglike = 0;
        int numSitesInGroup;

        for(int i = 0; i < maxGroupCount; i++) {
            numSitesInGroup = 0;
            groupSpecificCatCounts[i] = 0;
            for (int j = 0; j < maxCategoryCount; j++) {
                int eta = cts[i][j];
                numSitesInGroup = numSitesInGroup + eta;
                if (eta > 0) {
                    //loglike += getLogFactorial(eta - 1);
                    System.err.println("getLogGammaRatio(eta,i,j): " + getLogGammaRatio(eta, i, j));
                    loglike += getLogGammaRatio(eta, i, j);
                    groupSpecificCatCounts[i]++;
                }
            }
            //System.err.println("Math.log(getMassParam(i)): " + Math.log(getMassParam(i)));
            loglike = loglike + groupSpecificCatCounts[i]*Math.log(getMassParam(i));

            for(int k = 1; k <= numSitesInGroup; k++){
                //System.err.println("Math.log(getMassParam(i)+k-1): " + Math.log(getMassParam(i)+k-1));
                loglike = loglike - Math.log(getMassParam(i) + k - 1);
            }
        }
        System.err.println("categoriesLogDensity: " + loglike);
        return loglike;
    }
    */

    // Need to override this for certain models with base distributions assigning nonzero probability
    // to singletons
    public double getLogGammaRatio(int eta, int groupNum, int catNum){
        return getLogFactorial(eta - 1);
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {

        this.fireModelChanged();
        likelihoodKnown = false;
        if (!likelihoodKnown) {

            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public double calculateLogLikelihood() {
        //double loglike = getCategoriesLogDensity() + getRealizedValuesLogDensity();
        //return loglike;

        int[][] cts = getCounts();

        Arrays.fill(isCatActive,0);

        //if (VERBOSE) {
        //    Utils.printArray(counts);
        //}

        // These values correspond to the K_{g_i}
        int[] groupSpecificCatCounts = new int[maxGroupCount];

        double loglike = 0;
        int numSitesInGroup;

        for(int i = 0; i < maxGroupCount; i++) {
            numSitesInGroup = 0;
            groupSpecificCatCounts[i] = 0;
            for (int j = 0; j < maxCategoryCount; j++) {
                //int eta = cts[i][j];
                numSitesInGroup = numSitesInGroup + cts[i][j];
                if (cts[i][j] > 0) {
                    isCatActive[j] = 1;
                    loglike += getLogGammaRatio(cts[i][j],i,j);
                    groupSpecificCatCounts[i]++;

                    Parameter param = uniquelyRealizedParameters.getParameter(j);
                    loglike = loglike + getLogDensity(param,i);
                }
            }
            //System.err.println("Math.log(getMassParam(i)): " + Math.log(getMassParam(i)));
            loglike = loglike + groupSpecificCatCounts[i]*Math.log(getMassParam(i));

            for(int k = 1; k <= numSitesInGroup; k++){
                //System.err.println("Math.log(getMassParam(i)+k-1): " + Math.log(getMassParam(i)+k-1));
                loglike = loglike - Math.log(getMassParam(i) + k - 1);
            }
        }

        //for (int i = 0; i < maxGroupCount; i++){
        //    for (int j = 0; j < maxCategoryCount; j++) {
        //        if (cts[i][j] > 0) {
        //            Parameter param = uniquelyRealizedParameters.getParameter(j);
        //            loglike = loglike + getLogDensity(param,i);
        //        }
        //    }
        //}

        //System.err.println("loglike from dppolyaurn: " + loglike);

        isCatActiveKnown = true;

        return loglike;
    }

    @Override
    public void makeDirty() {
        //likelihoodKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Variable.ChangeType type) {

        if (variable == categoriesParameter) {
            isCatActiveKnown = false;
            this.fireModelChanged();
        } else if (variable == groupAssignments) {
            this.fireModelChanged();
        } else if (variable == tableCounts) {
            this.fireModelChanged();
            //} else if (variable == commonBaseDistWeights) {
            //    this.fireModelChanged();
        }else if (variable == stickProportions){
            this.fireModelChanged();
        } else if (variable == massParameterList) {
            this.fireModelChanged();
        } else if(massParameterList.contains(variable)) {
            this.fireModelChanged();
        } else if(variable == commonBaseDistMass){
            this.fireModelChanged();
        } else if (variable == uniquelyRealizedParameters) {
            likelihoodKnown = false;
            this.fireModelChanged();
        } else {
            System.err.println("variable name: " + variable.getVariableName());
            throw new IllegalArgumentException("Unknown parameter");
        }
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {

    }


}