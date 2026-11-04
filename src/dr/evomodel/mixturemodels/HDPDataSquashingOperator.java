package dr.evomodel.mixturemodels;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.mixturemodels.HDPDataSquashingOperatorParser;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import org.apache.commons.math.MathException;
import java.util.*;

public class HDPDataSquashingOperator extends SimpleMCMCOperator {

    private static final boolean DEBUG = false;
    private static final boolean DEBUG_DIST = false;
    private static final boolean DEBUG_NEW = true;
    private boolean newCatMethodOne = true;
    private boolean newCatMethodTwo = false;
    private HDPPolyaUrn hdp;
    private PatternList patternList;
    private boolean orderOneAdditions = true;
    public static final double ACCURACY_THRESHOLD = 1E-12;

    private int realizationCount;
    private int mhSteps;
    private int uniquePatternCount;

    private Parameter groupAssignments;
    private Parameter categoriesParameter;
    private CountableRealizationsParameter allParameters;
    private CompoundParameter uniquelyRealizedParameters;
    private CompoundParameter tableCounts;
    //private CompoundParameter commonBaseDistWeights;
    private CompoundParameter stickProportions;
    private List<Parameter> massParameterList ;
    private Parameter commonBaseDistMass;
    private List<ParametricMultivariateDistributionModel> commonBaseDistBaseDist;
    private TreeDataLikelihood tdl;
    private CompoundLikelihood cl;
    private int order;
    private double epsilon;
    private double sampleProportion;
    private boolean isIHMM;
    private Parameter currentGroupAssignments;
    private Parameter currentCategoriesParameter;
    private CompoundParameter currentUniquelyRealizedParameters;
    private CompoundParameter currentTableCounts;
    //private CompoundParameter currentCommonBaseDistWeights;
    private int[] sitesStillNeedingUpdate;
    private int numSitesStillNeedingUpdate;
    private boolean cyclical;
    private int distMethod;
    private int fixedNumber;
    private boolean strictCutoff;
    private int parameterDimension;
    // BEGIN ADDED
    private int maxNewCat;
    // END ADDED
    private boolean old;
    private List<SiteRateModel> siteRateModelList;
    private Tree treeModel;

    public HDPDataSquashingOperator(HDPPolyaUrn hdp,
                                    TreeDataLikelihood tdl,
                                    CompoundLikelihood cl,
                                    List<SiteRateModel> siteRateModelList,
                                    Tree treeModel,
                                    PatternList patternList,
                                    int mhSteps,
                                    double weight,
                                    boolean cyclical,
                                    int distMethod,
                                    double eps,
                                    double sampleProp,
                                    int fixedNumber,
                                    boolean strictCutoff,
                                    int maxNewCat,
                                    boolean old
    ) {
        this.old = old;

        this.hdp = hdp;

        // g_i, i=0,...,n-1
        this.groupAssignments = hdp.getGroupAssignments();

        // n
        this.realizationCount = this.groupAssignments.getDimension();

        this.categoriesParameter = hdp.getCategoriesParameter();

        this.tableCounts = hdp.getTableCounts();

        //this.commonBaseDistWeights = hdp.getCommonBaseDistWeights();

        this.stickProportions = hdp.getStickProportions();

        this.uniquelyRealizedParameters = hdp.getUniquelyRealizedParameters();

        this.allParameters = hdp.getAllParameters();

        // M_g for g in G
        this.massParameterList = hdp.getMassParameterList();

        this.commonBaseDistMass = hdp.getCommonBaseDistMass();

        this.commonBaseDistBaseDist = hdp.getCommonBaseDistBaseDist();

        // 0 or 1
        //this.order = gpupp.getOrder();
        this.order = hdp.getOrder();

        this.isIHMM = hdp.isIHMM();

        this.sampleProportion = sampleProp;
        //System.err.println("sampleProportion: " + sampleProportion);
        this.epsilon = eps;
        //System.err.println("epsilon: " + epsilon);
        this.distMethod = distMethod;
        //System.err.println("distMethod: " + distMethod);
        this.fixedNumber = fixedNumber;
        //System.err.println("fixedNumber: " + fixedNumber);
        this.cyclical = cyclical;
        //System.err.println("cyclical: " + cyclical);
        this.strictCutoff = strictCutoff;
        //System.err.println("strictCutoff: " + strictCutoff);

        this.tdl = tdl;
        this.cl = cl;
        this.siteRateModelList = siteRateModelList;
        this.treeModel = treeModel;
        this.patternList = patternList;

        this.uniquePatternCount = patternList.getPatternCount();

        this.mhSteps = mhSteps;

        currentGroupAssignments = new Parameter.Default(groupAssignments.getParameterValues());
        currentCategoriesParameter = new Parameter.Default(categoriesParameter.getParameterValues());
        currentUniquelyRealizedParameters = new CompoundParameter("current");
        for(int k = 0; k < uniquelyRealizedParameters.getParameterCount();k++){
            currentUniquelyRealizedParameters.addParameter(new Parameter.Default(uniquelyRealizedParameters.getParameter(k).getParameterValues()));
        }

        currentTableCounts = new CompoundParameter("currentTableCounts");
        //currentCommonBaseDistWeights = new CompoundParameter("currentCommonBaseDistWeights");
        for(int g = 0; g < hdp.maxGroupCount; g++){
            currentTableCounts.addParameter(new Parameter.Default(tableCounts.getParameter(g).getParameterValues()));
            //currentCommonBaseDistWeights.addParameter(new Parameter.Default(commonBaseDistWeights.getParameter(g).getParameterValues()));
        }

        sitesStillNeedingUpdate = new int[realizationCount];

        if(cyclical) {
            for (int i = 0; i < realizationCount; i++) {
                sitesStillNeedingUpdate[i] = 1;
            }
        }

        numSitesStillNeedingUpdate = realizationCount;

        parameterDimension = uniquelyRealizedParameters.getParameter(0).getSize();

        // BEGIN ADDED
        this.maxNewCat = maxNewCat;
        // END ADDED

        setWeight(weight);
    }

    @Override
    public double doOperation() {
        double returnValue = 0.0;
        try {
            returnValue = doOp();
        } catch (MathException e) {
            e.printStackTrace();
        }
        return returnValue;
    }

    private double doOp() throws MathException {

        for(int i = 0; i < categoriesParameter.getSize(); i++) {
            currentGroupAssignments.setParameterValue(i, groupAssignments.getParameterValue(i));
            currentCategoriesParameter.setParameterValue(i, categoriesParameter.getParameterValue(i));
        }

        // Under this implementation, uniquelyRealizedParameters should have the same parameter count
        // as the max category count
        if(uniquelyRealizedParameters.getParameterCount() != hdp.maxCategoryCount){
            throw new RuntimeException("uniquelyRealizedParameters should have same number of parameters" +
                    "as maxCategoryCount!");
        }

        for(int k = 0; k < uniquelyRealizedParameters.getParameterCount(); k++){
            for(int l = 0; l < uniquelyRealizedParameters.getParameter(k).getSize(); l++){
                currentUniquelyRealizedParameters.getParameter(k).setParameterValue(l,
                        uniquelyRealizedParameters.getParameter(k).getParameterValue(l));
            }
        }

        for(int k = 0; k < hdp.maxGroupCount; k++){
            //currentCommonBaseDistWeights.getParameter(k).setParameterValue(0, commonBaseDistWeights.getParameter(k).getParameterValue(0));
            for(int l = 0; l < tableCounts.getParameter(k).getSize(); l++){
                currentTableCounts.getParameter(k).setParameterValue(l,
                        tableCounts.getParameter(k).getParameterValue(l));
            }
        }

        //for(int g = 0; g < hdp.maxGroupCount; g++){
        //    currentTableCounts.addParameter(new Parameter.Default(tableCounts.getParameter(g).getParameterValues()));
        //    proposedTableCounts.addParameter(new Parameter.Default(tableCounts.getParameter(g).getParameterValues()));
        //}

        if(DEBUG){
            System.err.println("Number of sites: " + categoriesParameter.getSize());
            System.err.println("uniquelyRealizedParameters.getParameterCount(): " + uniquelyRealizedParameters.getParameterCount());
            System.err.println("currentUniqelyRealizedParameters.getParameterCount(): " + currentUniquelyRealizedParameters.getParameterCount());
            for(int c = 0; c < currentCategoriesParameter.getSize(); c++){
                System.err.println("site " + c + " has category " + currentCategoriesParameter.getParameterValue(c) +
                        " and group " + currentGroupAssignments.getParameterValue(c));
            }

            for(int r = 0; r < hdp.maxGroupCount; r++){
                for(int c = 0; c < hdp.maxCategoryCount; c++){
                    System.err.println("Current Table count for group " + r + " and category " + c
                            + " is: " + currentTableCounts.getParameter(r).getParameterValue(c));
                }
            }

        }


        // n_{ij}, entry (i,j) is number of sites associated with group i and category j
        // K categories that run from 0,...,K-1 across all groups
        // Right now, j runs from 0 to maxCategoryCount-1, so some counts may be zero
        // "Empty" categories will be at end, as a result of moving and relabelling done in this algorithm
        // So counts[i][j] == 0 when j >= currentNumCat
        int counts[][] = new int[hdp.maxGroupCount][currentUniquelyRealizedParameters.getParameterCount()];

        for (int i = 0; i < realizationCount; i++) {
            int group = (int) currentGroupAssignments.getParameterValue(i);
            int category = (int) currentCategoriesParameter.getParameterValue(i);
            counts[group][category]++;
        }

        // Current number of active categories
        int currentNumCat = 0;

        int[] countsForEachCat = new int[currentUniquelyRealizedParameters.getParameterCount()];

        // categories that currently have occupancy of 0
        ArrayList<Integer> emptyCats = new ArrayList<Integer>(hdp.maxCategoryCount);

        // BEGIN ADDED
        int[] currentOccupiedCats = new int[hdp.maxCategoryCount];
        // END ADDED

        double[] currentComBaseDistWeights = new double[hdp.maxCategoryCount+1];

        // initialize w_u = 1
        currentComBaseDistWeights[currentComBaseDistWeights.length-1] = 1;

        for(int j = 0; j < counts[0].length; j++){
            for(int i = 0; i < counts.length; i++){
                countsForEachCat[j] = countsForEachCat[j] + counts[i][j];
            }
            if(countsForEachCat[j] > 0){
                currentNumCat++;
                // BEGIN ADDED
                currentOccupiedCats[j] = 1;
                // END ADDED

                // w_k = \beta_k*w_u
                currentComBaseDistWeights[j] = stickProportions.getParameter(j).getParameterValue(0) * currentComBaseDistWeights[currentComBaseDistWeights.length - 1];
                // w_u = (1-\beta_k)*w_u
                currentComBaseDistWeights[currentComBaseDistWeights.length - 1] = (1 - stickProportions.getParameter(j).getParameterValue(0)) * currentComBaseDistWeights[currentComBaseDistWeights.length - 1];
            } else {
                emptyCats.add(j);
                // BEGIN ADDED
                currentOccupiedCats[j] = 0;
                // END ADDED
            }
        }

        double currentWeightsSum = 0;
        for(int c = 0; c < (hdp.maxCategoryCount+1); c++){
            currentWeightsSum = currentWeightsSum + currentComBaseDistWeights[c];
        }

        if(Math.abs(currentWeightsSum-1) > ACCURACY_THRESHOLD){
            throw new RuntimeException("common base dist weights must sum to 1");
        }


        // BEGIN ADDED
        if(emptyCats.size() < maxNewCat){
            throw new RuntimeException("There are not enough unoccupied categories");
        }
        // END ADDED

        if(DEBUG){
            for(int r = 0; r < hdp.maxGroupCount; r++){
                for(int c = 0; c < hdp.maxCategoryCount; c++){
                    if(counts[r][c] > 0) {
                        System.err.println("Current Table count for group " + r + " and category " + c
                                + " is: " + currentTableCounts.getParameter(r).getParameterValue(c));
                        System.err.println("Current count for group " + r + " and category " + c
                                + " is: " + counts[r][c]);
                    }
                }
            }
        }

        if(DEBUG) {
            System.err.println("currentNumCat: " + currentNumCat);
            for(int k = 0; k < currentUniquelyRealizedParameters.getParameterCount(); k++){
                System.err.println("count for cat " + k + " is: " + countsForEachCat[k]);
            }
            System.err.println("currentUniquelyRealizedParameters for active categories: ");
            for(int k = 0; k < currentNumCat; k++){
                System.err.println("cat: " + k + " parameter value: " + currentUniquelyRealizedParameters.getParameter(k).getParameterValue(0));
            }
        }

        // Fill proposedUniquelyRealizedParameters with parameter values corresponding to occupied
        // categories 0,...,K-1
        //for(int k = 0; k < currentNumCat;k++){
        //    proposedUniquelyRealizedParameters.addParameter(new Parameter.Default(currentUniquelyRealizedParameters.getParameter(k).getParameterValues()));
        //}

        // count of sites in each group
        // n_i = sum_{j=1}^{K-1} n_{ij}
        int[] numSitesInGroup = new int[hdp.maxGroupCount];

        // loop over different groups
        for(int k = 0; k < hdp.maxGroupCount; k++){
            // loop over different categories
            for(int l = 0; l < currentNumCat; l++) {
                numSitesInGroup[k] = numSitesInGroup[k] + counts[k][l];
            }

            if(DEBUG){
                System.err.println("Number of sites in group " + k + " is: " + numSitesInGroup[k]);
            }
        }

        //int[] sitesToUpdate = new int[realizationCount];
        // Don't use arraylist, redo this if it is indeed ok to do without the
        // looping of while(numSitesUpdated < realizationCount)
        ArrayList<Integer> sitesToUpdate = new ArrayList<Integer>(0);

        int numSitesUpdated = 0;

        // 1. Initialize siteToUpdate = {0,...,n-1}.

        if(cyclical) {
            for (int i = 0; i < realizationCount; i++) {
                if (sitesStillNeedingUpdate[i] == 1) {
                    sitesToUpdate.add(i);
                }
            }
        }else{
            for (int i = 0; i < realizationCount; i++) {
                sitesToUpdate.add(i);
            }
        }

        //System.err.println("sitesToUpdate.size(): " + sitesToUpdate.size());

        //numSitesStillNeedingUpdate = numSitesStillNeedingUpdate - sitesToUpdate.size();


        // 2. Compute threshold

        // Compute threshold later to use already computed mass functions
        //double threshold = computeThreshold();

        // 3. Repeat until numSitesUpdated == n

        // commented out for the moment
        //while(numSitesUpdated < realizationCount){

        // 3. a) Get current groupVariables and current paramValues

        // 3. b) randomly choose a site from sitesToUpdate

        // BEGIN ADDED
        if(sitesToUpdate.size() == 0){
            for (int i = 0; i < realizationCount; i++) {
                sitesToUpdate.add(i);
            }
        }
        // END ADDED

        // index of randomly chosen site in terms of order in sitesToUpdate ArrayList
        int index = MathUtils.nextInt(sitesToUpdate.size());
        // actual index of randomly chosen site
        int siteIndex = sitesToUpdate.get(index);

        if(DEBUG){
            System.err.println("randomly chosen site has index: " + siteIndex);
        }

        // 3. c) compute mass functions

        // BEGIN ADDED

        // newCatInProposal[j] == 1 if cat j is currently empty but will be used as one of the potential new
        // categories in the proposal
        int[] newCatInProposal = new int[hdp.maxCategoryCount];

        int currentNumDelegates = 0;
        for(int d = 0; d < cl.getLikelihoodCount(); d++){
            if (((TreeDataLikelihood)cl.getLikelihood(d)).getDataLikelihoodDelegate() != null){
                currentNumDelegates++;
            }else {
                break;
            }
        }


        if(old) {
            int idx = 0;

            // Change values of uniquelyRealizedParameter corresponding to empty categories
            // to do: can we replace emptyCats.size() in line below with something else? Careful about needing empty cats in backward move
            for (int k = 0; k < emptyCats.size(); k++) {
                if (idx < maxNewCat) {
                    newCatInProposal[emptyCats.get(k)] = 1;
                }
                //double[] valueForNewCat = commonBaseDistBaseDist.nextRandom();
                double[] valueForNewCat = new double[parameterDimension];

                int counter = 0;
                for (int b = 0; b < commonBaseDistBaseDist.size(); b++) {
                    double[] distVal = commonBaseDistBaseDist.get(b).nextRandom();
                    for (int d = 0; d < distVal.length; d++) {
                        valueForNewCat[counter] = distVal[d];
                        counter++;
                    }
                }

                for (int l = 0; l < parameterDimension; l++) {
                    uniquelyRealizedParameters.getParameter(emptyCats.get(k)).setParameterValue(l, valueForNewCat[l]);
                }
                idx++;
            }

            tdl.makeDirty();
        }else {

            // Change values of uniquelyRealizedParameter corresponding to empty categories
            for (int k = 0; k < maxNewCat; k++) {
                newCatInProposal[emptyCats.get(k)] = 1;

                double[] valueForNewCat = new double[parameterDimension];

                int counter = 0;
                for (int b = 0; b < commonBaseDistBaseDist.size(); b++) {
                    double[] distVal = commonBaseDistBaseDist.get(b).nextRandom();
                    for (int d = 0; d < distVal.length; d++) {
                        valueForNewCat[counter] = distVal[d];
                        counter++;
                    }
                }

                for (int l = 0; l < parameterDimension; l++) {
                    uniquelyRealizedParameters.getParameter(emptyCats.get(k)).setParameterValue(l, valueForNewCat[l]);
                }

                if(emptyCats.get(k) >= currentNumDelegates){

                    if(emptyCats.get(k) > currentNumDelegates){
                        throw new RuntimeException("next empty category should not be greater than currentNumDelegates");
                    }

                    System.out.println("New data likelihood delegate being created");
                    GammaSiteRateModel srm = (GammaSiteRateModel) siteRateModelList.get(currentNumDelegates);
                    BranchModel bm = new HomogeneousBranchModel(srm.getSubstitutionModel(), null);
                    Parameter newCatParam = new Parameter.Default(1);
                    newCatParam.setParameterValue(0,currentNumDelegates);

                    DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                            treeModel,
                            patternList,
                            bm,
                            srm,
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getUseAmbiguities(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getPreferGPU(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getRescalingScheme(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getDelayRescalingUntilUnderflow(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getPreOrderSettings(),
                            categoriesParameter,
                            newCatParam);

                    ((TreeDataLikelihood)cl.getLikelihood(emptyCats.get(k))).setDataLikelihoodDelegate(dataLikelihoodDelegate);

                    currentNumDelegates++;
                }
                cl.getLikelihood(emptyCats.get(k)).makeDirty();
            }
        }

        // END ADDED

        double[] patternLogLikelihoods = null;
        if(old) {
            patternLogLikelihoods = tdl.getDataLikelihoodDelegate().getSiteLogLikelihoods();
        }

        // entry (i,j) is probability of mass function sitesToUpdate.get(i) assuming value j in {0,...,K-1}
        // Where K is number of unique categories, across all groups
        // mass function needed for each site in sitesToUpdate

        double[][] logMassFunctions = null;
        if(old) {
            logMassFunctions = computeLogMassFunctions(sitesToUpdate, currentNumCat, counts, numSitesInGroup,
                    currentGroupAssignments, currentComBaseDistWeights, currentCategoriesParameter, patternLogLikelihoods);
        }else{
            logMassFunctions = computeLogMassFunctionsNew(sitesToUpdate, currentNumCat, counts, numSitesInGroup,
                    currentGroupAssignments, currentComBaseDistWeights, currentCategoriesParameter);
        }


        // Identify sites delta-similar to site at siteIndex
        ArrayList<Integer> sitesForCurrentUpdatePrelim = new ArrayList<Integer>(sitesToUpdate.size());

        if(distMethod == 1) {

            if(DEBUG_DIST) {
                System.err.println("sitesToUpdate before computeThreshold is called: ");
                for (int s = 0; s < sitesToUpdate.size(); s++) {
                    System.err.println("site: " + sitesToUpdate.get(s));

                    for (int t = 0; t < currentNumCat; t++) {
                        System.err.println("logMassFunctions[" + sitesToUpdate.get(s) + "][" + t + "]: " + logMassFunctions[sitesToUpdate.get(s)][t]);
                    }
                }
            }

            if(sitesToUpdate.size() >= 2) {

                // Compute threshold, using logMassFunctions computed above
                double threshold = computeThreshold(epsilon, sampleProportion, sitesToUpdate, logMassFunctions, currentNumCat);
                //double threshold = 0.0000000000005;

                //System.err.println("threshold: " + threshold);

                // Evaluate |R| difference measures

                double[] diffMeasures = computeHellingerDist(logMassFunctions, siteIndex, sitesToUpdate, currentNumCat);

                if (DEBUG_DIST) {
                    for (int k = 0; k < diffMeasures.length; k++) {
                        System.err.println("diffMeasures[" + k + "]: " + diffMeasures[k]);
                    }
                }

                for (int k = 0; k < sitesToUpdate.size(); k++) {
                    if (diffMeasures[k] <= threshold) {
                        sitesForCurrentUpdatePrelim.add(sitesToUpdate.get(k));
                        if (cyclical && order == 0) {
                            sitesStillNeedingUpdate[sitesToUpdate.get(k)] = 0;
                        }
                    }
                }

                if (DEBUG_DIST) {
                    System.err.println("Preliminary number of sites for current update is " +
                            sitesForCurrentUpdatePrelim.size() + " out of " + categoriesParameter.getSize() + " total sites");
                }
            }else{

                for (int k = 0; k < sitesToUpdate.size(); k++) {
                    sitesForCurrentUpdatePrelim.add(sitesToUpdate.get(k));
                    if (cyclical && order == 0) {
                        sitesStillNeedingUpdate[sitesToUpdate.get(k)] = 0;
                    }
                }
            }

        }

        if(distMethod == 2){

            int numSitesInUpdate = (int) (epsilon*sitesToUpdate.size());

            if(fixedNumber > 0){
                numSitesInUpdate = fixedNumber;
                if(fixedNumber > sitesToUpdate.size()){
                    numSitesInUpdate = sitesToUpdate.size();
                }
            }

            //System.err.println("numSitesInUpdate: " + numSitesInUpdate);

            if(numSitesInUpdate == 0){
                // System.err.println("numSitesInUpdate initially 0");

                numSitesInUpdate = sitesToUpdate.size();
                //System.err.println("numSitesInUpdate after adjustment: " + numSitesInUpdate);
            }

            if(currentNumCat > 1) {

                /*
                System.err.println("sitesToUpdate before computeHellingerDist is called: " );
                for(int s = 0; s < sitesToUpdate.size(); s++){
                    System.err.println("site: " + sitesToUpdate.get(s));

                    for(int t = 0; t < currentNumCat; t++){
                        System.err.println("logMassFunctions[" + sitesToUpdate.get(s) + "][" + t + "]: " + logMassFunctions[sitesToUpdate.get(s)][t]);
                    }

                }
                */

                double[] diffMeasures = computeHellingerDist(logMassFunctions, siteIndex, sitesToUpdate, currentNumCat);
                double[] sortedDiffMeasures = new double[diffMeasures.length];
                System.arraycopy(diffMeasures, 0, sortedDiffMeasures, 0, diffMeasures.length);
                Arrays.sort(sortedDiffMeasures);
                double threshold = sortedDiffMeasures[numSitesInUpdate - 1];

                if (DEBUG) {
                    System.err.println("sortedDiffMeasures to compute threshold: ");
                    System.err.println(Arrays.toString(sortedDiffMeasures));
                    System.err.println("threshold index: " + (numSitesInUpdate - 1));
                    System.err.println("threshold: " + threshold);
                    System.err.println("sortedDiffMeasures.length: " + sortedDiffMeasures.length);
                    System.err.println("sitesToUpdate.size(): " + sitesToUpdate.size());
                    //System.err.println("diffMeasures[siteIndex]: " + diffMeasures[siteIndex]);
                }

                if(fixedNumber == 1){
                    sitesForCurrentUpdatePrelim.add(siteIndex);
                }else {

                    if(strictCutoff) {
                        for (int k = 0; k < sitesToUpdate.size(); k++) {
                            if (diffMeasures[k] <= threshold && sitesForCurrentUpdatePrelim.size() < numSitesInUpdate) {
                                sitesForCurrentUpdatePrelim.add(sitesToUpdate.get(k));
                                if (cyclical && order == 0) {
                                    sitesStillNeedingUpdate[sitesToUpdate.get(k)] = 0;
                                }
                            }
                        }
                    }else{

                        if(threshold < 0 ){
                            threshold = 0;
                        }

                        for (int k = 0; k < sitesToUpdate.size(); k++) {
                            if (diffMeasures[k] <= threshold) {
                                sitesForCurrentUpdatePrelim.add(sitesToUpdate.get(k));
                                if (cyclical && order == 0) {
                                    sitesStillNeedingUpdate[sitesToUpdate.get(k)] = 0;
                                }
                            }
                        }
                    }
                }
            } else {

                ArrayList<Integer> shuffledSitesToUpdate = new ArrayList<>(sitesToUpdate.size());

                for (int i = 0; i < sitesToUpdate.size(); i++) {
                    shuffledSitesToUpdate.add(i);
                }

                Random rand = new Random(MathUtils.getSeed());

                Collections.shuffle(shuffledSitesToUpdate, rand);

                //ArrayList<Integer> sampledIndices = new ArrayList<>(numSitesInUpdate);

                if(fixedNumber == 1){
                    sitesForCurrentUpdatePrelim.add(siteIndex);
                }else {
                    for (int i = 0; i < numSitesInUpdate; i++) {
                        sitesForCurrentUpdatePrelim.add(shuffledSitesToUpdate.get(i));
                        if (cyclical && order == 0) {
                            sitesStillNeedingUpdate[shuffledSitesToUpdate.get(i)] = 0;
                        }
                    }
                }
            }
        }

        if(sitesForCurrentUpdatePrelim.size() == 0){
            System.err.println("0 sites in update");
        }

        //System.err.println("number of sites in update: " + sitesForCurrentUpdatePrelim.size());


        // 3. d)
        // for first order processes, remove indices of sites that are not separated
        // from site corresponding to siteIndex

        ArrayList<Integer> sitesForCurrentUpdate = new ArrayList<Integer>(0);

        //System.err.println("siteIndex: " + siteIndex);

        if(order == 1){
            if(cyclical) {
                for (int k = 0; k < sitesForCurrentUpdatePrelim.size(); k++) {
                    //System.err.println("sitesForCurrentUpdatePrelim.get(" + k +"): " + sitesForCurrentUpdatePrelim.get(k));
                    if (Math.abs(sitesForCurrentUpdatePrelim.get(k) - siteIndex) % 2 == 0) {
                        sitesForCurrentUpdate.add(sitesForCurrentUpdatePrelim.get(k));
                        sitesStillNeedingUpdate[sitesForCurrentUpdatePrelim.get(k)] = 0;
                    }
                }
            }else{
                for (int k = 0; k < sitesForCurrentUpdatePrelim.size(); k++) {
                    if (Math.abs(sitesForCurrentUpdatePrelim.get(k) - siteIndex) % 2 == 0) {
                        sitesForCurrentUpdate.add(sitesForCurrentUpdatePrelim.get(k));
                    }
                }
            }
        }else{
            sitesForCurrentUpdate = sitesForCurrentUpdatePrelim;
        }

        if(cyclical) {
            numSitesStillNeedingUpdate = numSitesStillNeedingUpdate - sitesForCurrentUpdate.size();
        }

        //System.err.println("sitesForCurrentUpdate.size(): " + sitesForCurrentUpdate.size());
        //System.err.println("numSitesStillNeedingUpdate: " + numSitesStillNeedingUpdate);

        if(DEBUG){
            System.err.println("Number of sites for current update is " +
                    sitesForCurrentUpdate.size() + " out of " + categoriesParameter.getSize() + " total sites");
            System.err.println("siteIndex: " + siteIndex + " sitesForCurrentUpdate.get(0): " + sitesForCurrentUpdate.get(0));
        }

        // 3. e)
        // compute group variables for sites in sitesForCurrentUpdate
        // right now, current group assignments are computed for all sites at beginning of operation
        // Is this necessary?
        // can we save time by computing group assignments only for sitesForCurrentUpdate?

        // 3. f)
        // Compute (log) allocation variables

        // Compute D_* from D (sitesForCurrentUpdate)
        ArrayList<Integer> sitesForCurrentUpdateStar = new ArrayList<Integer>(0);

        if(order == 1) {
            for (int k = 0; k < sitesForCurrentUpdate.size(); k++) {
                sitesForCurrentUpdateStar.add(sitesForCurrentUpdate.get(k));
                int indexPlusOne = sitesForCurrentUpdate.get(k) + 1;
                if (indexPlusOne < realizationCount) {
                    sitesForCurrentUpdateStar.add(indexPlusOne);
                }
            }
        }else{
            sitesForCurrentUpdateStar = sitesForCurrentUpdate;
        }

        if(DEBUG){
            System.err.println("Number of sites in sitesForCurrentUpdateStar is " +
                    sitesForCurrentUpdateStar.size() + " out of " + categoriesParameter.getSize() + " total sites");
        }

        // K overall categories numbered from 0,...,K-1
        // We adopt convention that alloc var is 0 for categories not associated with sites
        // belonging to group g_{siteIndex}

        // BEGIN ADDED
        double[] logAllocVar = new double [hdp.maxCategoryCount];

        if((currentNumCat+maxNewCat) > hdp.maxCategoryCount){
            throw new RuntimeException("There are not enough empty categories");
        }
        // END ADDED

        // BEGIN REMOVED
        // logAllocVar[K] corresponds to log of what a_{siteIndex,0}^{(D)} is in the text
        // for j = 0,...,K-1, logAllocVar[j] corresponds to log of a_{siteIndex,j}^{(D)}
        //double[] logAllocVar = new double[hdp.maxCategoryCount+1];
        // END REMOVED

        // entry (i,j) is number of sites associated with group i and category j, excluding sites belonging to D
        double[][] countsModified = new double[hdp.maxGroupCount][hdp.maxCategoryCount];

        // entry (i,j) is number of sites associated with group i and category j, excluding sites belonging to D_*
        double[][] countsModifiedStar = new double[hdp.maxGroupCount][hdp.maxCategoryCount];

        // entry k number of sites associated with group k, excluding sites belonging to D_*
        double[] sumCountsModifiedStar = new double[hdp.maxGroupCount];

        // g_{siteIndex}
        int groupForSiteIndex = (int) currentGroupAssignments.getParameterValue(siteIndex);

        // g_{siteIndex+1}
        int groupForSiteAfterSiteIndex;

        for(int k = 0; k < realizationCount; k++){
            int groupNum = (int) currentGroupAssignments.getParameterValue(k);
            int catNum = (int) currentCategoriesParameter.getParameterValue(k);

            // check if site is not in D
            if(!sitesForCurrentUpdate.contains(k)){
                countsModified[groupNum][catNum]++;

                // check additionally if site is not in D_*
                if(!sitesForCurrentUpdateStar.contains(k)){
                    countsModifiedStar[groupNum][catNum]++;
                    sumCountsModifiedStar[groupNum]++;
                }
            }

        }

        if(DEBUG){
            for(int g = 0; g < hdp.maxGroupCount; g++){
                for(int k = 0; k < hdp.maxCategoryCount; k++){
                    System.err.println("counts for group " + g + " and category " + k + " is: " + counts[g][k]);
                    System.err.println("countsModified for group " + g + " and category " + k + " is: " + countsModified[g][k]);
                    System.err.println("countsModifiedStar for group " + g + " and category " + k + " is: " + countsModifiedStar[g][k]);
                }
            }
        }

        int indexWithCategoryOffset;

        // store logLikelihood values so we don't need to recompute them later
        // for k=0,...,K-1, index k is P(Y_{siteIndex} | parameters corresponding to cat. k)
        double[] storedLogLikelihoodValues = new double[hdp.maxCategoryCount];

        // alpha_0 or M_0
        double commonMass = massParameterList.get(0).getParameterValue(0);
        // gamma
        double comBaseDistMass = commonBaseDistMass.getParameterValue(0);

        /*
        int[] tablesForEachCat = new int[hdp.maxCategoryCount];
        int numTables = 0;

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            for(int g = 0; g < hdp.maxGroupCount; g++) {

                if(currentTableCounts.getParameter(g).getParameterValue(k) == 0 && counts[g][k] > 0){
                    throw new RuntimeException("no table in group associated with category, but there are sites" +
                            "in the group associated with the category");
                }

                tablesForEachCat[k] = tablesForEachCat[k] + (int) currentTableCounts.getParameter(g).getParameterValue(k);
            }

            if(tablesForEachCat[k] == 0 && countsForEachCat[k]>0){
                throw new RuntimeException("occupied category with no table");
            }

            if(DEBUG){
                System.err.println("tablesForEachCat[" + k + "]: " + tablesForEachCat[k]);
            }
            numTables = numTables + tablesForEachCat[k];
        }

        if(DEBUG){
            System.err.println("numTables: " + numTables);
        }
        */

        // Now we have what is needed to compute alloc var
        for(int k = 0; k < hdp.maxCategoryCount; k++) {

            // BEGIN ADDED
            if(newCatInProposal[k] == 1) {

                if(old) {
                    indexWithCategoryOffset = uniquePatternCount * k + patternList.getPatternIndex(siteIndex);
                    logAllocVar[k] = patternLogLikelihoods[indexWithCategoryOffset];
                }else{
                    logAllocVar[k] = ((TreeDataLikelihood) cl.getLikelihood(k)).getDataLikelihoodDelegate().getSiteLogLikelihoods()[patternList.getPatternIndex(siteIndex)];
                }

                if(logAllocVar[k] == Double.POSITIVE_INFINITY || Double.isNaN(logAllocVar[k])){
                    logAllocVar[k] = Double.NEGATIVE_INFINITY;
                }

                logAllocVar[k] =
                        logAllocVar[k]
                                + Math.log((commonMass*currentComBaseDistWeights[hdp.maxCategoryCount]/maxNewCat)/(commonMass + sumCountsModifiedStar[groupForSiteIndex]));

                if(order == 1 && orderOneAdditions){
                    if((siteIndex + 1) < realizationCount) {

                        int catSiteAfterSiteIndex = (int) currentCategoriesParameter.getParameterValue(siteIndex + 1);

                        // countsModifiedStar[groupForSiteAfterSiteIndex][] should be zero since we are dealing with a new group
                        logAllocVar[k] = logAllocVar[k] + Math.log(
                                currentComBaseDistWeights[catSiteAfterSiteIndex]);
                    }
                }


            }else {
                // END ADDED

                if (currentOccupiedCats[k] == 1) {

                    if(old) {
                        indexWithCategoryOffset = uniquePatternCount * k + patternList.getPatternIndex(siteIndex);

                        logAllocVar[k] = patternLogLikelihoods[indexWithCategoryOffset];
                    }else{
                        logAllocVar[k] = ((TreeDataLikelihood) cl.getLikelihood(k)).getDataLikelihoodDelegate().getSiteLogLikelihoods()[patternList.getPatternIndex(siteIndex)];
                    }

                    //if(Double.isNaN(patternLogLikelihoods[indexWithCategoryOffset]) || patternLogLikelihoods[indexWithCategoryOffset] == Double.POSITIVE_INFINITY){
                    //    System.out.println("patternLogLikelihoods for cat " + k + " is: " + patternLogLikelihoods[indexWithCategoryOffset]);
                    //}

                    if (DEBUG) {
                        System.err.println("Computing allocation variable for category: " + k);
                        System.err.println("patternLogLikelihoods[k] is: " + logAllocVar[k]);
                        System.err.println("commonMass is: " + commonMass + " comBaseDistMass is: " + comBaseDistMass);
                        //System.err.println("tablesForEachCat[k] is " + tablesForEachCat[k] + " numTables is " + numTables);
                        System.err.println("countsModifiedStar for group " + groupForSiteIndex + " and category "
                                + k + " is: " + countsModifiedStar[groupForSiteIndex][k]);
                        System.err.println("currentComBaseDistWeights is: " + currentComBaseDistWeights[k]);
                        System.err.println("sumCountsModifiedStar for group " + groupForSiteIndex + " is: " + sumCountsModifiedStar[groupForSiteIndex]);
                    }

                    logAllocVar[k] = logAllocVar[k] +
                            Math.log((countsModifiedStar[groupForSiteIndex][k] + commonMass*currentComBaseDistWeights[k])/ (commonMass + sumCountsModifiedStar[groupForSiteIndex]));


                    if (order == 1 && orderOneAdditions) {

                        if ((siteIndex + 1) < realizationCount) {
                            //groupForSiteAfterSiteIndex = (int) currentGroupAssignments.getParameterValue(siteIndex + 1);
                            groupForSiteAfterSiteIndex = 1 + k;

                            int d1 = 0;
                            int d2 = 0;
                            if (groupForSiteIndex == groupForSiteAfterSiteIndex) {
                                d1 = 1;
                            }
                            int catSiteIndex = k;
                            int catSiteAfterSiteIndex = (int) currentCategoriesParameter.getParameterValue(siteIndex + 1);
                            if (catSiteIndex == catSiteAfterSiteIndex) {
                                d2 = 1;
                            }

                            logAllocVar[k] = logAllocVar[k] + Math.log(
                                    (commonMass*currentComBaseDistWeights[catSiteAfterSiteIndex] + countsModifiedStar[groupForSiteAfterSiteIndex][catSiteAfterSiteIndex] + d1 * d2) / (commonMass + sumCountsModifiedStar[groupForSiteAfterSiteIndex] + d1));

                        }
                    }

                    if (DEBUG) {
                        System.err.println("logAllocVar[k] for category " + k + " is " + logAllocVar[k]);
                    }
                } else {
                    logAllocVar[k] = Double.NEGATIVE_INFINITY;
                }
                // BEGIN ADDED
            }
            // END ADDED
        }


        // Done computing logAllocVar elements

        // 3. g)
        // Generate iid realizations for set D according to allocation variables

        double maxVal = logAllocVar[0];

        for(int k = 1; k < logAllocVar.length; k++){
            if(logAllocVar[k] > maxVal){
                maxVal = logAllocVar[k];
            }
        }

        double[] allocVar = new double[logAllocVar.length];
        double allocVarNormConst = 0;

        if(DEBUG){
            System.err.println("subtracting max logAllocVar from each logAllocVar value for numerical stability");
        }

        for(int k = 0; k < allocVar.length; k++) {
            if(DEBUG) {
                System.err.println("logAllocVar for category " + k + " is: " + logAllocVar[k]);
                System.err.println("Corrected logAllocVar is: " + (logAllocVar[k] - maxVal));
            }
            allocVar[k] = Math.exp(logAllocVar[k]-maxVal);
            if(DEBUG){
                System.err.println("unnormalized allocVar for category " + k + " is: " + allocVar[k]);
            }
            allocVarNormConst = allocVarNormConst + allocVar[k];
        }

        if(DEBUG) {
            System.err.println("allocVarNormConst: " + allocVarNormConst);
        }

        //Normalize
        for(int k = 0; k < allocVar.length; k++) {
            allocVar[k] = allocVar[k]/allocVarNormConst;
            //if(DEBUG_NEW && k < 5 && allocVar[0]>0.95) {
            //    System.err.println("normalized allocVar for category: " + k + " is: " + allocVar[k]);
            //}
        }

        // realizations assumes values in {0,1,...,K-1} where K is currentNumCat
        int[] realizations = new int[sitesForCurrentUpdate.size()];

        for(int k = 0; k < realizations.length; k++) {
            realizations[k] = MathUtils.randomChoicePDF(allocVar);
            if(DEBUG){
                System.err.println("realizations for site " + sitesForCurrentUpdate.get(k) + " is: " + realizations[k]);
            }
        }

        // 3. h)
        // Inspect the realizations and identify the sites with realizations from step g) equal to K

        // log Hastings ratio
        double hRatio = 0;

        // BEGIN REMOVED
        //ArrayList<Integer> sitesStillNeedingProposals = new ArrayList<Integer>(0);
        // END REMOVED

        //double currentll = tdl.getLogLikelihood();
        //System.err.println("currentll: " + currentll);


        // BEGIN ADDED

        // equal to 1 if occupied new cat in proposal
        int[] newCatInProposalOccupied = new int[hdp.maxCategoryCount];

        for(int k = 0; k < realizations.length; k++) {

            if(newCatInProposal[realizations[k]] == 1){
                newCatInProposalOccupied[realizations[k]] = 1;
            }

            categoriesParameter.setParameterValue(sitesForCurrentUpdate.get(k), realizations[k]);

            // allocVar used for hRatio instead of logAllocVar, because latter is not normalized
            // Also, allocVar is what is actually used to generate new cat assignments
            hRatio = hRatio - Math.log(allocVar[realizations[k]]);
        }

        for(int c = 0; c < hdp.maxCategoryCount; c++){
            if(newCatInProposalOccupied[c] == 1) {
                //hRatio = hRatio - commonBaseDistBaseDist.logPdf(uniquelyRealizedParameters.getParameter(c).getParameterValues());
                hRatio = hRatio - hdp.getLogDensity(uniquelyRealizedParameters.getParameter(c),0);

                // adjust commonBaseDistWeights

                // draw b from Beta(1,commonBastDistMass)
                //double b = MathUtils.nextBeta(1,comBaseDistMass);
                //double oldWeightNewCat = commonBaseDistWeights.getParameter(hdp.maxCategoryCount).getParameterValue(0);
                //commonBaseDistWeights.getParameter(c).setParameterValue(0,b*oldWeightNewCat);
                //commonBaseDistWeights.getParameter(hdp.maxCategoryCount).setParameterValue(0, (1-b)*oldWeightNewCat);

            }
        }

        // END ADDED

        // 3. j)



        // 3. k)
        // for first order processes, update group assignments for sites that need to have them updated
        // these correspond to sites with indices one greater than the indices of sites in current update
        if(order == 1){
            for(int k = 0; k < sitesForCurrentUpdate.size(); k++){
                if((sitesForCurrentUpdate.get(k)+1) < realizationCount){
                    int catPrecedingSite = (int) categoriesParameter.getParameterValue(sitesForCurrentUpdate.get(k));
                    hdp.setGroup(groupAssignments,sitesForCurrentUpdate.get(k)+1,catPrecedingSite);
                }
            }
            if(DEBUG){
                for(int k = 0; k < categoriesParameter.getSize(); k++){
                    System.err.println("site " + k + " has proposed category " + categoriesParameter.getParameterValue(k)
                            + " and proposed group " + groupAssignments.getParameterValue(k));
                }
            }
        }

        int countsProposed[][] = new int[hdp.maxGroupCount][uniquelyRealizedParameters.getParameterCount()];
        int[] countsForEachProposedCat = new int[uniquelyRealizedParameters.getParameterCount()];

        for (int i = 0; i < realizationCount; i++) {
            countsProposed[(int) groupAssignments.getParameterValue(i)][(int) categoriesParameter.getParameterValue(i)]++;
            countsForEachProposedCat[(int) categoriesParameter.getParameterValue(i)]++;
        }

        // Update proposedTableCounts, if necessary, after group assignments have been updated
        // if proposedTableCounts[g][k] > countsProposed[g][k], set proposedTableCounts[g][k] = countsProposed[g][k]
        // if countsProposed[g][k] > 0 and proposedTableCounts[g][k] == 0, set proposedTableCounts[g][k] = 1
        //if(order == 1){
        for(int g = 0; g < hdp.maxGroupCount; g++){
            for(int k = 0; k < uniquelyRealizedParameters.getParameterCount(); k++){
                if(tableCounts.getParameter(g).getParameterValue(k) > countsProposed[g][k]){
                    tableCounts.getParameter(g).setParameterValue(k, countsProposed[g][k]);
                }
                if( (countsProposed[g][k] > 0) && (tableCounts.getParameter(g).getParameterValue(k) == 0) ){
                    tableCounts.getParameter(g).setParameterValue(k,1.0);
                }
            }
        }

        // update the stickProportions and compute proposedCommonBaseDistWeights

        //System.err.println("Before updating commonBaseDistWeights: " + commonBaseDistWeights.getParameter(hdp.maxCategoryCount).getParameterValue(0));

        double[] proposedComBaseDistWeights = new double[hdp.maxCategoryCount+1];

        // initialize w_u = 1
        proposedComBaseDistWeights[proposedComBaseDistWeights.length-1] = 1;

        double sumProposedWeights = 0;

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            // set weights for newly empty categories to zero
            if(countsForEachProposedCat[k] == 0){
                stickProportions.getParameter(k).setParameterValue(0,0);
            }else{
                // if newly occupied category, drawn no stick proportion from Beta(1, comBaseDistMass)
                if(newCatInProposalOccupied[k] == 1){
                    double b = MathUtils.nextBeta(1,comBaseDistMass);
                    // we are interested in simulating from the version of the Beta distribution with support (0,1)
                    while(b == 1 || b == 0){
                        b = MathUtils.nextBeta(1,comBaseDistMass);
                    }
                    stickProportions.getParameter(k).setParameterValue(0,b);
                }else{
                    if(stickProportions.getParameter(k).getParameterValue(0) == 0){
                        throw new RuntimeException("stick proportion corresponding to occupied category should not be zero");
                    }
                }

                //compute proposedComBaseDistWeights
                proposedComBaseDistWeights[k] = stickProportions.getParameter(k).getParameterValue(0)*proposedComBaseDistWeights[proposedComBaseDistWeights.length-1];
                proposedComBaseDistWeights[proposedComBaseDistWeights.length-1] = (1 - stickProportions.getParameter(k).getParameterValue(0))*proposedComBaseDistWeights[proposedComBaseDistWeights.length-1];
            }

            sumProposedWeights = sumProposedWeights + proposedComBaseDistWeights[k];
        }

        sumProposedWeights = sumProposedWeights + proposedComBaseDistWeights[proposedComBaseDistWeights.length-1];

        if(Math.abs(sumProposedWeights - 1) > ACCURACY_THRESHOLD){
            throw new RuntimeException("commonBaseDistWeights must sum to 1");
        }


        if(DEBUG){
            System.err.println("hRatio before accounting for backward move: " + hRatio);
        }

        //double forwardRatio = hRatio;

        //System.err.println("hRatio for forward move: " + Math.exp(-forwardRatio));

        //double proposedll = tdl.getLogLikelihood();
        //System.err.println("proposedll: " + proposedll);

        // 3. l)
        // finish computing hRatio by adding log probability of current value given proposal

        // We need to compute allocation variables for backward move

        // BEGIN ADDED

        int backwardCounts[][] = new int[hdp.maxGroupCount][uniquelyRealizedParameters.getParameterCount()];

        double[] backwardAllocVar = new double[hdp.maxCategoryCount];
        double[] backwardLogAllocVar = new double[hdp.maxCategoryCount];

        // entry (i,j) is number of sites associated with group i and category j, excluding sites belonging to D
        double[][] backwardCountsModified = new double[hdp.maxGroupCount][hdp.maxCategoryCount];

        // entry (i,j) is number of sites associated with group i and category j, excluding sites belonging to D_*
        double[][] backwardCountsModifiedStar = new double[hdp.maxGroupCount][hdp.maxCategoryCount];

        // entry k number of sites associated with group k, excluding sites belonging to D_*
        double[] backwardSumCountsModifiedStar = new double[hdp.maxGroupCount];

        // END ADDED

        // g_{siteIndex}
        int backwardGroupForSiteIndex = (int) groupAssignments.getParameterValue(siteIndex);
        // g_{siteIndex+1}
        int backwardGroupForSiteAfterSiteIndex;

        for(int k = 0; k < realizationCount; k++){
            int gNum = (int) groupAssignments.getParameterValue(k);
            int cNum = (int) categoriesParameter.getParameterValue(k);

            // BEGIN ADDED
            backwardCounts[gNum][cNum]++;
            // END ADDED

            // check if site is not in D
            if(!sitesForCurrentUpdate.contains(k)){
                backwardCountsModified[gNum][cNum]++;
                // check additionally if site is not in D_*
                if(!sitesForCurrentUpdateStar.contains(k)){
                    backwardCountsModifiedStar[gNum][cNum]++;
                    backwardSumCountsModifiedStar[gNum]++;
                }
            }
        }

        // BEGIN ADDED

        // number of occupied categories in proposal
        int proposedNumCat = 0;

        int[] proposedCountsForEachCat = new int[uniquelyRealizedParameters.getParameterCount()];

        // equal to 1 if cat is occupied in proposal, otherwise 0
        int[] proposedOccupiedCats = new int[hdp.maxCategoryCount];

        // list of cats that are empty in proposal
        ArrayList<Integer> proposedEmptyCats = new ArrayList<Integer>(hdp.maxCategoryCount);


        // TO DO: we can probably extract this information in the loop above, so this is unnecessary

        for(int j = 0; j < backwardCounts[0].length; j++){
            for(int i = 0; i < backwardCounts.length; i++){
                proposedCountsForEachCat[j] = proposedCountsForEachCat[j] + backwardCounts[i][j];
            }
            if(proposedCountsForEachCat[j] > 0){
                proposedNumCat++;
                proposedOccupiedCats[j] = 1;
            }else{
                proposedEmptyCats.add(j);
            }
        }

        // END ADDED


        if(DEBUG){
            System.err.println("Now doing backward move computations");

            for(int g = 0; g < hdp.maxGroupCount; g++){
                for(int k = 0; k < hdp.maxCategoryCount; k++){
                    System.err.println("backwardCountsModified for group " + g + " and category " + k + " is: " + backwardCountsModified[g][k]);
                    System.err.println("backwardCountsModifiedStar for group " + g + " and category " + k + " is: " + backwardCountsModifiedStar[g][k]);
                }
            }
        }

        //int backwardIndexWithCategoryOffset;


        // BEGIN REMOVED
        // need to get logLiklihoods after changing parameter values
        //tdl.makeDirty();
        // END REMOVED

        /*
        int[] bTablesForEachCat = new int[uniquelyRealizedParameters.getParameterCount()];
        int bNumTables = 0;

        for(int k = 0; k < uniquelyRealizedParameters.getParameterCount(); k++){
            for(int g = 0; g < hdp.maxGroupCount; g++) {
                bTablesForEachCat[k] = bTablesForEachCat[k] + (int) tableCounts.getParameter(g).getParameterValue(k);
            }
            bNumTables = bNumTables + bTablesForEachCat[k];
            if(DEBUG){
                System.err.println("bTablesForEachCat[" + k +"]: " + bTablesForEachCat[k]);
            }
        }

        if(DEBUG){
            System.err.println("bNumTables: " + bNumTables);
        }
        */


        // BEGIN ADDED

        // equal to 1 if cat is occupied in current value but not in proposal (so "new" with respect to proposal)
        int[] newCatInCurrent = new int[hdp.maxCategoryCount];

        // newPotentialCatInCurrent[k] == 1 if cat k is one of the maxNewCat potential new cats to be used in backward move
        int[] newPotentialCatInCurrent = new int[hdp.maxCategoryCount];

        // number of categories that are occupied in current value but not in proposal
        int numNewCatsInCurrent = 0;

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            if((currentOccupiedCats[k] == 1) && (proposedOccupiedCats[k] != 1)){
                newCatInCurrent[k] = 1;
                // cat k is a newPotentialCatInCurrent if (but not only if) it is not occupied in proposed state, but is occupied in current state
                newPotentialCatInCurrent[k] = 1;
                numNewCatsInCurrent++;
            }
        }

        // backward move needs maxNewCat potential new categories. Some of these potential new categories will
        // correspond to categories that are occupied in the current state but not the proposed state. We
        // need additional potential new categories if the number of categories occupied in the current state but
        // not proposed state (numNewCatsInCurrent) is less than maxNewCat
        if(numNewCatsInCurrent < maxNewCat && maxNewCat <= proposedEmptyCats.size()){
            int numNewCatsInCurrentStillNeeded = maxNewCat-numNewCatsInCurrent;

            //int counterIndex = 0;

            //for(int k = 0; k < proposedEmptyCats.size(); k++){
            //    if(newPotentialCatInCurrent[proposedEmptyCats.get(k)] != 1 && counterIndex < numNewCatsInCurrentStillNeeded){
            //        newPotentialCatInCurrent[proposedEmptyCats.get(k)] = 1;
            //        counterIndex++;
            //    }
            //}


            for(int k = 0; k < numNewCatsInCurrentStillNeeded; k++){
                if(newPotentialCatInCurrent[proposedEmptyCats.get(k)] != 1){
                    newPotentialCatInCurrent[proposedEmptyCats.get(k)] = 1;
                }

                double[] valForNewCat = new double[parameterDimension];

                int count = 0;
                for (int b = 0; b < commonBaseDistBaseDist.size(); b++) {
                    double[] val = commonBaseDistBaseDist.get(b).nextRandom();
                    for (int d = 0; d < val.length; d++) {
                        valForNewCat[count] = val[d];
                        count++;
                    }
                }

                for (int l = 0; l < parameterDimension; l++) {
                    uniquelyRealizedParameters.getParameter(proposedEmptyCats.get(k)).setParameterValue(l, valForNewCat[l]);
                }

                if(proposedEmptyCats.get(k) >= currentNumDelegates){

                    if(proposedEmptyCats.get(k) > currentNumDelegates){
                        throw new RuntimeException("next empty category should not be greater than currentNumDelegates");
                    }

                    System.out.println("New data likelihood delegate being created");
                    GammaSiteRateModel srm = (GammaSiteRateModel) siteRateModelList.get(currentNumDelegates);
                    BranchModel bm = new HomogeneousBranchModel(srm.getSubstitutionModel(), null);
                    Parameter newCatParam = new Parameter.Default(1);
                    newCatParam.setParameterValue(0,currentNumDelegates);

                    DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                            treeModel,
                            patternList,
                            bm,
                            srm,
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getUseAmbiguities(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getPreferGPU(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getRescalingScheme(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getDelayRescalingUntilUnderflow(),
                            ((TreeDataLikelihood)cl.getLikelihood(0)).getDataLikelihoodDelegate().getPreOrderSettings(),
                            categoriesParameter,
                            newCatParam);

                    ((TreeDataLikelihood)cl.getLikelihood(proposedEmptyCats.get(k))).setDataLikelihoodDelegate(dataLikelihoodDelegate);

                    currentNumDelegates++;
                }else{
                    cl.getLikelihood(proposedEmptyCats.get(k)).makeDirty();
                }
            }
        }

        // Impossible to reach current state via backward move given the number of "new" categories
        // so hastings ratio is zero
        if(numNewCatsInCurrent > maxNewCat || maxNewCat > proposedEmptyCats.size()){
            hRatio = Double.NEGATIVE_INFINITY;
        }else {

            for(int k = 0; k < hdp.maxCategoryCount; k++) {

                if(newPotentialCatInCurrent[k] == 1) {

                    //if(k == 0){
                    //    System.err.println("cat 0 in newPotentialCatInCurrent");
                    //}
                    if(old) {
                        indexWithCategoryOffset = uniquePatternCount * k + patternList.getPatternIndex(siteIndex);
                        backwardLogAllocVar[k] = patternLogLikelihoods[indexWithCategoryOffset];
                    }else{
                        //System.out.println("currentNumDelegates: " + currentNumDelegates);
                        //System.out.println("k: " + k);
                        //System.out.println("proposedNumCat: " + proposedNumCat);
                        //System.out.println("currentNumCat: " + currentNumCat);
                        backwardLogAllocVar[k] = ((TreeDataLikelihood) cl.getLikelihood(k)).getDataLikelihoodDelegate().getSiteLogLikelihoods()[patternList.getPatternIndex(siteIndex)];
                    }

                    backwardLogAllocVar[k] =
                            backwardLogAllocVar[k]
                                    + Math.log((commonMass*proposedComBaseDistWeights[hdp.maxCategoryCount]/maxNewCat)/(commonMass + backwardSumCountsModifiedStar[backwardGroupForSiteIndex]));

                    if(order == 1 && orderOneAdditions){
                        if((siteIndex + 1) < realizationCount) {

                            int catSiteAfterSiteIndex = (int) categoriesParameter.getParameterValue(siteIndex+1);

                            backwardLogAllocVar[k] = backwardLogAllocVar[k]
                                    + Math.log(proposedComBaseDistWeights[catSiteAfterSiteIndex]);
                        }
                    }

                    /*
                    if(order == 1 && orderOneAdditions){
                        if((siteIndex +1) < realizationCount) {
                            backwardGroupForSiteAfterSiteIndex = (int) groupAssignments.getParameterValue(siteIndex+1);

                            int bd1 = 0;
                            int bd2 = 0;
                            if(backwardGroupForSiteIndex == backwardGroupForSiteAfterSiteIndex){
                                bd1 = 1;
                            }
                            int catSiteIndex = k;
                            int catSiteAfterSiteIndex = (int) categoriesParameter.getParameterValue(siteIndex+1);
                            if(catSiteIndex == catSiteAfterSiteIndex){
                                bd2 = 1;
                            }

                            backwardLogAllocVar[k] = backwardLogAllocVar[k]
                                    + Math.log(
                                    (backwardCountsModifiedStar[backwardGroupForSiteAfterSiteIndex][catSiteAfterSiteIndex] + bd1*bd2)/(commonMass + backwardSumCountsModifiedStar[backwardGroupForSiteAfterSiteIndex] + bd1)
                                            + (commonMass/(commonMass + backwardSumCountsModifiedStar[backwardGroupForSiteAfterSiteIndex] + bd1))
                                            *(bTablesForEachCat[catSiteAfterSiteIndex]/(bNumTables + comBaseDistMass)));
                        }
                    }
                    */
                }else {

                    if (proposedOccupiedCats[k] == 1) {

                        // if(k == 0){
                        //     System.err.println("cat 0 in proposedOccupiedCats");
                        //}
                        if(old) {
                            indexWithCategoryOffset = uniquePatternCount * k + patternList.getPatternIndex(siteIndex);
                            backwardLogAllocVar[k] = patternLogLikelihoods[indexWithCategoryOffset];
                        }else{
                            backwardLogAllocVar[k] = ((TreeDataLikelihood) cl.getLikelihood(k)).getDataLikelihoodDelegate().getSiteLogLikelihoods()[patternList.getPatternIndex(siteIndex)];
                        }

                        if (DEBUG) {
                            System.err.println("Computing allocation variable for category: " + k);
                            System.err.println("storedLogLikelihoodValues[k] is: " + backwardLogAllocVar[k]);
                            System.err.println("commonMass is: " + commonMass + " comBaseDistMass is: " + comBaseDistMass);
                            //System.err.println("bTablesForEachCat[k] is " + bTablesForEachCat[k] + " bNumTables is " + bNumTables);
                            System.err.println("countsModifiedStar for group " + backwardGroupForSiteIndex + " and category "
                                    + k + " is: " + backwardCountsModifiedStar[backwardGroupForSiteIndex][k]);
                            System.err.println("sumCountsModifiedStar for group " + backwardGroupForSiteIndex + " is: " + backwardSumCountsModifiedStar[backwardGroupForSiteIndex]);
                        }

                        backwardLogAllocVar[k] = backwardLogAllocVar[k] +
                                Math.log((backwardCountsModifiedStar[backwardGroupForSiteIndex][k] + commonMass*proposedComBaseDistWeights[k])/ (commonMass + backwardSumCountsModifiedStar[backwardGroupForSiteIndex]));

                        if (order == 1 && orderOneAdditions) {
                            if ((siteIndex + 1) < realizationCount) {
                                //backwardGroupForSiteAfterSiteIndex = (int) groupAssignments.getParameterValue(siteIndex + 1);
                                backwardGroupForSiteAfterSiteIndex = k + 1;

                                int bd1 = 0;
                                int bd2 = 0;
                                if (backwardGroupForSiteIndex == backwardGroupForSiteAfterSiteIndex) {
                                    bd1 = 1;
                                }
                                int cSiteIndex = k;
                                int cSiteAfterSiteIndex = (int) categoriesParameter.getParameterValue(siteIndex + 1);
                                if (cSiteIndex == cSiteAfterSiteIndex) {
                                    bd2 = 1;
                                }

                                backwardLogAllocVar[k] = backwardLogAllocVar[k] + Math.log(
                                        (commonMass*proposedComBaseDistWeights[cSiteAfterSiteIndex] + backwardCountsModifiedStar[backwardGroupForSiteAfterSiteIndex][cSiteAfterSiteIndex] + bd1 * bd2) / (commonMass + backwardSumCountsModifiedStar[backwardGroupForSiteAfterSiteIndex] + bd1));

                            }
                        }

                    } else {
                        backwardLogAllocVar[k] = Double.NEGATIVE_INFINITY;
                    }
                }

                if(DEBUG){
                    System.err.println("backwardLogAllocVar[k] for category " + k + " is " + backwardLogAllocVar[k]);
                }
            }



            // Done computing backwardLogAllocVar elements

            // Now need to compute and normalize backwardAllocVar elements

            double backwardMaxVal = backwardLogAllocVar[0];

            for(int k = 1; k < backwardLogAllocVar.length; k++){
                if(backwardLogAllocVar[k] > backwardMaxVal){
                    backwardMaxVal = backwardLogAllocVar[k];
                }
            }

            double backwardAllocVarNormConst = 0;

            for(int k = 0; k < backwardAllocVar.length; k++) {
                if(DEBUG) {
                    System.err.println("backwardLogAllocVar for category " + k + " is: " + backwardLogAllocVar[k]);
                    System.err.println("Corrected backwardLogAllocVar is: " + (backwardLogAllocVar[k] - backwardMaxVal));
                }
                backwardAllocVar[k] = Math.exp(backwardLogAllocVar[k]-backwardMaxVal);
                backwardAllocVarNormConst = backwardAllocVarNormConst + backwardAllocVar[k];
                if(DEBUG){
                    System.err.println("unnormalized backwardAllocVar for category " + k + " is: " + backwardAllocVar[k]);
                }
            }

            if(DEBUG) {
                System.err.println("backwardAllocVarNormConst: " + backwardAllocVarNormConst);
            }

            //Normalize
            for(int k = 0; k < backwardAllocVar.length; k++) {
                backwardAllocVar[k] = backwardAllocVar[k]/backwardAllocVarNormConst;
                if(DEBUG) {
                    System.err.println("normalized backwardAllocVar for category: " + k + " is: " + backwardAllocVar[k]);
                }
            }

            int currentCatForSite;

            // This should be identical to newCatInCurrent, but it is here for sanity check
            int[] newCatInCurrentOccupied = new int[hdp.maxCategoryCount];

            for(int k=0; k < sitesForCurrentUpdate.size(); k++){
                currentCatForSite = (int) currentCategoriesParameter.getParameterValue(sitesForCurrentUpdate.get(k));
                hRatio = hRatio + Math.log(backwardAllocVar[currentCatForSite]);

                if(newCatInCurrent[currentCatForSite] == 1){
                    newCatInCurrentOccupied[currentCatForSite] = 1;
                }
            }

            for(int c = 0; c < hdp.maxCategoryCount; c++){
                if(newCatInCurrent[c] == 1) {

                    if(newCatInCurrentOccupied[c] != 1){
                        throw new RuntimeException("problem in keeping track of current categories that are unoccupied in proposal");
                    }

                    //hRatio = hRatio + commonBaseDistBaseDist.logPdf(uniquelyRealizedParameters.getParameter(c).getParameterValues());
                    hRatio = hRatio + hdp.getLogDensity(uniquelyRealizedParameters.getParameter(c),0);
                }
            }


        }

        // END ADDED




        if(cyclical) {
            if (numSitesStillNeedingUpdate == 0) {
                Arrays.fill(sitesStillNeedingUpdate, 1);
                numSitesStillNeedingUpdate = sitesStillNeedingUpdate.length;
            }
        }

        // Done computing hRatio


        if(DEBUG) {
            System.err.println("uniquelyRealizedParameters.getParameterCount(): " + uniquelyRealizedParameters.getParameterCount());
            System.err.println("hRatio: " + hRatio);
        }

        /*
        System.err.println("hRatio for backward move: " + (Math.exp(hRatio)*Math.exp(-forwardRatio)));
        System.err.println("total hRatio: " + Math.exp(hRatio));
        System.err.println("log hRatio: " + hRatio);
        System.err.println("total likelihood ratio: " + Math.exp(proposedll-currentll));
        System.err.println("full mh ratio: " + ((Math.exp(proposedll-currentll))*Math.exp(hRatio)));
        */


        //System.err.println("hRatio: " + hRatio);

        return hRatio;
        // }

    }



    private int drawRandomElementFromSet(ArrayList<Integer> elements){
        int randomInt = MathUtils.nextInt(elements.size());
        return elements.get(randomInt);
    }


    // returns array diffMeasures, where diffMeasures[j] is difference measure delta_{siteIndex,j}
    /*
    private double[] computeHellingerDist(double[][] logMassFunctions, int siteIndex, ArrayList<Integer> sites, int currentNumCat){

        double[] diffMeasures = new double[sites.size()];

        for(int j = 0; j < sites.size(); j++){

            double diff = 0.0;

            for(int k = 0; k < currentNumCat; k++){
                diff = diff + Math.sqrt(Math.exp(logMassFunctions[siteIndex][k])*Math.exp(logMassFunctions[j][k]));
            }
            diffMeasures[j] = 2*(1-diff);
        }

        return diffMeasures;
    }
    */

    // returns array diffMeasures, where diffMeasures[j] is difference measure delta_{siteIndex,sites.get(j)}
    private double[] computeHellingerDist(double[][] logMassFunctions, int siteIndex, ArrayList<Integer> sites, int currentNumCat){

        double[] diffMeasures = new double[sites.size()];

        for(int j = 0; j < sites.size(); j++) {

            if(DEBUG) {
                System.err.println("mass function for site siteIndex == " + siteIndex);
                for (int s = 0; s < currentNumCat; s++) {
                    System.err.println("Math.exp(logMassFunctions[siteIndex][" + s + "]: " + Math.exp(logMassFunctions[siteIndex][s]));
                }

                System.err.println("mass function for site " + sites.get(j));
                for (int s = 0; s < currentNumCat; s++) {
                    System.err.println("Math.exp(logMassFunctions[" + sites.get(j) + "][" + s + "]: " + Math.exp(logMassFunctions[sites.get(j)][s]));
                }
            }

            if (siteIndex == sites.get(j)) {
                diffMeasures[j] = 0;
            } else {

                double diff = 0.0;

                for (int k = 0; k < currentNumCat; k++) {
                    //if (currentNumCat == 1) {
                    //    System.err.println("logMassFunctions[siteIndex][" + k + "]: " + logMassFunctions[siteIndex][k]);
                    //    System.err.println("logMassFunctions[j][" + k + "]): " + logMassFunctions[sites.get(j)][k]);
                    //}
                    //System.err.println("Math.sqrt(Math.exp(logMassFunctions[siteIndex][k]) * Math.exp(logMassFunctions[sites.get(j)][k])): " + Math.sqrt(Math.exp(logMassFunctions[siteIndex][k]) * Math.exp(logMassFunctions[sites.get(j)][k])));
                    diff = diff + Math.sqrt(Math.exp(logMassFunctions[siteIndex][k]) * Math.exp(logMassFunctions[sites.get(j)][k]));
                }

                //System.err.println("diff should be at most 1, but it is: " + diff);
                // if (siteIndex == j) {
                //     System.err.println("Computing Hellinger Dist between site at siteIndex and itself. diff should be 1 and is:" +
                //             diff + " while diffMeasure is: " + 2 * (1 - diff));
                // }

                diffMeasures[j] = 2 * (1 - diff);
            }
        }
        return diffMeasures;
    }


    // returns array with values delta_{sites.get(i), sites.get(j)} for i!=j
    /*
    private double[] getUnsortedHellingerDistances(double[][] logMassFunctions, ArrayList<Integer> sites, int currentNumCat){

        double[] distances = new double[sites.size()*(sites.size()-1)/2];
        int index = 0;

        if(DEBUG){
            System.err.println("Getting unsorted Hellinger Distances from logMassFunctions ");
            System.err.println("sites.size(): " + sites.size() + " currentNumCat: " + currentNumCat);
            System.err.println("distances.length: " + distances.length);

            for(int i = 0; i < sites.size(); i++) {
                System.err.println("site: " + sites.get(i));
                for (int k = 0; k < currentNumCat; k++) {
                    System.err.println("mass function: " + Math.exp(logMassFunctions[sites.get(i)][k]) + " logMassFunction: " + logMassFunctions[sites.get(i)][k]);
                }
            }
        }

        for(int i = 0; i < sites.size(); i++){
            for(int j = i; j < sites.size(); j++){

                if(i!=j) {
                    double dist = 0.0;

                    for (int k = 0; k < currentNumCat; k++) {
                        dist = dist + Math.sqrt(Math.exp(logMassFunctions[sites.get(i)][k]) * Math.exp(logMassFunctions[sites.get(j)][k]));
                    }

                    distances[index] = 2*(1 - dist);
                    index++;
                }
            }
        }

        if(DEBUG){
            System.err.println("unsorted distances to compute threshold: ");
            System.err.println(Arrays.toString(distances));
        }

        return distances;
    }
    */

    // returns array with values delta_{sites.get(i), sites.get(j)} for i!=j
    private double[] getUnsortedHellingerDistances(double[][] logMassFunctions, ArrayList<Integer> sites, int currentNumCat){

        double[] distances = new double[sites.size()*(sites.size()-1)/2];
        int index = 0;

        if(DEBUG_DIST){
            System.err.println("Getting unsorted Hellinger Distances from logMassFunctions ");
            System.err.println("sites.size(): " + sites.size() + " currentNumCat: " + currentNumCat);
            System.err.println("logMassFunctions.length: " + logMassFunctions.length);

            for(int i = 0; i < sites.size(); i++) {
                System.err.println("site: " + sites.get(i));
                for (int k = 0; k < currentNumCat; k++) {
                    System.err.println("mass function: " + Math.exp(logMassFunctions[sites.get(i)][k]) + " logMassFunction: " + logMassFunctions[sites.get(i)][k]);
                }
            }
        }

        for(int i = 0; i < sites.size(); i++){
            for(int j = i+1; j < sites.size(); j++){

                double dist = 0.0;

                for (int k = 0; k < currentNumCat; k++) {
                    dist = dist + Math.sqrt(Math.exp(logMassFunctions[sites.get(i)][k]) * Math.exp(logMassFunctions[sites.get(j)][k]));
                }

                distances[index] = 2*(1 - dist);

                index++;
            }
        }

        if(DEBUG){
            System.err.println("unsorted distances to compute threshold: ");
            System.err.println(Arrays.toString(distances));
        }

        return distances;
    }

    private double[][] computeLogMassFunctions(ArrayList<Integer> sitesToUpdate,
                                               int currentNumCat,
                                               int[][] counts,
                                               int[] numSitesInGroup,
                                               Parameter currentGroupAssign,
                                               //CompoundParameter tblCts,
                                               double[] weights,
                                               Parameter currentCatParam,
                                               double[] patternLogLikelihoods){

        double[][] logMassFunctions = new double[realizationCount][currentNumCat];

        double commonMass = massParameterList.get(0).getParameterValue(0);
        //double comBaseDistMass = commonBaseDistMass.getParameterValue(0);

        /*
        int[] tablesForEachCat = new int[hdp.maxCategoryCount];
        int numTables = 0;

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            for(int g = 0; g < hdp.maxGroupCount; g++) {
                tablesForEachCat[k] = (int) tblCts.getParameter(g).getParameterValue(k);
            }
            numTables = numTables + tablesForEachCat[k];
        }
        */

        for(int i = 0; i < sitesToUpdate.size(); i++) {

            double normalizingConst = 0.0;

            int activeCatCounter = 0;

            for (int j = 0; j < hdp.maxCategoryCount; j++) {

                int groupNum = (int) currentGroupAssign.getParameterValue(sitesToUpdate.get(i));

                if (counts[groupNum][j] != 0) {
                    int indexWithCatOffset = uniquePatternCount*j + patternList.getPatternIndex(sitesToUpdate.get(i));
                    // System.err.println("indexwithcatoffset: " + indexWithCatOffset);

                    //logMassFunctions[i][j] = tdl.getDataLikelihoodDelegate().getSiteLogLikelihoods()[indexWithCatOffset];
                    logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = patternLogLikelihoods[indexWithCatOffset];

                    // b_g
                    logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = logMassFunctions[sitesToUpdate.get(i)][activeCatCounter]
                            +  Math.log((counts[groupNum][j] + commonMass*weights[j])/(commonMass + numSitesInGroup[groupNum]));

                    if(order == 1 && orderOneAdditions){

                        if((sitesToUpdate.get(i)+1) < realizationCount) {

                            // int groupForSiteAfterSiteIndex = (int) currentGroupAssign.getParameterValue(sitesToUpdate.get(i)+1);
                            int groupForSiteAfterSiteIndex = j + 1;

                            int d1 = 0;
                            int d2 = 0;
                            if(groupNum == groupForSiteAfterSiteIndex){
                                d1 = 1;
                            }
                            int catSiteIndex = j;
                            int catSiteAfterSiteIndex = (int) currentCatParam.getParameterValue(sitesToUpdate.get(i)+1);
                            if(catSiteIndex == catSiteAfterSiteIndex){
                                d2 = 1;
                            }
                            //System.err.println("activeCatCounter: " + activeCatCounter);
                            //System.err.println("sitesToUpdate.get(i): " + sitesToUpdate.get(i));
                            //System.err.println("groupForSiteAfterSiteIndex: " + groupForSiteAfterSiteIndex);
                            //System.err.println("catSiteAfterSiteIndex: " + catSiteAfterSiteIndex);

                            logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] + Math.log(
                                    (commonMass*weights[catSiteAfterSiteIndex] + counts[groupForSiteAfterSiteIndex][catSiteAfterSiteIndex] + d1*d2)/(commonMass + numSitesInGroup[groupForSiteAfterSiteIndex] + d1));

                        }
                    }

                    activeCatCounter++;
                }
            }

            // compute logNormalizingConstant

            // first, subtract logMassFunctions[i][0],...,logMassFunctions[i][currentNumCat-1]
            // by maxValue to ensure numerical stability
            double maxValue = logMassFunctions[sitesToUpdate.get(i)][0];

            for(int j = 1; j < currentNumCat; j++){
                if(logMassFunctions[sitesToUpdate.get(i)][j] > maxValue){
                    maxValue = logMassFunctions[sitesToUpdate.get(i)][j];
                }
            }

            for(int j = 0; j < currentNumCat; j++){
                logMassFunctions[sitesToUpdate.get(i)][j] = logMassFunctions[sitesToUpdate.get(i)][j] - maxValue;
                normalizingConst = normalizingConst + Math.exp(logMassFunctions[sitesToUpdate.get(i)][j]);
            }

            // normalize the entries for mass function i
            for (int j = 0; j < currentNumCat; j++) {
                // System.err.println("log normalizing const: " + Math.log(normalizingConst));
                // System.err.println(" log mass: " + logMassFunctions[i][j]);

                logMassFunctions[sitesToUpdate.get(i)][j] = logMassFunctions[sitesToUpdate.get(i)][j] - Math.log(normalizingConst);
            }
        }
        return logMassFunctions;
    }

    private double[][] computeLogMassFunctionsNew(ArrayList<Integer> sitesToUpdate,
                                                  int currentNumCat,
                                                  int[][] counts,
                                                  int[] numSitesInGroup,
                                                  Parameter currentGroupAssign,
                                                  double[] weights,
                                                  Parameter currentCatParam
                                                  //,
                                                  //double[] patternLogLikelihoods
    ){

        double[][] logMassFunctions = new double[realizationCount][currentNumCat];

        double commonMass = massParameterList.get(0).getParameterValue(0);
        //double comBaseDistMass = commonBaseDistMass.getParameterValue(0);

        /*
        int[] tablesForEachCat = new int[hdp.maxCategoryCount];
        int numTables = 0;

        for(int k = 0; k < hdp.maxCategoryCount; k++){
            for(int g = 0; g < hdp.maxGroupCount; g++) {
                tablesForEachCat[k] = (int) tblCts.getParameter(g).getParameterValue(k);
            }
            numTables = numTables + tablesForEachCat[k];
        }
        */

        for(int i = 0; i < sitesToUpdate.size(); i++) {

            double normalizingConst = 0.0;

            int activeCatCounter = 0;

            for (int j = 0; j < hdp.maxCategoryCount; j++) {

                int groupNum = (int) currentGroupAssign.getParameterValue(sitesToUpdate.get(i));

                if (counts[groupNum][j] != 0) {
                    //int indexWithCatOffset = uniquePatternCount*j + patternList.getPatternIndex(sitesToUpdate.get(i));
                    // System.err.println("indexwithcatoffset: " + indexWithCatOffset);

                    //logMassFunctions[i][j] = tdl.getDataLikelihoodDelegate().getSiteLogLikelihoods()[indexWithCatOffset];
                    //logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = patternLogLikelihoods[indexWithCatOffset];
                    logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = ((TreeDataLikelihood)cl.getLikelihood(j)).getDataLikelihoodDelegate().getSiteLogLikelihoods()[patternList.getPatternIndex(sitesToUpdate.get(i))];

                    // b_g
                    logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = logMassFunctions[sitesToUpdate.get(i)][activeCatCounter]
                            +  Math.log((counts[groupNum][j] + commonMass*weights[j])/(commonMass + numSitesInGroup[groupNum]));

                    if(order == 1 && orderOneAdditions){

                        if((sitesToUpdate.get(i)+1) < realizationCount) {

                            // int groupForSiteAfterSiteIndex = (int) currentGroupAssign.getParameterValue(sitesToUpdate.get(i)+1);
                            int groupForSiteAfterSiteIndex = j + 1;

                            int d1 = 0;
                            int d2 = 0;
                            if(groupNum == groupForSiteAfterSiteIndex){
                                d1 = 1;
                            }
                            int catSiteIndex = j;
                            int catSiteAfterSiteIndex = (int) currentCatParam.getParameterValue(sitesToUpdate.get(i)+1);
                            if(catSiteIndex == catSiteAfterSiteIndex){
                                d2 = 1;
                            }
                            //System.err.println("activeCatCounter: " + activeCatCounter);
                            //System.err.println("sitesToUpdate.get(i): " + sitesToUpdate.get(i));
                            //System.err.println("groupForSiteAfterSiteIndex: " + groupForSiteAfterSiteIndex);
                            //System.err.println("catSiteAfterSiteIndex: " + catSiteAfterSiteIndex);

                            logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] = logMassFunctions[sitesToUpdate.get(i)][activeCatCounter] + Math.log(
                                    (commonMass*weights[catSiteAfterSiteIndex] + counts[groupForSiteAfterSiteIndex][catSiteAfterSiteIndex] + d1*d2)/(commonMass + numSitesInGroup[groupForSiteAfterSiteIndex] + d1));

                        }
                    }

                    activeCatCounter++;
                }
            }

            // compute logNormalizingConstant

            // first, subtract logMassFunctions[i][0],...,logMassFunctions[i][currentNumCat-1]
            // by maxValue to ensure numerical stability
            double maxValue = logMassFunctions[sitesToUpdate.get(i)][0];

            for(int j = 1; j < currentNumCat; j++){
                if(logMassFunctions[sitesToUpdate.get(i)][j] > maxValue){
                    maxValue = logMassFunctions[sitesToUpdate.get(i)][j];
                }
            }

            for(int j = 0; j < currentNumCat; j++){
                logMassFunctions[sitesToUpdate.get(i)][j] = logMassFunctions[sitesToUpdate.get(i)][j] - maxValue;
                normalizingConst = normalizingConst + Math.exp(logMassFunctions[sitesToUpdate.get(i)][j]);
            }

            // normalize the entries for mass function i
            for (int j = 0; j < currentNumCat; j++) {
                // System.err.println("log normalizing const: " + Math.log(normalizingConst));
                // System.err.println(" log mass: " + logMassFunctions[i][j]);

                logMassFunctions[sitesToUpdate.get(i)][j] = logMassFunctions[sitesToUpdate.get(i)][j] - Math.log(normalizingConst);
            }
        }
        return logMassFunctions;
    }

    private double computeThreshold(double epsilon, double sampleProportion, ArrayList<Integer> sitesToUpdate,
                                    double[][] logMassFunctions, int currentNumCat) {

        // First, take a random sample of size sampleProportion*sitesToUpdate.size()
        int sampleSize = (int) (sampleProportion * sitesToUpdate.size());

        if(sampleSize < 2){
            sampleSize = 2;
        }

        if(DEBUG) {
            System.err.println("Now computing threshold");
            System.err.println("sampleSize for computing threshold: " + sampleSize);
        }

        ArrayList<Integer> shuffledSitesToUpdate = new ArrayList<>(sitesToUpdate.size());

        for (int i = 0; i < sitesToUpdate.size(); i++) {
            shuffledSitesToUpdate.add(sitesToUpdate.get(i));
        }

        Random rand = new Random(MathUtils.getSeed());

        Collections.shuffle(shuffledSitesToUpdate, rand);

        ArrayList<Integer> sampledIndices = new ArrayList<>(sampleSize);

        for (int i = 0; i < sampleSize; i++) {
            sampledIndices.add(shuffledSitesToUpdate.get(i));
            if(DEBUG){
                System.err.println("Index " + shuffledSitesToUpdate.get(i) + " is in sampledIndices.");
            }
        }

        double[] distances = getUnsortedHellingerDistances(logMassFunctions, sampledIndices, currentNumCat);

        // Sort entries of distances in ascending order
        Arrays.sort(distances, 0, distances.length - 1);

        if(DEBUG){
            System.err.println("sorted distances to compute threshold: ");
            System.err.println(Arrays.toString(distances));
        }

        // take threshold to be greatest value such that the number of elements of distances[] with
        // distances[i] <= threshold is less than or equal to maxElements
        int maxElements = (int) (epsilon * sampleSize * (sampleSize - 1));
        // first, consider setting threshold = distances[maxElements-1];
        // distances[i] < distances[maxElements-1] for each i <= maxElements-1, because distances[] is sorted
        // but distances[maxElements-1] could be equal to distances[maxElements], which would put the number
        // of elements of distances[] satisfying distances[i] <= threshold at maxElements+1
        // In this case, take threshold = distances[maxIndex], where maxIndex is greatest index such that
        // distances[maxIndex] < distances[maxElements-1]

        if(maxElements == 0){
            maxElements = 1;
        }

        double threshold = distances[maxElements - 1];

        if(DEBUG){
            System.err.println("maxElements: " + maxElements);
            System.err.println("distances[maxElements-1]: " + distances[maxElements-1]);
        }

        if (threshold < 0){
            threshold = 0.0;
            //System.err.println("sorted distances: " + Arrays.toString(distances));
            //throw new RuntimeException("Threshold is 0. Something has gone wrong");
        }
        /*
        else{
            if(distances[maxElements - 1] == distances[maxElements]) {
                int maxIndex = maxElements - 2;
                while((maxIndex > 0) && distances[maxIndex] == distances[maxElements - 1]) {
                    //System.err.println("maxIndex: " + maxIndex + " distances[maxIndex]: " + distances[maxIndex] + " distances[maxElements-1]: " + distances[maxElements-1]);
                    maxIndex = maxIndex - 1;
                }
                threshold = distances[maxIndex];
            }
        }
        */

        if(DEBUG) {
            System.err.println("threshold: " + threshold);
        }
        return threshold;
    }

    @Override
    public String getOperatorName() {
        return HDPDataSquashingOperatorParser.HDP_DATA_SQUASHING_OPERATOR;
    }

}// END: class