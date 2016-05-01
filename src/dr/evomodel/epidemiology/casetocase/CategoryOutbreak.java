/*
 * CategoryOutbreak.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ProductStatistic;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Outbreak class for within-case coalescent.
 *
 * Each case belongs to an infectious (and latent) category which corresponds to one of a list of probability
 * distributions (most likely gamma or exponential) for the length of the infectious (latent) period. The XML rules for
 * the outbreak class ask for at least one ParametricDistributionModel.
 * Assignment of outbreak to distributions should be handled in whatever script or GUI writes the XML.
 *
 * Intended for situations where no data on infection times exists.
 *
 * User: Matthew Hall
 * Date: 17/12/2013
 */

public class CategoryOutbreak extends AbstractOutbreak {

    public static final String CATEGORY_OUTBREAK = "categoryOutbreak";
    private final HashSet<String> latentCategories;
    private final HashSet<String> infectiousCategories;
    private final HashMap<String, Parameter> latentMap;
    private final HashMap<String, AbstractPeriodPriorDistribution> infectiousMap;
    private final HashMap<AbstractCase, Double> weightMap;
    private double[][] distances;


    public CategoryOutbreak(String name, Taxa taxa, boolean hasGeography, boolean hasLatentPeriods,
                            HashMap<String, AbstractPeriodPriorDistribution> infectiousMap,
                            HashMap<String, Parameter> latentMap){
        super(name, taxa, hasLatentPeriods, hasGeography);
        cases = new ArrayList<AbstractCase>();
        latentCategories = new HashSet<String>();
        infectiousCategories = new HashSet<String>();
        this.latentMap = latentMap;
        this.infectiousMap = infectiousMap;
        for(AbstractPeriodPriorDistribution hyperprior : infectiousMap.values()){
            addModel(hyperprior);
        }
        for(Parameter hyperprior : latentMap.values()){
            addVariable(hyperprior);
        }
        weightMap = new HashMap<AbstractCase, Double>();
    }


    private void addCase(String caseID, double endTime, Parameter coords,
                         Parameter infectionPosition, Taxa associatedTaxa, double indexPriorWeight,
                         String infectiousCategory, String latentCategory){
        CategoryCase thisCase;

        if(latentCategory==null){
            thisCase =  new CategoryCase(caseID, endTime, coords, infectionPosition, associatedTaxa,
                    indexPriorWeight, infectiousCategory);
        } else {
            thisCase =
                    new CategoryCase(caseID, endTime, coords, infectionPosition, associatedTaxa,
                            indexPriorWeight, infectiousCategory, latentCategory);
            latentCategories.add(latentCategory);
        }
        weightMap.put(thisCase, indexPriorWeight);

        infectiousCategories.add(infectiousCategory);
        cases.add(thisCase);
        infectedSize++;
        addModel(thisCase);
    }

    private void addNoninfectedCase(String caseID, Parameter coords){
        CategoryCase thisCase = new CategoryCase(caseID, Double.POSITIVE_INFINITY, coords, null, null, 0.0, null);
        thisCase.setEverInfected(false);

        cases.add(thisCase);
        addModel(thisCase);
    }


    public HashMap<AbstractCase, Double> getWeightMap(){
        return weightMap;
    }

    public HashSet<String> getLatentCategories(){
        return latentCategories;
    }

    public HashSet<String> getInfectiousCategories(){
        return infectiousCategories;
    }

    public int getLatentCategoryCount(){
        return latentCategories.size();
    }

    public int getInfectiousCategoryCount(){
        return infectiousCategories.size();
    }

    public HashMap<String, Parameter> getLatentMap(){
        return latentMap;
    }

    public HashMap<String, AbstractPeriodPriorDistribution> getInfectiousMap(){
        return infectiousMap;
    }

    public Parameter getLatentPeriod(String category){
        return latentMap.get(category);
    }

    public AbstractPeriodPriorDistribution getInfectiousCategoryPrior(String category){
        return infectiousMap.get(category);
    }

    public String getInfectiousCategory(AbstractCase aCase){
        return ((CategoryCase)aCase).getInfectiousCategory();
    }

    public String getLatentCategory(AbstractCase aCase){
        return ((CategoryCase)aCase).getLatentCategory();
    }

    public double getLatentPeriod(AbstractCase aCase) {
        return latentMap.get(((CategoryCase)aCase).getLatentCategory()).getParameterValue(0);
    }

    public double getDistance(AbstractCase a, AbstractCase b) {
        if(distances==null){
            throw new RuntimeException("Distance matrix has not been initialised");
        }
        return distances[getCaseIndex(a)][getCaseIndex(b)];
    }

    private void setDistanceMatrix(double[][] distances){
        this.distances = distances;
    }

    private void buildDistanceMatrix(){
        distances = new double[cases.size()][cases.size()];

        if(hasGeography){

            for(int i=0; i<cases.size(); i++){
                for(int j=0; j<cases.size(); j++){
                    distances[i][j]= SpatialKernel.EuclideanDistance(getCase(i).getCoords(), getCase(j).getCoords());
                }
            }
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(!(model instanceof AbstractPeriodPriorDistribution)) {
            fireModelChanged(object);
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
        //nothing to do
    }

    protected void restoreState() {
        //nothing to do
    }

    protected void acceptState() {
        //nothing to do
    }

    private class CategoryCase extends AbstractCase{

        public static final String CATEGORY_CASE = "categoryCase";
        private String infectiousCategory;
        private String latentCategory;
        private Parameter coords;
        private double indexPriorWeight;
        private final ArrayList<Date> examinationTimes;

        private CategoryCase(String name, String caseID, double endTime, Parameter coords,
                             Parameter infectionBranchPosition, Taxa associatedTaxa, double indexPriorWeight,
                             String infectiousCategory, String latentCategory){
            super(name);

            wasEverInfected = associatedTaxa != null;

            this.caseID = caseID;
            this.infectiousCategory = infectiousCategory;
            this.infectionBranchPosition = infectionBranchPosition;
            if(infectionBranchPosition!=null) {
                addVariable(infectionBranchPosition);
            }
            endOfInfectiousTime = endTime;
            this.associatedTaxa = associatedTaxa;
            this.coords = coords;
            this.indexPriorWeight = indexPriorWeight;
            this.latentCategory = latentCategory;

            examinationTimes = new ArrayList<Date>();

            if(wasEverInfected) {
                for (Taxon taxon : associatedTaxa) {
                    examinationTimes.add(taxon.getDate());
                }
            }

        }

        private CategoryCase(String name, String caseID, double endTime, Parameter coords,
                             Parameter infectionBranchPosition, Taxa associatedTaxa,
                             String infectiousCategory){
            this(name, caseID, endTime, coords, infectionBranchPosition, associatedTaxa, 1.0,
                    infectiousCategory, null);

        }



        private CategoryCase(String caseID, double endTime, Parameter coords,
                             Parameter infectionBranchPosition, Taxa associatedTaxa, double indexPriorWeight,
                             String infectiousCategory){
            this(CATEGORY_CASE, caseID, endTime, coords, infectionBranchPosition, associatedTaxa,
                    indexPriorWeight, infectiousCategory, null);
        }


        private CategoryCase(String caseID, double endTime, Parameter coords,
                             Parameter infectionBranchPosition, Taxa associatedTaxa, double indexPriorWeight,
                             String infectiousCategory, String latentCategory){
            this(CATEGORY_CASE, caseID, endTime, coords, infectionBranchPosition, associatedTaxa,
                    indexPriorWeight, infectiousCategory, latentCategory);
        }

        //noninfected susceptible constructor

        public String getLatentCategory(){
            return latentCategory;
        }

        public String getInfectiousCategory(){
            return infectiousCategory;
        }

        public double getIndexPriorWeight() { return indexPriorWeight;}

        public boolean noninfectiousYet(double time) {
            return time>endOfInfectiousTime;
        }

        public ArrayList<Date> getExaminationTimes(){
            return examinationTimes;
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) {
            fireModelChanged();
        }

        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireModelChanged(this);
        }

        protected void storeState() {
            // nothing to do?
        }

        protected void restoreState() {
            // nothing to do?
        }

        protected void acceptState() {
            //nothing to do?
        }

        public double[] getCoords() {
            return new double[]{coords.getParameterValue(0), coords.getParameterValue(1)};
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        //for the outbreak

        public static final String HAS_GEOGRAPHY = "hasGeography";
        public static final String HAS_LATENT_PERIODS = "hasLatentPeriods";
        public static final String INFECTIOUS_PERIOD_PRIOR = "infectiousPeriodPrior";
        public static final String LATENT_PERIODS = "latentPeriods";
        public static final String DISTANCE_MATRIX = "distanceMatrix";

        //for the cases

        public static final String CASE_ID = "hostID";
        public static final String END_TIME = "endTime";
        public static final String COORDINATES = "spatialCoordinates";
        public static final String INFECTION_TIME_BRANCH_POSITION = "infectionTimeBranchPosition";
        public static final String LATENT_CATEGORY = "latentCategory";
        public static final String INFECTIOUS_CATEGORY = "infectiousCategory";
        public static final String WAS_EVER_INFECTED = "wasEverInfected";
        public static final String INDEX_PRIOR_WEIGHT = "indexPriorWeight";

        //for the normal-gamma priors

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final boolean hasGeography = xo.hasAttribute(HAS_GEOGRAPHY) && xo.getBooleanAttribute(HAS_GEOGRAPHY);
            final boolean hasLatentPeriods = Boolean.parseBoolean((String)xo.getAttribute(HAS_LATENT_PERIODS));

            final Taxa taxa = (Taxa) xo.getChild(Taxa.class);

            HashMap<String, AbstractPeriodPriorDistribution> infMap
                    = new HashMap<String, AbstractPeriodPriorDistribution>();

            HashMap<String, Parameter> latMap
                    = new HashMap<String, Parameter>();

            for(int i=0; i<xo.getChildCount(); i++){
                Object cxo = xo.getChild(i);
                if(cxo instanceof XMLObject){
                    if (((XMLObject)cxo).getName().equals(INFECTIOUS_PERIOD_PRIOR)){
                        AbstractPeriodPriorDistribution hyperprior
                                = (AbstractPeriodPriorDistribution)
                                ((XMLObject)cxo).getChild(AbstractPeriodPriorDistribution.class);
                        infMap.put(hyperprior.getModelName(), hyperprior);
                    } else if ((((XMLObject)cxo).getName().equals(LATENT_PERIODS))){
                        Parameter latentPeriod
                                = (Parameter)
                                ((XMLObject)cxo).getChild(Parameter.class);
                        latMap.put(latentPeriod.getParameterName(), latentPeriod);
                    }
                }
            }


            CategoryOutbreak cases = new CategoryOutbreak(null, taxa,
                    hasGeography, hasLatentPeriods, infMap, latMap);
            for(int i=0; i<xo.getChildCount(); i++){
                Object cxo = xo.getChild(i);
                if(cxo instanceof XMLObject && ((XMLObject)cxo).getName()
                        .equals(CategoryCase.CATEGORY_CASE)){
                    parseCase((XMLObject)cxo, cases, hasLatentPeriods);
                }
            }

            if(xo.hasChildNamed(DISTANCE_MATRIX)){
                if(!hasGeography){
                    throw new XMLParseException("Told there is no geography but given a distance matrix");
                }

                Parameter matrixParameter = (Parameter)xo.getElementFirstChild(DISTANCE_MATRIX);

                int size = cases.size();

                if(matrixParameter.getDimension()!=size*size){
                    throw new XMLParseException("Wrong number of distance matrix entries");
                }

                double[][] distances = new double[size][size];
                int count=0;

                for(int i=0; i<size; i++){
                    for(int j=0; j<size; j++){
                        if(i==j){
                            distances[i][j]=0;
                        } else {
                            distances[i][j]=matrixParameter.getParameterValue(count);
                        }
                        count++;
                    }
                }

                cases.setDistanceMatrix(distances);

            } else if(hasGeography){

                for(AbstractCase aCase : cases.getCases()){
                    if(aCase.getCoords()==null){
                        throw new XMLParseException("Some cases have no geographical information");
                    }
                }

                cases.buildDistanceMatrix();
            }


            return cases;
        }

        public void parseCase(XMLObject xo, CategoryOutbreak outbreak, boolean expectLatentPeriods)
                throws XMLParseException {


            String farmID = (String) xo.getAttribute(CASE_ID);

            final Parameter coords = xo.hasChildNamed(COORDINATES) ?
                    (Parameter) xo.getElementFirstChild(COORDINATES) : null;
            boolean wasEverInfected = xo.getBooleanAttribute(WAS_EVER_INFECTED);
            if(wasEverInfected) {

                if(!xo.hasAttribute(INFECTIOUS_CATEGORY)
                        || !xo.hasAttribute(END_TIME)
                        || !xo.hasChildNamed(INFECTION_TIME_BRANCH_POSITION)){
                    throw new XMLParseException("Case " + farmID + " wasEverInfected but lacks infection-related data");
                }


                String infectiousCategory = (String) xo.getAttribute(INFECTIOUS_CATEGORY);
                final double endTime = Double.parseDouble((String) xo.getAttribute(END_TIME));
                String latentCategory = null;
                if (xo.hasAttribute(LATENT_CATEGORY)) {
                    latentCategory = (String) xo.getAttribute(LATENT_CATEGORY);
                } else if (expectLatentPeriods) {
                    throw new XMLParseException("Case " + farmID + " not assigned a latent periods distribution");
                }

                double indexPriorWeight = 1;

                if(xo.hasAttribute(INDEX_PRIOR_WEIGHT)){
                    indexPriorWeight = Double.parseDouble((String)xo.getAttribute(INDEX_PRIOR_WEIGHT));
                }

                final Parameter ibp = (Parameter) xo.getElementFirstChild(INFECTION_TIME_BRANCH_POSITION);


                Taxa taxa = new Taxa();
                for (int i = 0; i < xo.getChildCount(); i++) {
                    if (xo.getChild(i) instanceof Taxon) {
                        taxa.addTaxon((Taxon) xo.getChild(i));
                    }
                }
                outbreak.addCase(farmID, endTime, coords, ibp, taxa, indexPriorWeight, infectiousCategory, latentCategory);
            } else {
                outbreak.addNoninfectedCase(farmID, coords);


            }
        }


        public String getParserDescription(){
            return "Parses a set of 'category' farm outbreak and the information that they all share";
        }

        public Class getReturnType(){
            return CategoryOutbreak.class;
        }

        public String getParserName(){
            return CATEGORY_OUTBREAK;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] caseRules = {
                AttributeRule.newBooleanRule(WAS_EVER_INFECTED),
                new StringAttributeRule(CASE_ID, "The unique identifier for this host"),
                new StringAttributeRule(END_TIME, "The time of noninfectiousness of this host", true),
                new ElementRule(Taxon.class, 0, Integer.MAX_VALUE),
                new ElementRule(INFECTION_TIME_BRANCH_POSITION, Parameter.class, "The exact position on the branch" +
                        " along which the infection of this case occurs that it actually does occur", true),
                new ElementRule(COORDINATES, Parameter.class, "The spatial coordinates of this case", true),
                new StringAttributeRule(LATENT_CATEGORY, "The category of latent period", true),
                new StringAttributeRule(INFECTIOUS_CATEGORY, "The category of infectious period", true),
                new StringAttributeRule(INDEX_PRIOR_WEIGHT, "The weight of this case in the prior probabilty for the" +
                        "index case", true)
        };


        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(HAS_LATENT_PERIODS, "Whether to include a latent period in the model"),
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(CategoryCase.CATEGORY_CASE, caseRules, 1, Integer.MAX_VALUE),
                new ElementRule(Taxa.class),
                new ElementRule(INFECTIOUS_PERIOD_PRIOR, AbstractPeriodPriorDistribution.class, "A prior " +
                        "distribution for the length of infectious periods", 1, Integer.MAX_VALUE),
                new ElementRule(LATENT_PERIODS, Parameter.class, "A prior distribution for the length of latent" +
                        " periods", 0, Integer.MAX_VALUE),
                AttributeRule.newBooleanRule(HAS_GEOGRAPHY, true),
                new ElementRule(DISTANCE_MATRIX, Parameter.class, "A matrix of distances between the cases in this " +
                        "outbreak", true)
        };
    };

}
