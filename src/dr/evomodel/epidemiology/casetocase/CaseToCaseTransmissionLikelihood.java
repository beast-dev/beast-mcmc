/*
 * CaseToCaseTransmissionLikelihood.java
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

import dr.app.tools.NexusExporter;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.xml.*;

import java.io.PrintStream;
import java.util.*;

/**
 * A likelihood function for transmission between identified epidemiological outbreak
 *
 * Timescale must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist;
 * contact MH.
 *
 * Latent periods are not implemented currently
 *
 * @author Matthew Hall
 * @version $Id: $
 */

public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable {

    private static final boolean DEBUG = false;

    private CategoryOutbreak outbreak;
    private CaseToCaseTreeLikelihood treeLikelihood;
    private SpatialKernel spatialKernel;
    private Parameter transmissionRate;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;
    private boolean transProbKnown;
    private boolean storedTransProbKnown;
    private boolean periodsProbKnown;
    private boolean storedPeriodsProbKnown;
    private boolean treeProbKnown;
    private boolean storedTreeProbKnown;
    private double logLikelihood;
    private double storedLogLikelihood;

    private double transLogProb;
    private double storedTransLogProb;
    private double periodsLogProb;
    private double storedPeriodsLogProb;
    private double treeLogProb;
    private double storedTreeLogProb;

    private ParametricDistributionModel initialInfectionTimePrior;
    private HashMap<AbstractCase, Double> indexCasePrior;

    private final boolean hasGeography;
    private final boolean hasLatentPeriods;
    private ArrayList<TreeEvent> sortedTreeEvents;
    private ArrayList<TreeEvent> storedSortedTreeEvents;

    private AbstractCase indexCase;
    private AbstractCase storedIndexCase;

//    private F f;

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    public CaseToCaseTransmissionLikelihood(String name, CategoryOutbreak outbreak,
                                            CaseToCaseTreeLikelihood treeLikelihood, SpatialKernel spatialKernal,
                                            Parameter transmissionRate,
                                            ParametricDistributionModel intialInfectionTimePrior){
        super(name);
        this.outbreak = outbreak;
        this.treeLikelihood = treeLikelihood;
        this.spatialKernel = spatialKernal;
        if(spatialKernal!=null){
            this.addModel(spatialKernal);
        }
        this.transmissionRate = transmissionRate;
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        likelihoodKnown = false;
        hasGeography = spatialKernal!=null;
        this.hasLatentPeriods = treeLikelihood.hasLatentPeriods();

        this.initialInfectionTimePrior = intialInfectionTimePrior;


        HashMap<AbstractCase, Double> weightMap = outbreak.getWeightMap();

        double totalWeights = 0;

        for(AbstractCase aCase : weightMap.keySet()){
            if(aCase.wasEverInfected) {
                totalWeights += weightMap.get(aCase);
            }
        }

        indexCasePrior = new HashMap<AbstractCase, Double>();

        for(AbstractCase aCase : outbreak.getCases()){
            if(aCase.wasEverInfected) {
                indexCasePrior.put(aCase, weightMap.get(aCase) / totalWeights);
            }
        }

        sortEvents();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){

            treeProbKnown = false;
            if(!(object instanceof DemographicModel)){
                transProbKnown = false;
                periodsProbKnown = false;
                sortedTreeEvents = null;
                indexCase = null;
            }


        } else if(model instanceof SpatialKernel){

            transProbKnown = false;

        } else if(model instanceof AbstractOutbreak){

            transProbKnown = false;
            periodsProbKnown = false;
            sortedTreeEvents = null;
            indexCase = null;


        }
        likelihoodKnown = false;
    }

    // no need to change the RNG queue unless the normalisation will have changed

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable==transmissionRate){
            transProbKnown = false;
        }
        likelihoodKnown = false;
    }

    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedPeriodsLogProb = periodsLogProb;
        storedPeriodsProbKnown = periodsProbKnown;
        storedTransLogProb = transLogProb;
        storedTransProbKnown = transProbKnown;
        storedTreeLogProb = treeLogProb;
        storedTreeProbKnown = treeProbKnown;
        storedSortedTreeEvents = new ArrayList<TreeEvent>(sortedTreeEvents);
        storedIndexCase = indexCase;
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        transLogProb = storedTransLogProb;
        transProbKnown = storedTransProbKnown;
        treeLogProb = storedTreeLogProb;
        treeProbKnown = storedTreeProbKnown;
        periodsLogProb = storedPeriodsLogProb;
        periodsProbKnown = storedPeriodsProbKnown;
        sortedTreeEvents = storedSortedTreeEvents;
        indexCase = storedIndexCase;
    }

    protected void acceptState() {
        // nothing to do
    }

    public SpatialKernel getSpatialKernel(){
        return spatialKernel;
    }

    public Model getModel() {
        return this;
    }

    public CaseToCaseTreeLikelihood getTreeLikelihood(){
        return treeLikelihood;
    }

    public double getLogLikelihood() {

        if(!likelihoodKnown) {
            if (!treeProbKnown) {
                treeLikelihood.prepareTimings();
            }
            if (!transProbKnown) {

                try {

                    transLogProb = 0;

                    if (sortedTreeEvents == null) {
                        sortEvents();
                    }

                    double rate = transmissionRate.getParameterValue(0);

                    ArrayList<AbstractCase> previouslyInfectious = new ArrayList<AbstractCase>();

                    double currentEventTime;
                    boolean first = true;

                    for (TreeEvent event : sortedTreeEvents) {
                        currentEventTime = event.getTime();

                        AbstractCase thisCase = event.getCase();

                        if (event.getType() == EventType.INFECTION) {
                            if (first) {
                                // index infection

                                if (indexCasePrior != null) {
                                    transLogProb += Math.log(indexCasePrior.get(thisCase));
                                }
                                if (initialInfectionTimePrior != null) {
                                    transLogProb += initialInfectionTimePrior.logPdf(currentEventTime);
                                }

                                if (!hasLatentPeriods) {
                                    previouslyInfectious.add(thisCase);
                                }

                                first = false;

                            } else {

                                AbstractCase infector = event.getInfector();

                                if(thisCase.wasEverInfected()) {


                                    if (previouslyInfectious.contains(thisCase)){
                                        throw new BadPartitionException(thisCase.caseID +
                                                " infected after it was infectious");
                                    }

                                    if (event.getTime() > thisCase.endOfInfectiousTime){
                                        throw new BadPartitionException(thisCase.caseID +
                                                " ceased to be infected before it was infected");
                                    }
                                    if (infector.endOfInfectiousTime < event.getTime()){
                                        throw new BadPartitionException(thisCase.caseID + " infected by "
                                                + infector.caseID + " after the latter ceased to be infectious");
                                    }
                                    if (treeLikelihood.getInfectiousTime(infector) > event.getTime()) {
                                        throw new BadPartitionException(thisCase.caseID + " infected by "
                                                + infector.caseID + " before the latter became infectious");
                                    }

                                    if(!previouslyInfectious.contains(infector)){
                                        throw new RuntimeException("Infector not previously infected");
                                    }
                                }

                                // no other previously infectious case has infected this case...

                                for (AbstractCase nonInfector : previouslyInfectious) {



                                    double timeDuringWhichNoInfection;
                                    if (nonInfector.endOfInfectiousTime < event.getTime()) {
                                        timeDuringWhichNoInfection = nonInfector.endOfInfectiousTime
                                                - treeLikelihood.getInfectiousTime(nonInfector);
                                    } else {
                                        timeDuringWhichNoInfection = event.getTime()
                                                - treeLikelihood.getInfectiousTime(nonInfector);
                                    }

                                    if(timeDuringWhichNoInfection<0){
                                        throw new RuntimeException("negative time");
                                    }

                                    double transRate = rate;
                                    if (hasGeography) {
                                        transRate *= outbreak.getKernelValue(thisCase, nonInfector, spatialKernel);
                                    }

                                    transLogProb += -transRate * timeDuringWhichNoInfection;


                                }

                                // ...until the end

                                if(thisCase.wasEverInfected()) {
                                    double transRate = rate;


                                    if (hasGeography) {
                                        transRate *= outbreak.getKernelValue(thisCase, infector, spatialKernel);
                                    }


                                    transLogProb += Math.log(transRate);
                                }

                                if (!hasLatentPeriods) {
                                    previouslyInfectious.add(thisCase);
                                }


                            }


                        } else if (event.getType() == EventType.INFECTIOUSNESS) {
                            if (event.getTime() < Double.POSITIVE_INFINITY) {

                                if(event.getTime() > event.getCase().endOfInfectiousTime){
                                    throw new BadPartitionException(event.getCase().caseID + " noninfectious before" +
                                            "infectious");
                                }

                                if (first) {
                                    throw new RuntimeException("First event is not an infection");
                                }

                                previouslyInfectious.add(thisCase);
                            }
                        }
                    }

                    transProbKnown = true;
                } catch (BadPartitionException e) {

                    transLogProb = Double.NEGATIVE_INFINITY;
                    transProbKnown = true;
                    logLikelihood = Double.NEGATIVE_INFINITY;
                    likelihoodKnown = true;
                    return logLikelihood;

                }


            }

            if(!periodsProbKnown){

                periodsLogProb = 0;

                HashMap<String, ArrayList<Double>> infectiousPeriodsByCategory
                        = new HashMap<String, ArrayList<Double>>();


                for (AbstractCase aCase : outbreak.getCases()) {
                    if(aCase.wasEverInfected()) {

                        String category = (outbreak).getInfectiousCategory(aCase);

                        if (!infectiousPeriodsByCategory.keySet().contains(category)) {
                            infectiousPeriodsByCategory.put(category, new ArrayList<Double>());
                        }

                        ArrayList<Double> correspondingList
                                = infectiousPeriodsByCategory.get(category);

                        correspondingList.add(treeLikelihood.getInfectiousPeriod(aCase));
                    }
                }


                for (String category : outbreak.getInfectiousCategories()) {

                    Double[] infPeriodsInThisCategory = infectiousPeriodsByCategory.get(category)
                            .toArray(new Double[infectiousPeriodsByCategory.get(category).size()]);

                    AbstractPeriodPriorDistribution hyperprior = outbreak.getInfectiousCategoryPrior(category);

                    double[] values = new double[infPeriodsInThisCategory.length];

                    for (int i = 0; i < infPeriodsInThisCategory.length; i++) {
                        values[i] = infPeriodsInThisCategory[i];
                    }

                    periodsLogProb += hyperprior.getLogLikelihood(values);

                }

                periodsProbKnown = true;

            }


            if(!treeProbKnown){
                treeLogProb = treeLikelihood.getLogLikelihood();
                treeProbKnown = true;
            }

            // just reject states where these round to +INF

            if(transLogProb == Double.POSITIVE_INFINITY){
                System.out.println("TransLogProb +INF");
                return Double.NEGATIVE_INFINITY;
            }
            if(periodsLogProb == Double.POSITIVE_INFINITY){
                System.out.println("PeriodsLogProb +INF");
                return Double.NEGATIVE_INFINITY;
            }
            if(treeLogProb == Double.POSITIVE_INFINITY){
                System.out.println("TreeLogProb +INF");
                return Double.NEGATIVE_INFINITY;
            }

            logLikelihood =  treeLogProb + periodsLogProb + transLogProb;
            likelihoodKnown = true;
        }

        return logLikelihood;
    }


    public void makeDirty() {
        likelihoodKnown = false;
        transProbKnown = false;
        periodsProbKnown = false;
        treeProbKnown = false;
        sortedTreeEvents = null;
        treeLikelihood.makeDirty();
        indexCase = null;
    }

    private class EventComparator implements Comparator<TreeEvent> {
        public int compare(TreeEvent treeEvent1, TreeEvent treeEvent2) {
            return Double.compare(treeEvent1.getTime(),
                    treeEvent2.getTime());
        }
    }


    private enum EventType{
        INFECTION,
        INFECTIOUSNESS,
        END
    }

    private void sortEvents(){
        ArrayList<TreeEvent> out = new ArrayList<TreeEvent>();
        for(AbstractCase aCase : outbreak.getCases()){


            double infectionTime = treeLikelihood.getInfectionTime(aCase);
            out.add(new TreeEvent(infectionTime, aCase, treeLikelihood.getInfector(outbreak.getCaseIndex(aCase))));

            if(aCase.wasEverInfected()) {

                double endTime = aCase.endOfInfectiousTime;

                out.add(new TreeEvent(EventType.END, endTime, aCase));

                if (hasLatentPeriods) {
                    double infectiousnessTime = treeLikelihood.getInfectiousTime(aCase);
                    out.add(new TreeEvent(EventType.INFECTIOUSNESS, infectiousnessTime, aCase));

                }
            }
        }

        Collections.sort(out, new EventComparator());



        indexCase = out.get(0).getCase();

        sortedTreeEvents = out;

    }

    private class TreeEvent{

        private EventType type;
        private double time;
        private AbstractCase aCase;
        private AbstractCase infectorCase;

        private TreeEvent(EventType type, double time, AbstractCase aCase){
            this.type = type;
            this.time = time;
            this.aCase = aCase;
            this.infectorCase = null;
        }

        private TreeEvent(double time, AbstractCase aCase, AbstractCase infectorCase){
            this.type = EventType.INFECTION;
            this.time = time;
            this.aCase = aCase;
            this.infectorCase = infectorCase;
        }

        public double getTime(){
            return time;
        }

        public EventType getType(){
            return type;
        }

        public AbstractCase getCase(){
            return aCase;
        }

        public AbstractCase getInfector(){
            return infectorCase;
        }

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String TRANSMISSION_RATE = "transmissionRate";
        public static final String INITIAL_INFECTION_TIME_PRIOR = "initialInfectionTimePrior";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)
                    xo.getChild(CaseToCaseTreeLikelihood.class);
            SpatialKernel kernel = (SpatialKernel) xo.getChild(SpatialKernel.class);
            Parameter transmissionRate = (Parameter) xo.getElementFirstChild(TRANSMISSION_RATE);

            ParametricDistributionModel iitp = null;

            if(xo.hasChildNamed(INITIAL_INFECTION_TIME_PRIOR)){
                iitp = (ParametricDistributionModel)xo.getElementFirstChild(INITIAL_INFECTION_TIME_PRIOR);
            }


            return new CaseToCaseTransmissionLikelihood(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD,
                    (CategoryOutbreak)c2cTL.getOutbreak(), c2cTL, kernel, transmissionRate, iitp);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        public String getParserDescription() {
            return "This element represents a probability distribution for epidemiological parameters of an outbreak" +
                    "given a phylogenetic tree";
        }

        public Class getReturnType() {
            return CaseToCaseTransmissionLikelihood.class;
        }

        public String getParserName() {
            return CASE_TO_CASE_TRANSMISSION_LIKELIHOOD;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(CaseToCaseTreeLikelihood.class, "The tree likelihood"),
                new ElementRule(SpatialKernel.class, "The spatial kernel", 0, 1),
                new ElementRule(TRANSMISSION_RATE, Parameter.class, "The transmission rate"),
                new ElementRule(INITIAL_INFECTION_TIME_PRIOR, ParametricDistributionModel.class, "The prior " +
                        "probability distibution of the first infection", true)
        };

    };

    // Not the most elegant solution, but you want two types of log out of this model, one for numerical parameters
    // (which Tracer can read) and one for the transmission tree (which it cannot). This is set up so that C2CTransL
    // is the numerical log and C2CTreeL the TT one.

    public LogColumn[] getColumns(){
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>();

        columns.add(new LogColumn.Abstract("trans_LL"){
            protected String getFormattedValue() {
                return String.valueOf(transLogProb);
            }
        });

        columns.add(new LogColumn.Abstract("period_LL") {
            protected String getFormattedValue() {
                return String.valueOf(periodsLogProb);
            }
        });

        columns.addAll(Arrays.asList(treeLikelihood.passColumns()));

        for (AbstractPeriodPriorDistribution hyperprior : (outbreak).getInfectiousMap().values()) {
            columns.addAll(Arrays.asList(hyperprior.getColumns()));
        }

        columns.add(new LogColumn.Abstract("FirstInfectionTime") {
            protected String getFormattedValue() {
                if(sortedTreeEvents==null){
                    sortEvents();
                }
                return String.valueOf(treeLikelihood.getInfectionTime(indexCase));
            }
        });

        columns.add(new LogColumn.Abstract("IndexCaseIndex") {
            protected String getFormattedValue() {
                return String.valueOf(treeLikelihood.getOutbreak().getCaseIndex(indexCase));
            }
        });


        return columns.toArray(new LogColumn[columns.size()]);
    }



}
