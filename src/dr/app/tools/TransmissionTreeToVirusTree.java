/*
 * TransmissionTreeToVirusTree.java
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

package dr.app.tools;

import dr.app.util.Arguments;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.tree.*;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.epidemiology.LogisticGrowthN0;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Simulated a virus tree given a transmission tree and dates of sampling
 *
 * @author mhall
 */

public class TransmissionTreeToVirusTree {

    protected static PrintStream progressStream = System.out;

    private enum ModelType{CONSTANT, EXPONENTIAL, LOGISTIC}

    public static final String HELP = "help";

    public static final String DEMOGRAPHIC_MODEL = "demoModel";

    public static final String[] demographics = {"Constant", "Exponential", "Logistic"};

    public static final String STARTING_POPULATION_SIZE = "N0";
    public static final String GROWTH_RATE = "growthRate";
    public static final String T50 = "t50";

    private DemographicFunction demFunct;
    private ArrayList<InfectedUnit> units;
    private HashMap<String, InfectedUnit> idMap;
    private String outputFileRoot;

    private double coalescentProbability;

    public TransmissionTreeToVirusTree(String fileName,
                                       DemographicFunction demFunct, String outputFileRoot){
        this.demFunct = demFunct;
        units = new ArrayList<InfectedUnit>();
        idMap = new HashMap<String, InfectedUnit>();
        this.outputFileRoot = outputFileRoot;
        coalescentProbability = 1;
        try {
            readInfectionEvents(fileName);
            readSamplingEvents(fileName);
        } catch(IOException e){
            e.printStackTrace();
        }
    }

    public TransmissionTreeToVirusTree(String sampFileName, String transFileName,
                                       DemographicFunction demFunct, String outputFileRoot){
        this.demFunct = demFunct;
        units = new ArrayList<InfectedUnit>();
        idMap = new HashMap<String, InfectedUnit>();
        this.outputFileRoot = outputFileRoot;
        try {
            readInfectionEvents(transFileName);
            readSamplingEvents(sampFileName);
        } catch(IOException e){
            e.printStackTrace();
        }
    }



    private enum EventType{
        INFECTION, SAMPLE
    }


    private void run() throws IOException{
        ArrayList<FlexibleTree> detailedTrees = makeTrees();
        ArrayList<FlexibleTree> simpleTrees = new ArrayList<FlexibleTree>();

        for(FlexibleTree tree : detailedTrees) {
            FlexibleTree wbTree = makeWellBehavedTree(tree);
            wbTree.setAttribute("firstCase", tree.getAttribute("firstCase"));

            simpleTrees.add(wbTree);
        }


        for(FlexibleTree tree: detailedTrees){
            NexusExporter exporter = new NexusExporter(new PrintStream(outputFileRoot
                    + tree.getAttribute("firstCase") + "_detailed.nex"));
            exporter.exportTree(tree);
        }

        for(FlexibleTree tree: simpleTrees){
            NexusExporter exporter = new NexusExporter(new PrintStream(outputFileRoot
                    + tree.getAttribute("firstCase") + "_simple.nex"));
            exporter.exportTree(tree);
        }

    }

    private void readInfectionEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        ArrayList<String[]> keptLines = new ArrayList<String[]>();

        reader.readLine();

        String line = reader.readLine();

        while(line!=null){

            String[] entries = line.split(",");
            keptLines.add(entries);

            InfectedUnit unit = new InfectedUnit("ID_"+entries[1]);
            units.add(unit);
            idMap.put("ID_"+entries[1], unit);

            line = reader.readLine();
        }

        for(String[] repeatLine: keptLines){

            InfectedUnit infectee = idMap.get("ID_"+repeatLine[1]);
            if(!repeatLine[2].equals("-1")){
                InfectedUnit infector = idMap.get("ID_"+repeatLine[2]);
                Event infection = new Event(EventType.INFECTION, Double.parseDouble(repeatLine[3]), infector, infectee);

                infector.addInfectionEvent(infection);
                infectee.setInfectionEvent(infection);

                infectee.parent = infector;
            } else {
                Event infection = new Event(EventType.INFECTION, Double.parseDouble(repeatLine[3]), null, infectee);
                infectee.setInfectionEvent(infection);
            }
        }
    }

    private void readSamplingEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        reader.readLine();

        String line = reader.readLine();

        while(line!=null){

            String[] entries = line.split(",");

            if(!entries[7].equals("NA")) {

                if (!idMap.containsKey("ID_"+entries[1])) {
                    throw new RuntimeException("Trying to add a sampling event to unit " + entries[2] + " but this " +
                            "unit not previously defined");
                }

                InfectedUnit unit = idMap.get("ID_"+entries[1]);

                unit.addSamplingEvent(Double.parseDouble(entries[7]));
            }

            line = reader.readLine();
        }
    }

    // events are only relevant if there is a sampling event somewhere further up the tree

    private FlexibleTree makeTreelet(InfectedUnit unit, ArrayList<Event> relevantEvents){

        if(relevantEvents.size()==0){
            return null;
        }

        ArrayList<SimpleNode> nodes = new ArrayList<SimpleNode>();

        unit.sortEvents();

        double lastRelevantEventTime = Double.NEGATIVE_INFINITY;

        for(Event event : relevantEvents){
            if(event.time > lastRelevantEventTime){
                lastRelevantEventTime = event.time;
            }
        }

        double activeTime = lastRelevantEventTime - unit.infectionEvent.time;

        for(Event event : relevantEvents){
            Taxon taxon;
            if(event.type == EventType.INFECTION){
                taxon = new Taxon(event.infectee.id+"_infected_by_"+event.infector.id+"_"+event.time);
            } else {
                taxon = new Taxon(unit.id+"_sampled_"+event.time);
            }
            taxon.setDate(new Date(event.time - unit.infectionEvent.time, Units.Type.YEARS, false));
            SimpleNode node = new SimpleNode();
            node.setTaxon(taxon);
            nodes.add(node);
            node.setHeight(unit.infectionEvent.time - event.time);
            node.setAttribute("Event", event);
        }

        FlexibleNode treeletRoot;

        if(nodes.size()>1){
            treeletRoot = simulateCoalescent(nodes, demFunct, activeTime);
        } else {
            treeletRoot = new FlexibleNode(new SimpleTree(nodes.get(0)), nodes.get(0), true);
            treeletRoot.setHeight(0);
        }

        // add the root branch length

        FlexibleNode infectionNode = new FlexibleNode();
        infectionNode.setHeight(activeTime);
        infectionNode.addChild(treeletRoot);
        treeletRoot.setLength(activeTime - treeletRoot.getHeight());
        infectionNode.setAttribute("Event", unit.infectionEvent);

        FlexibleTree outTree = new FlexibleTree(infectionNode);

        for(int i=0; i<outTree.getNodeCount(); i++){
            FlexibleNode node = (FlexibleNode)outTree.getNode(i);
            node.setAttribute("Unit", unit.id);
        }

        return outTree;
    }

    private ArrayList<FlexibleTree> makeTrees(){

        // find the first case

        ArrayList<InfectedUnit> introducedCases = new ArrayList<InfectedUnit>();

        for(InfectedUnit unit : units){
            if(unit.parent==null){
                introducedCases.add(unit);
            }
        }

        if(introducedCases.size()==0){
            throw new RuntimeException("Can't find a first case");
        }

        ArrayList<FlexibleTree> out = new ArrayList<FlexibleTree>();

        for(InfectedUnit introduction : introducedCases) {
            coalescentProbability = 1;

            System.out.println("Building tree for descendants of " + introduction.id);
            FlexibleNode outTreeRoot = makeSubtree(introduction);

            if (outTreeRoot != null) {
                FlexibleTree finalTree = new FlexibleTree(outTreeRoot, false, true);
                finalTree.setAttribute("firstCase", introduction.id);
                out.add(finalTree);

                if(coalescentProbability<0.9){
                    progressStream.println("WARNING: any phylogeny for descendants of "+introduction.id+" is quite " +
                            "improbable (p<"+(coalescentProbability)+") given this demographic function. Consider " +
                            "another.");
                }
            } else {
                progressStream.println("This individual has no sampled descendants");
            }
            System.out.println();


        }
        return out;
    }

    // make the tree from this unit up

    private FlexibleNode makeSubtree(InfectedUnit unit){

        HashMap<Event, FlexibleNode> eventToSubtreeRoot = new HashMap<Event, FlexibleNode>();

        ArrayList<Event> relevantEvents = new ArrayList<Event>();

        for(Event event : unit.childEvents){

            if(event.type == EventType.INFECTION){

                FlexibleNode childSubtreeRoot = makeSubtree(event.infectee);

                if(childSubtreeRoot!=null){
                    relevantEvents.add(event);
                    eventToSubtreeRoot.put(event, childSubtreeRoot);
                }

            } else if(event.type==EventType.SAMPLE) {
                relevantEvents.add(event);
            }
        }

        FlexibleTree unitTreelet = makeTreelet(unit, relevantEvents);

        if(unitTreelet==null){
            return null;
        }

        for(int i=0; i<unitTreelet.getExternalNodeCount(); i++){

            FlexibleNode tip = (FlexibleNode)unitTreelet.getExternalNode(i);

            Event tipEvent = (Event)unitTreelet.getNodeAttribute(tip, "Event");

            if(tipEvent.type == EventType.INFECTION){
                FlexibleNode subtreeRoot = eventToSubtreeRoot.get(tipEvent);

                FlexibleNode firstSubtreeSplit = subtreeRoot.getChild(0);

                subtreeRoot.removeChild(firstSubtreeSplit);
                tip.addChild(firstSubtreeSplit);

            }
        }

        return (FlexibleNode)unitTreelet.getRoot();
    }

    private FlexibleNode simulateCoalescent(ArrayList<SimpleNode> nodes, DemographicFunction demogFunct,
                                            double maxHeight){

        double earliestNodeHeight = Double.NEGATIVE_INFINITY;

        for(SimpleNode node : nodes){
            if(node.getHeight()>earliestNodeHeight){
                earliestNodeHeight = node.getHeight();
            }
        }
        double maxLastInterval = earliestNodeHeight;
        double probNoCoalesenceInTime = Math.exp(demogFunct.getIntensity(maxLastInterval));

        coalescentProbability *= (1-probNoCoalesenceInTime);



        CoalescentSimulator simulator = new CoalescentSimulator();
        SimpleNode root;

        SimpleNode[] simResults;
        int failCount = 0;

        do {
            simResults = simulator.simulateCoalescent(nodes.toArray(new SimpleNode[nodes.size()]),
                    demogFunct, -maxHeight, 0, true);
            if(simResults.length>1){
                failCount++;
                System.out.println("Failed to coalesce lineages: "+failCount);
            }
        } while(simResults.length!=1);

        root = simResults[0];

        SimpleTree simpleTreelet = new SimpleTree(root);


        for (int i=0; i<simpleTreelet.getNodeCount(); i++) {
            SimpleNode node = (SimpleNode)simpleTreelet.getNode(i);
            node.setHeight(node.getHeight() + maxHeight);
        }

        return new FlexibleNode(simpleTreelet, root, true);
    }

    private FlexibleTree makeWellBehavedTree(FlexibleTree tree){
        FlexibleTree newPhylogeneticTree = new FlexibleTree(tree, false);

        newPhylogeneticTree.beginTreeEdit();
        for(int i=0; i<newPhylogeneticTree.getInternalNodeCount(); i++){
            FlexibleNode node = (FlexibleNode)newPhylogeneticTree.getInternalNode(i);
            if(newPhylogeneticTree.getChildCount(node)==1){
                FlexibleNode parent = (FlexibleNode)newPhylogeneticTree.getParent(node);
                FlexibleNode child = (FlexibleNode)newPhylogeneticTree.getChild(node, 0);
                if(parent!=null){
                    double childHeight = newPhylogeneticTree.getNodeHeight(child);
                    newPhylogeneticTree.removeChild(parent, node);
                    newPhylogeneticTree.addChild(parent, child);
                    newPhylogeneticTree.setNodeHeight(child, childHeight);
                } else {
                    child.setParent(null);
                    newPhylogeneticTree.setRoot(child);
                }
            }
        }
        newPhylogeneticTree.endTreeEdit();


        return new FlexibleTree(newPhylogeneticTree, true);
    }

    private class InfectedUnit{
        private String id;
        private ArrayList<Event> childEvents;
        private Event infectionEvent;
        private InfectedUnit parent;

        private InfectedUnit(String id){
            this.id = id;
            parent = null;
            childEvents = new ArrayList<Event>();
        }

        private void addSamplingEvent(double time){
            if(time < infectionEvent.time){
                throw new RuntimeException("Adding an event to case "+id+" before its infection time");
            }
            childEvents.add(new Event(EventType.SAMPLE, time));
        }

        private void setInfectionEvent(double time, InfectedUnit infector){
            setInfectionEvent(new Event(EventType.INFECTION, time, infector, this));
        }

        private void setInfectionEvent(Event event){
            for(Event childEvent : childEvents){
                if(event.time > childEvent.time){
                    throw new RuntimeException("Setting infection time for case "+id+" after an existing child event");
                }
            }

            infectionEvent = event;
        }

        private void addChildInfectionEvent(double time, InfectedUnit infectee){
            addInfectionEvent(new Event(EventType.INFECTION, time, this, infectee));
        }

        private void addInfectionEvent(Event event){
            if(infectionEvent!=null && event.time < infectionEvent.time){
                throw new RuntimeException("Adding an infection event to case "+id+" at "+event.time+" before its " +
                        "infection time at "+infectionEvent.time);
            }
            childEvents.add(event);
        }

        private void sortEvents(){
            Collections.sort(childEvents);
            Collections.reverse(childEvents);
        }

    }

    private class Event implements Comparable<Event>{

        private EventType type;
        private double time;
        private InfectedUnit infector;
        private InfectedUnit infectee;


        private Event(EventType type, double time){
            this.type = type;
            this.time = time;
        }

        private Event(EventType type, double time, InfectedUnit infector, InfectedUnit infectee){
            this.type = type;
            this.time = time;
            this.infector = infector;
            this.infectee = infectee;
        }

        public int compareTo(Event event) {
            return Double.compare(time, event.time);
        }
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("virusTreeBuilder", "<infections-file-name> <sample-file-name> <output-file-name-root>");
    }



    public static void main(String[] args){

        ModelType model = ModelType.CONSTANT;
        double startNe = 1;
        double growthRate = 0;
        double t50 = 0;

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption(DEMOGRAPHIC_MODEL, demographics, false, "The type of within-host" +
                                " demographic function to use, default = constant"),
                        new Arguments.RealOption(STARTING_POPULATION_SIZE,"The effective population size at time zero" +
                                " (used in all models), default = 1"),
                        new Arguments.RealOption(GROWTH_RATE,"The effective population size growth rate (used in" +
                                " exponential and logistic models), default = 0"),
                        new Arguments.RealOption(T50,"The time point, relative to the time of infection in backwards " +
                                "time, at which the population is equal to half its final asymptotic value, in the " +
                                "logistic model default = 0")
                });


        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption(HELP)) {
            printUsage(arguments);
            System.exit(0);
        }

        if (arguments.hasOption(DEMOGRAPHIC_MODEL)) {
            String modelString = arguments.getStringOption(DEMOGRAPHIC_MODEL);
            if(modelString.toLowerCase().startsWith("c")){
                model = ModelType.CONSTANT;
            } else if(modelString.toLowerCase().startsWith("e")){
                model = ModelType.EXPONENTIAL;
            } else if(modelString.toLowerCase().startsWith("l")){
                model = ModelType.LOGISTIC;
            } else {
                progressStream.print("Unrecognised demographic model type");
                System.exit(1);
            }
        }


        if(arguments.hasOption(STARTING_POPULATION_SIZE)){
            startNe = arguments.getRealOption(STARTING_POPULATION_SIZE);
        }

        if(arguments.hasOption(GROWTH_RATE) && model!=ModelType.CONSTANT){
            growthRate = arguments.getRealOption(GROWTH_RATE);
        }

        if(arguments.hasOption(T50) && model==ModelType.LOGISTIC){
            t50 = arguments.getRealOption(T50);
        }

        DemographicFunction demoFunction = null;

        switch(model){
            case CONSTANT: {
                demoFunction = new ConstantPopulation(Units.Type.YEARS);
                ((ConstantPopulation)demoFunction).setN0(startNe);
            }
            case EXPONENTIAL: {
                demoFunction = new ExponentialGrowth(Units.Type.YEARS);
                ((ExponentialGrowth)demoFunction).setN0(startNe);
                ((ExponentialGrowth)demoFunction).setGrowthRate(growthRate);
            }
            case LOGISTIC: {
                demoFunction = new LogisticGrowthN0(Units.Type.YEARS);
                ((LogisticGrowthN0)demoFunction).setN0(startNe);
                ((LogisticGrowthN0)demoFunction).setGrowthRate(growthRate);
                ((LogisticGrowthN0)demoFunction).setT50(t50);
            }
        }

        final String[] args2 = arguments.getLeftoverArguments();

        if(args2.length!=3){
            printUsage(arguments);
            System.exit(1);
        }

        String infectionsFileName = args2[0];
        String samplesFileName = args2[1];
        String outputFileRoot = args2[2];


        TransmissionTreeToVirusTree instance = new TransmissionTreeToVirusTree(samplesFileName,
                infectionsFileName, demoFunction, outputFileRoot);

        try{
            instance.run();
        } catch (IOException e){
            e.printStackTrace();
        }
    }
}
