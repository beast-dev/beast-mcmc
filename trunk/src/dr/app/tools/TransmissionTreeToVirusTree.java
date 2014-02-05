package dr.app.tools;

import dr.app.util.Arguments;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.tree.*;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;

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

    public static final String INFECTIONS_FILE_NAME = "infectionsfile";
    public static final String SAMPLES_FILE_NAME = "samplesfile";
    public static final String OUTPUT_NAME_ROOT = "outputroot";

    private DemographicFunction demFunct;
    private ArrayList<InfectedUnit> units;
    private HashMap<String, InfectedUnit> idMap;
    private String outputFileRoot;

    public TransmissionTreeToVirusTree(String infectionEventsFileName, String samplingEventsFileName,
                                       DemographicFunction demFunct, String outputFileRoot){
        this.demFunct = demFunct;
        units = new ArrayList<InfectedUnit>();
        idMap = new HashMap<String, InfectedUnit>();
        this.outputFileRoot = outputFileRoot;
        try {
            readInfectionEvents(infectionEventsFileName);
            readSamplingEvents(samplingEventsFileName);
        } catch(IOException e){
            e.printStackTrace();
        }
    }

    private enum EventType{
        INFECTION, SAMPLE
    }


    private void run() throws IOException{
        FlexibleTree detailedTree = makeTree();
        FlexibleTree simpleTree = makeWellBehavedTree(detailedTree);

        NexusExporter exporter1 = new NexusExporter(new PrintStream(outputFileRoot+"_detailed.nex"));
        exporter1.exportTree(detailedTree);

        NexusExporter exporter2 = new NexusExporter(new PrintStream(outputFileRoot+"_simple.nex"));
        exporter2.exportTree(simpleTree);
    }

    private void readInfectionEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        ArrayList<String[]> lines = new ArrayList<String[]>();
        String line = reader.readLine();

        while(line!=null){
            String[] entries = line.split(",");
            lines.add(entries);

            InfectedUnit unit = new InfectedUnit(entries[0]);
            units.add(unit);
            idMap.put(entries[0], unit);

            line = reader.readLine();
        }

        for(String[] repeatLine: lines){

            InfectedUnit infectee = idMap.get(repeatLine[0]);
            if(!repeatLine[1].equals("none")){
                InfectedUnit infector = idMap.get(repeatLine[1]);
                Event infection = new Event(EventType.INFECTION, Double.parseDouble(repeatLine[2]), infector, infectee);

                infector.addInfectionEvent(infection);
                infectee.setInfectionEvent(infection);

                infectee.parent = infector;
            } else {
                Event infection = new Event(EventType.INFECTION, Double.parseDouble(repeatLine[2]), null, infectee);
                infectee.setInfectionEvent(infection);
            }



        }
    }

    private void readSamplingEvents(String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        String line = reader.readLine();

        while(line!=null){
            String[] entries = line.split(",");

            if(!idMap.containsKey(entries[0])){
                throw new RuntimeException("Trying to add a sampling event to unit "+entries[0]+" but this unit" +
                        "not previously defined");
            }

            InfectedUnit unit = idMap.get(entries[0]);

            unit.addSamplingEvent(Double.parseDouble(entries[1]));

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

        double activeTime = unit.childEvents.get(0).time - unit.infectionEvent.time;

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

    private FlexibleTree makeTree(){

        // find the first case

        InfectedUnit firstCase = null;

        for(InfectedUnit unit : units){
            if(unit.parent==null){
                firstCase = unit;
            }
        }

        if(firstCase==null){
            throw new RuntimeException("Can't find a first case");
        }

        FlexibleNode outTreeRoot = makeSubtree(firstCase);

        FlexibleTree finalTree = new FlexibleTree(outTreeRoot, false, true);

        return finalTree;

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


        CoalescentSimulator simulator = new CoalescentSimulator();
        SimpleNode root;

        SimpleNode[] simResults = simulator.simulateCoalescent(nodes.toArray(new SimpleNode[nodes.size()]),
                demogFunct, -maxHeight, 0, true);
        if(simResults.length==1){
            root = simResults[0];
        } else {
            throw new RuntimeException("Multiple lineages still exist at time of infection (probably a rounding " +
                    "problem)");
        }

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
                throw new RuntimeException("Adding an event to case "+id+" before its infection time");
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

    public static void main(String[] args){

        ExponentialGrowth testGrowth = new ExponentialGrowth(Units.Type.YEARS);
        testGrowth.setN0(1);
        testGrowth.setGrowthRate(1);

        TransmissionTreeToVirusTree instance = new TransmissionTreeToVirusTree("inf_test.csv", "samp_test.csv",
                testGrowth, "test_out");

        try{
            instance.run();
        } catch (IOException e){
            e.printStackTrace();
        }
    }
}
