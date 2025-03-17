package dr.inferencexml.model;

import dr.evolution.datatype.DataType;
import dr.inference.model.StronglyLumpableGenerator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import dr.inference.model.StateSet;
import java.util.*;

public class StronglyLumpableGeneratorParser extends AbstractXMLObjectParser {

    public static final String STATE_SET_PARAMETER = "stronglyLumpableGenerator";
    public static final String STATE_SET = "stateSet";
//    public static final String STATE = "state";
//    public static final String CODE = "code";

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DataType.class),
            new ElementRule(StateSet.class, 1, Integer.MAX_VALUE),

//                    new ElementRule(Identifiable.class),
//            new ElementRule(STATE_SET, new XMLSyntaxRule[]{
//                    new ElementRule(StateSet.class),
//            }, true),

     //       }, 1, Integer.MAX_VALUE),
////get from GuassianProcessFieldParser and GeneralDataTypeParser
//            new ElementRule(STATE_SET, new XMLSyntaxRule[] {
////                    new ElementRule(Identifiable.class),
//                    new ElementRule(Identifiable.class, 0, Integer.MAX_VALUE),
//                    new ContentRule("<state code=\"X\"/>"),
//            }, 1, Integer.MAX_VALUE),
//            new ElementRule(dr.inference.model.StateSet.class)
    };




    public StateSet[] parseXMLObject(XMLObject xo) throws XMLParseException {

        // get from glmSubstitutionModel parser
        DataType dataType = (DataType) xo.getChild(DataType.class);

        List<StateSet> states = new ArrayList<>();
        for (StateSet currentSet : xo.getAllChildren(StateSet.class)) {
            states.add(currentSet);
        }

        // Validate state sets for duplicate numbers
        validateStateSets(states, dataType);

        StateSet[] stateSetArray= states.toArray(new StateSet[0]);

        String new_xml = new StronglyLumpableGenerator(stateSetArray, dataType).getXmlOutput();


        //return new_xml;
        //return stateSetArray;
        return null;
    }

    /**
     * Validates if there are duplicate numbers across different StateSet objects.
     * Throws an IllegalArgumentException if duplicates are found.
     */
    private void validateStateSets(List<StateSet> stateSets, DataType dataType) throws XMLParseException {
        Map<Integer, String> seenNumbers = new HashMap<>();

        for (StateSet stateSet : stateSets) {
            for (int num : stateSet.states()) {
                if (seenNumbers.containsKey(num)) {
                    String name = dataType.getCode(num);
                    throw new XMLParseException("Conflict detected: " + stateSet.getId() + " and "
                            + seenNumbers.get(num) + " contain the same state " + name);
                }
                seenNumbers.put(num, stateSet.getId());
            }
        }
    }



    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "generator for strongly lumpable xml.";
    }

    public Class getReturnType() {
        //return LogCtmcRatesStronglyLumpabilityParameter.class;
        return null;
    }

    public String getParserName() {
        //  <ctmcRatesStronglyLumpability  id="hsha">; String LUMP_PARAMETER = "ctmcRatesStronglyLumpability";
        return STATE_SET_PARAMETER;
    }

}


// creat a new class StateSet to store info in a certain StateSet
//    public class StateSet {
//        // two fields
//        private String id;
//        private int[] statesIdArray;
//        private String[] statesNameArray;
//        //private String[] statesNameArray;
//        // constructor
////        public StateSet(){   Use the fields directly (i.e., method(){}) without passing them as parameters for instance method; pass the values as parameters for static method; but I will still keep them to make the code in the upcoming code clear
//        private StateSet(String id, int[] statesIdArray, String[] statesNameArray){
//            setId(id);
//            setStatesId(statesIdArray);
//            setStatesName(statesNameArray);
//        }
//        // getter and setter
//        private void setId(String id){
//            this.id = id;
//        }
//        private void setStatesId(int[] statesIdArray) {
//            this.statesIdArray = statesIdArray;
//        }
//        private void setStatesName(String[] statesNameArray) {
//            this.statesNameArray = statesNameArray;
//        }
//        public String getId() {
//            return id;
//        }
//        public int[] getStatesId() {
//            return statesIdArray;
//        }
//        public String[] statesNameArray() {
//            return statesNameArray;
//        }
//
//    }

//    // define new parser function parseLump(xo)
//    // from GuassianProcessFieldParser
//    private StateSet[] parseStateSet(XMLObject xo, DataType dataType) throws XMLParseException {
//        List<StateSet> lumpsList = new ArrayList<>();
//        // Use a Set to keep track of state ids already assigned to lumps.
//        Set<Integer> usedStateIds = new HashSet<>();
//
//        for (XMLObject cxo : xo.getAllChildren(STATE_SET)) {
//            // for each lump, get first field of lump class: id
//            String currentId = cxo.getId();
//
//            // get second lump field of lump class: states_id_array
//            List<Integer> currentStatesIdList = new ArrayList<Integer>();
//            List<String> currentStatesNameList = new ArrayList<String>();
//            for (XMLObject current_state : cxo.getAllChildren(STATE)) {
//                // from GeneralDataTypeParser
//                // read in the code info (state name) of a state block, like <state code="Cafornia"/>
//                String stateName = current_state.getStringAttribute(CODE);
//                // find the integer corresponding to that state name; this info is stored in dataType
//                int stateId = dataType.getState(stateName);
//
//                // Check if this stateId has already been added to a previous lump.
//                if (usedStateIds.contains(stateId)) {
//                    throw new XMLParseException("State id " + stateId + " (state '" + stateName +
//                            "') in lump '" + currentId + "' is already defined in a previous lump.");
//                }
//
//                // add the new state id to the array
//                currentStatesIdList.add(stateId);
//                currentStatesNameList.add(stateName);
//                usedStateIds.add(stateId);
//            }
//
//            // throw an error if a lump does not have any state
//            // Check if the list is empty
//            if (currentStatesIdList.isEmpty()) {
//                throw new IllegalStateException("Lump" + currentId + "does not have any state in it.");
//            }
//
//            // transfer List<Integer> to int[]
//            int[] currentStatesIdArray = convertListToIntArray(currentStatesIdList);
//            String[] currentStatesNameArray = convertListToStringArray(currentStatesNameList);
//
//            // creat a new lump instance based on the two fields we get
//            StateSet current_Lump = new StateSet(currentId, currentStatesIdArray, currentStatesNameArray);
//            // add current lump to lump array list
//            lumpsList.add(current_Lump);
//        }
//
//        // transfer List<Integer to int[]
//        StateSet[] lumpsArray = convertListToLumpArray(lumpsList);
//
//        return lumpsArray;
//    }

// function change List<Integer> to int[]
//    public static int[] convertListToIntArray(List<Integer> integerList) {
//        // Initialize the array with the size of the List
//        int[] intArray = new int[integerList.size()];
//        // Iterate over the List to convert Integer to int
//        int i = 0;
//        for (Integer value: integerList) {
//            intArray[i++] = value; // Autounboxing converts Integer to int
//        }
//        return intArray;
//    }
//    // function change List<String> to String[]
//    public static String[] convertListToStringArray(List<String> stringList) {
//        // Initialize the array with the size of the List
//        String[] stringArray = new String[stringList.size()];
//        // Iterate over the List to convert Integer to int
//        int i = 0;
//        for (String name: stringList) {
//            stringArray[i++] = name;
//        }
//        return stringArray;
//    }
//    public static StateSet[] convertListToLumpArray(List<StateSet> lumpsList) {
//        // Initialize the array with the size of the List
//        StateSet[] lumpsArray = new StateSet[lumpsList.size()];
//        // Iterate over the List to convert
//        int i = 0;
//        for (StateSet lump:lumpsList) {
//            lumpsArray[i++] = lump;
//        }
//        return lumpsArray;
//    }