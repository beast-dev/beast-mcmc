package dr.inference.model;


import dr.evolution.datatype.DataType;

import java.util.List;

public class StateSet {
    final String id;
    final int[] states;
    final DataType dataType;

    public StateSet(String id,
                    List<Integer> states,
                    DataType dataType) {
        this(id, states.stream().mapToInt(i->i).toArray(), dataType);
    }

    public StateSet(String id,
                    int[] states,
                    DataType dataType) {
        this.id = id;
        this.states = states;
        this.dataType = dataType;
    }

    public int size() { return states.length; }

    public int[] states() { return states; }

    // XT added
    public String getId() {
        return id;
    }
}

//public class StateSet {
//    // two fields
//    private String id;
//    private int[] statesIdArray;
//    private String[] statesNameArray;
//    //private String[] statesNameArray;
//    // constructor
////        public StateSet(){   Use the fields directly (i.e., method(){}) without passing them as parameters for instance method; pass the values as parameters for static method; but I will still keep them to make the code in the upcoming code clear
//    private StateSet(String id, int[] statesIdArray, String[] statesNameArray){
//        setId(id);
//        setStatesId(statesIdArray);
//        setStatesName(statesNameArray);
//    }
//    // getter and setter
//    private void setId(String id){
//        this.id = id;
//    }
//    private void setStatesId(int[] statesIdArray) {
//        this.statesIdArray = statesIdArray;
//    }
//    private void setStatesName(String[] statesNameArray) {
//        this.statesNameArray = statesNameArray;
//    }
//    public String getId() {
//        return id;
//    }
//    public int[] getStatesId() {
//        return statesIdArray;
//    }
//    public String[] statesNameArray() {
//        return statesNameArray;
//    }
//
//}
