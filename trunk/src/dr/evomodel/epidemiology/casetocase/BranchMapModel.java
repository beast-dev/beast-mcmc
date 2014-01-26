package dr.evomodel.epidemiology.casetocase;

import dr.evolution.alignment.Alignment;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

/**
 * Created by twoseventwo on 27/12/2013.
 */
public class BranchMapModel extends AbstractModel {

    private AbstractCase[] map;
    private AbstractCase[] storedMap;
    public final static String BRANCH_MAP_MODEL = "branchMapModel";

    public BranchMapModel(PartitionedTreeModel tree){
        super(BRANCH_MAP_MODEL);
        map = new AbstractCase[tree.getNodeCount()];
        storedMap = new AbstractCase[tree.getNodeCount()];
    }

    public void set(int index, AbstractCase aCase, boolean silent){

        AbstractCase oldCase = map[index];

        map[index] = aCase;

        if(!silent){
            pushMapChangedEvent(new BranchMapChangedEvent(index, oldCase, aCase));
        }

    }

    public AbstractCase get(int index){
        return map[index];
    }

    public void pushMapChangedEvents(ArrayList<BranchMapChangedEvent> events){
        listenerHelper.fireModelChanged(this, events);
    }

    public void pushMapChangedEvent(BranchMapChangedEvent event){
        ArrayList<BranchMapChangedEvent> out = new ArrayList<BranchMapChangedEvent>();
        out.add(event);

        pushMapChangedEvents(out);
    }

    // WARNING WARNING WARNING This is to be called ONLY at the start of the run

    public AbstractCase[] getArray(){
        return map;
    }

    public void setAll(AbstractCase[] newMap, boolean silent){
        ArrayList<BranchMapChangedEvent> changes = new ArrayList<BranchMapChangedEvent>();

        for(int i=0; i<map.length; i++){
            if(map[i]!=newMap[i]){
                changes.add(new BranchMapChangedEvent(i, map[i], newMap[i]));
            }
        }

        if(!silent){
            pushMapChangedEvents(changes);
        }

        map = newMap;
    }

    public AbstractCase[] getArrayCopy(){
        return Arrays.copyOf(map, map.length);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // nothing to do
    }

    public int size(){
        return map.length;
    }

    protected void storeState() {
        storedMap = Arrays.copyOf(map, map.length);
    }

    protected void restoreState() {
        map = storedMap;
    }

    protected void acceptState() {
        // nothing to do
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do
    }

    public class BranchMapChangedEvent {

        private final int nodeToRecalculate;
        private final AbstractCase oldCase;
        private final AbstractCase newCase;

        public BranchMapChangedEvent(int node, AbstractCase oldCase, AbstractCase newCase){
            this.nodeToRecalculate = node;
            this.oldCase = oldCase;
            this.newCase = newCase;
        }

        public int getNodeToRecalculate(){
            return nodeToRecalculate;
        }

        public AbstractCase getOldCase() {
            return oldCase;
        }

        public AbstractCase getNewCase() {
            return newCase;
        }

    }

}
