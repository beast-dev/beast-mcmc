package dr.evomodel.epidemiology.casetocase;

import dr.evolution.alignment.Alignment;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by twoseventwo on 27/12/2013.
 */
public class BranchMapModel extends AbstractModel {

    private TreeModel tree;
    private AbstractCase[] map;
    private AbstractCase[] storedMap;
    public final static String BRANCH_MAP_MODEL = "branchMapModel";

    public BranchMapModel(TreeModel tree){
        super(BRANCH_MAP_MODEL);
        this.tree = tree;
        map = new AbstractCase[tree.getNodeCount()];
        storedMap = new AbstractCase[tree.getNodeCount()];
    }

    public void set(int index, AbstractCase aCase){
        AbstractCase oldCase = map[index];
        ArrayList<AbstractCase> changes = new ArrayList<AbstractCase>();
        changes.add(oldCase);
        changes.add(aCase);

        map[index] = aCase;
        if(oldCase!=aCase){
            pushMapChangedEvent(changes);
        }
    }

    public AbstractCase get(int index){
        return map[index];
    }

    public void pushMapChangedEvent(ArrayList<AbstractCase> casesToRecalculate){
        listenerHelper.fireModelChanged(this, new BranchMapChangedEvent(casesToRecalculate));
    }

    // WARNING WARNING WARNING This is to be called ONLY at the start of the run

    public void setupMap(AbstractCase[] map){
        this.map = map;
    }

    public AbstractCase[] getArray(){
        return map;
    }

    public void setAll(AbstractCase[] newMap){
        ArrayList<AbstractCase> changes = new ArrayList<AbstractCase>();

        for(int i=0; i<map.length; i++){
            if(map[i]!=newMap[i]){
                if(map[i]!=null){
                    changes.add(map[i]);
                }
                changes.add(newMap[i]);
            }
        }

        pushMapChangedEvent(changes);

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

        private final ArrayList<AbstractCase> casesToRecalculate;

        public BranchMapChangedEvent(ArrayList<AbstractCase> casesToRecalculate){
            this.casesToRecalculate = casesToRecalculate;
        }

        public ArrayList<AbstractCase> getCasesToRecalculate(){
            return casesToRecalculate;
        }

    }

}
