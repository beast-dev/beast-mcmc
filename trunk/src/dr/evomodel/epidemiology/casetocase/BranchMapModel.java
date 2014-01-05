package dr.evomodel.epidemiology.casetocase;

import dr.evolution.alignment.Alignment;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * Created by twoseventwo on 27/12/2013.
 */
public class BranchMapModel extends AbstractModel {

    private AbstractCase[] map;
    private AbstractCase[] storedMap;
    public final static String BRANCH_MAP_MODEL = "branchMapModel";

    public BranchMapModel(TreeModel tree){
        super(BRANCH_MAP_MODEL);
        map = new AbstractCase[tree.getNodeCount()];
        storedMap = new AbstractCase[tree.getNodeCount()];
    }

    public void set(int index, AbstractCase aCase){
        AbstractCase oldCase = map[index];
        map[index] = aCase;
        if(oldCase!=aCase){
            fireModelChanged();
        }
    }

    public AbstractCase get(int index){
        return map[index];
    }

    public void setAll(AbstractCase[] newMap){
        map = newMap;
        fireModelChanged();
    }

    public AbstractCase[] getArrayCopy(){
        return Arrays.copyOf(map, map.length);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // nothing to do
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
}
