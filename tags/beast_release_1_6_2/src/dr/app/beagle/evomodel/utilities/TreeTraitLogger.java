package dr.app.beagle.evomodel.utilities;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;

import java.util.ArrayList;
import java.util.List;

/**
 * A utility class to format traits on a tree in a LogColumn-based file.  This class takes an array of TraitProvider
 * and its package location is bound to move
 *
 * @author Marc A. Suchard
 */

public class TreeTraitLogger implements Loggable {

    public TreeTraitLogger(Tree tree,
                           TreeTraitProvider[] traitProviders) {
        this.tree = tree;
        for (TreeTraitProvider provider : traitProviders) {
            addTraits(provider.getTreeTraits());
        }
    }

    public TreeTraitLogger(Tree tree,
                           TreeTrait[] traits) {
        this.tree = tree;
        addTraits(traits);
    }

    public void addTraits(TreeTrait[] traits) {
        if (loggableTreeTraits == null) {
            loggableTreeTraits = new ArrayList<TreeTrait>();
        }
        for (TreeTrait trait : traits) {
            if (trait.getIntent() == TreeTrait.Intent.WHOLE_TREE) {
                loggableTreeTraits.add(trait);
            }
        }
    }

    public LogColumn[] getColumns() {

        if (loggableTreeTraits.size() == 0) {
            return null;
        }

        LogColumn[] columns = new LogColumn[loggableTreeTraits.size()];

        for (int i = 0; i < loggableTreeTraits.size(); i++) {
            final TreeTrait trait = loggableTreeTraits.get(i);
            columns[i] = new LogColumn.Abstract(trait.getTraitName()) {
                @Override
                protected String getFormattedValue() {
                    return trait.getTraitString(tree, null);
                }
            };
        }
        return columns;
    }

    private Tree tree;
    private List<TreeTrait> loggableTreeTraits;
}
