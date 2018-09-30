package dr.evolution.tree;

import dr.util.Transform;

/**
 * @author Alex Fisher
 * @author Marc A. Suchard
 */
public class TransformedTreeTraitProvider implements TreeTraitProvider {

    public TransformedTreeTraitProvider(TreeTraitProvider treeTraitProvider,
                                        Transform transform) {

        for (TreeTrait trait : treeTraitProvider.getTreeTraits()) {
            if (trait instanceof TreeTrait.D) {
                transformedTreeTraits.addTrait(createD((TreeTrait.D) trait, transform));
            } else if (trait instanceof TreeTrait.DA) {
                transformedTreeTraits.addTrait(createDA((TreeTrait.DA) trait, transform));
            } else {
                throw new java.lang.IllegalArgumentException("Not transformable trait: " + trait.getTraitName());
            }
        }
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return transformedTreeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return transformedTreeTraits.getTreeTrait(key);
    }

    private TreeTrait.D createD(final TreeTrait.D originalTrait, final Transform transform) {

        return new TreeTrait.D() {
            @Override
            public String getTraitName() {
                return transform.getTransformName() + "." + originalTrait.getTraitName();
            }

            @Override
            public Intent getIntent() {
                return originalTrait.getIntent();
            }

            @Override
            public Double getTrait(Tree tree, NodeRef node) {
                return transform.transform(originalTrait.getTrait(tree, node));
            }

            @Override
            public boolean getLoggable() {
                return true;
            }
        };
    }

    private TreeTrait.DA createDA(final TreeTrait.DA originalTrait, final Transform transform) {

        return new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return transform.getTransformName() + "." + originalTrait.getTraitName();
            }

            @Override
            public Intent getIntent() {
                return originalTrait.getIntent();
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                double[] raw = originalTrait.getTrait(tree, node);
                return transform.transform(raw, 0, raw.length);
            }

            @Override
            public boolean getLoggable() {
                return true;
            }
        };
    }

    private Helper transformedTreeTraits = new Helper();
}
