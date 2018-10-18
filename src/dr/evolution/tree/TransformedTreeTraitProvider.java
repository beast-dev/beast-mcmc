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
            if (trait.getTraitClass() == Double.class) {
                transformedTreeTraits.addTrait(createD(trait, transform));
            } else if (trait.getTraitClass() == double[].class) {
                transformedTreeTraits.addTrait(createDA(trait, transform));
            } else {
                throw new IllegalArgumentException("Not transformable trait: " + trait.getTraitName());
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

    private TreeTrait.D createD(final TreeTrait originalTrait, final Transform transform) {

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
                return transform.transform((Double) originalTrait.getTrait(tree, node));
            }

            @Override
            public boolean getLoggable() {
                return true;
            }
        };
    }

    private TreeTrait.DA createDA(final TreeTrait originalTrait, final Transform transform) {

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
                double[] raw = (double[]) originalTrait.getTrait(tree, node);
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
