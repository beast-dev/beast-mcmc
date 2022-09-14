package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.HashMap;

/**
 * @author Gabriel Hassler
 * @author Marc Suchard
 */


public class ConditionalTraitSimulationHelper implements Reportable {

    private final TreeDataLikelihood treeLikelihood;
    private final TreeTrait treeTrait;
    private final ContinuousTraitPartialsProvider topDataModel;
    private final HashMap<ContinuousTraitPartialsProvider, ParentMapHelper> parentMap;


    public ConditionalTraitSimulationHelper(TreeDataLikelihood treeLikelihood) {
        this.treeLikelihood = treeLikelihood;
        ContinuousDataLikelihoodDelegate delegate = (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();
        this.topDataModel = delegate.getDataModel();

        this.parentMap = new HashMap<>();
        makeParentMap(topDataModel, parentMap);

        this.treeTrait = treeLikelihood.getTreeTrait(topDataModel.getTipTraitName());
    }

    public TreeTrait getTreeTrait() {
        return treeTrait;
    }

    private class ParentMapHelper {
        public final int traitOffset;
        public final int traitDimension;
        public final ContinuousTraitPartialsProvider parent;

        private ParentMapHelper(ContinuousTraitPartialsProvider parent, int traitOffset, int traitDimension) {
            this.traitOffset = traitOffset;
            this.traitDimension = traitDimension;
            this.parent = parent;
        }
    }

    private void makeParentMap(ContinuousTraitPartialsProvider model,
                               HashMap<ContinuousTraitPartialsProvider, ParentMapHelper> map) {

        int offset = 0;
        for (ContinuousTraitPartialsProvider child : model.getChildModels()) {

            map.put(child, new ParentMapHelper(model, offset, child.getTraitDimension()));
            makeParentMap(child, map);
            offset += child.getTraitDimension();
        }
    }

    public double[] drawTraitsAbove(ContinuousTraitPartialsProvider model) {
        int dimTrait = model.getTraitDimension();

        if (model == topDataModel) {
            treeLikelihood.fireModelChanged();
            return (double[]) treeTrait.getTrait(treeLikelihood.getTree(), null);
        }

        ParentMapHelper helper = parentMap.get(model);

        double[] fullTraitsAbove = drawTraitsBelow(helper.parent);
        if (helper.traitOffset == 0 && helper.traitDimension == helper.parent.getDataDimension()) {
            return fullTraitsAbove;
        }

        int nTaxa = treeLikelihood.getTree().getExternalNodeCount();
        double[] traitsAbove = new double[nTaxa * dimTrait];

        int fullOffset = helper.traitOffset;
        int thisOffset = 0;
        int dimAbove = helper.parent.getTraitDimension();
        for (int i = 0; i < nTaxa; i++) {
            System.arraycopy(fullTraitsAbove, fullOffset, traitsAbove, thisOffset, helper.traitDimension);
            fullOffset += dimAbove;
            thisOffset += dimTrait;
        }

        return traitsAbove;
    }

    public double[] drawTraitsBelow(ContinuousTraitPartialsProvider model) {
        double[] aboveTraits = drawTraitsAbove(model);
        return model.drawTraitsBelowConditionalOnDataAndTraitsAbove(aboveTraits);
    }

    public class JointSamples {

        private final double[] traitsAbove;
        private final double[] traitsBelow;

        public JointSamples(double[] traitsAbove, double[] traitsBelow) {
            this.traitsAbove = traitsAbove;
            this.traitsBelow = traitsBelow;
        }

        public double[] getTraitsAbove() {
            return traitsAbove;
        }

        public double[] getTraitsBelow() {
            return traitsBelow;
        }
    }

    public JointSamples drawTraitsAboveAndBelow(ContinuousTraitPartialsProvider model) {
        return drawTraitsAboveAndBelow(model, false);
    }

    public JointSamples drawTraitsAboveAndBelow(ContinuousTraitPartialsProvider model, boolean transformAbove) {
        double[] aboveTraits = drawTraitsAbove(model);
        double[] belowTraits = model.drawTraitsBelowConditionalOnDataAndTraitsAbove(aboveTraits);
        if (transformAbove) {
            aboveTraits = model.transformTreeTraits(aboveTraits); //TODO: this is probably done twice for something like a latent factor model. can be more efficient
        }
        return new JointSamples(aboveTraits, belowTraits);
    }


    @Override
    public String getReport() {
        int repeats = 10000;

        double[] mean = drawTraitsAbove(topDataModel);
        for (int i = 1; i < repeats; i++) {
            double[] draw = drawTraitsAbove(topDataModel);
            for (int j = 0; j < draw.length; j++) {
                mean[j] += draw[j];
            }
        }

        for (int i = 0; i < mean.length; i++) {
            mean[i] /= repeats;
        }

        StringBuilder sb = new StringBuilder("Trait simulation report:\n\ttree trait mean: ");
        sb.append(new Vector(mean));
        sb.append("\n");

        return sb.toString();
    }


}
