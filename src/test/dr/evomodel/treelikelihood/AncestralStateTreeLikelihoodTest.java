package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.TestCase;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class AncestralStateTreeLikelihoodTest extends TestCase {

    private FlexibleTree tree;

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        NewickImporter importer = new NewickImporter("(0:2.0,(1:1.0,2:1.0):1.0);");
        tree = (FlexibleTree) importer.importTree(null);
    }

    // transition prob of JC69
    private double t(boolean same, double time) {
        if (same) {
            return 0.25 + 0.75 * Math.exp(-4.0 / 3.0 * time);
        } else {
            return 0.25 - 0.25 * Math.exp(-4.0 / 3.0 * time);
        }
    }

    public void testJointLikelihood() {

        TreeModel treeModel = new TreeModel("treeModel", tree);

        Sequence[] sequence = new Sequence[3];

        sequence[0] = new Sequence(new Taxon("0"), "A");
        sequence[1] = new Sequence(new Taxon("1"), "C");
        sequence[2] = new Sequence(new Taxon("2"), "C");

        Taxa taxa = new Taxa();
        for (Sequence s : sequence) {
            taxa.addTaxon(s.getTaxon());
        }

        SimpleAlignment alignment = new SimpleAlignment();
        for (Sequence s : sequence) {
            alignment.addSequence(s);
        }

        Parameter mu = new Parameter.Default(1, 1.0);

        Parameter kappa = new Parameter.Default(1, 1.0);
        double[] pi = {0.25, 0.25, 0.25, 0.25};

        Parameter freqs = new Parameter.Default(pi);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        AncestralStateTreeLikelihood treeLikelihood = new AncestralStateTreeLikelihood(
                alignment, treeModel,
                new GammaSiteModel(hky), new StrictClockBranchRates(mu),
                false, true,
                Nucleotides.INSTANCE,
                "state",
                false,
                true, // useMap = true
                false);

        double logLike = treeLikelihood.getLogLikelihood();

        StringBuffer buffer = new StringBuffer();

        TreeUtils.newick(treeModel, treeModel.getRoot(), false, TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                null, null, new TreeTraitProvider[] { treeLikelihood }, null, buffer); 


        System.out.println(buffer);

        System.out.println("t_CA(2) = " + t(false, 2.0));
        System.out.println("t_CC(1) = " + t(true, 1.0));

        double trueValue = 0.25 * t(false, 2.0) * Math.pow(t(true, 1.0), 3.0);

        assertEquals(logLike, Math.log(trueValue), 1e-6);
    }
}
