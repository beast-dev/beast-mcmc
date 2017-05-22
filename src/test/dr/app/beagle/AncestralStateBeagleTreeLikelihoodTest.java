package test.dr.app.beagle;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import test.dr.inference.trace.TraceCorrelationAssert;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.evolution.util.Taxa;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.math.MathUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.Parameter;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodel.siteratemodel.GammaSiteRateModel;


/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class AncestralStateBeagleTreeLikelihoodTest extends TraceCorrelationAssert {

    private FlexibleTree tree;

    public AncestralStateBeagleTreeLikelihoodTest(String name) {
        super(name);
    }

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

        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gammaModel", mu, null, -1, null);
        siteRateModel.setSubstitutionModel(hky);

        BranchModel branchModel = new HomogeneousBranchModel(
                siteRateModel.getSubstitutionModel());

        BranchRateModel branchRateModel = null;


        AncestralStateBeagleTreeLikelihood treeLikelihood = new AncestralStateBeagleTreeLikelihood(
                alignment,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                null,
                false,
                PartialsRescalingScheme.DEFAULT,
                true,
                null,
                hky.getDataType(),
                "stateTag",
                true, // useMap = true
                false
        );

        double logLike = treeLikelihood.getLogLikelihood();

        StringBuffer buffer = new StringBuffer();

//        Tree.Utils.newick(treeModel, treeModel.getRoot(), false, Tree.BranchLengthType.LENGTHS_AS_TIME,
//                null, null, new NodeAttributeProvider[]{treeLikelihood}, null, null, buffer);
        TreeUtils.newick(treeModel, treeModel.getRoot(), false, TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                null, null, new TreeTraitProvider[] { treeLikelihood }, null, buffer);

        System.out.println(buffer);

        System.out.println("t_CA(2) = " + t(false, 2.0));
        System.out.println("t_CC(1) = " + t(true, 1.0));

        double trueValue = 0.25 * t(false, 2.0) * Math.pow(t(true, 1.0), 3.0);

        assertEquals(logLike, Math.log(trueValue), 1e-6);
    }
}