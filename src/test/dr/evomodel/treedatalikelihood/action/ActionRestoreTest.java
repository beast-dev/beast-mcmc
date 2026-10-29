package test.dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

/**
 * The action-based likelihood overwrites the rate matrix in BEAGLE in place, so a rejected substitution-model move
 * must put the restored matrix back.  This drives the store / propose / reject / restore sequence of an MCMC step
 * by hand, then checks the next evaluation against a likelihood built from scratch with the same parameter values.
 * (MarkovChain's own full-evaluation check cannot see this: it compares with makeDirty( ), which reads the same cached
 * matrix entries and then re-sends them.)  Rescaling is "always", as in the action XMLs.
 */
public class ActionRestoreTest extends TraceCorrelationAssert {

    private TreeModel treeModel;

    public ActionRestoreTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();
    }

    // Whether the BEAGLE library supports action-based likelihoods (release builds without them have no
    // setSparseMatrix( ), so the call throws UnsatisfiedLinkError).
    private static boolean actionAvailable() {
        try {
            Beagle beagle = BeagleFactory.loadBeagleInstance(2, 3, 0, 4, 1, 1, 2, 1, 0, null,
                    BeagleFlag.PRECISION_DOUBLE.getMask() | BeagleFlag.BEAGLE_FLAG_COMPUTATION_ACTION.getMask(), 0);
            if (beagle == null) {
                return false;
            }
            try {
                beagle.setSparseMatrix(0, new int[]{0, 1, 2, 3}, new int[]{0, 1, 2, 3}, new double[]{0, 0, 0, 0}, 4);
            } finally {
                beagle.finalize();
            }
            return true;
        } catch (Throwable e) {
            return false;
        }
    }

    private TreeDataLikelihood createLikelihood(Parameter kappa, Parameter clockRate) {
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        HKY hky = new HKY(kappa, frequencies);
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma", null, 1.0, new Parameter.Default(0.5), 4, null);
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(
                treeModel, patterns, new HomogeneousBranchModel(hky), siteRateModel,
                true, false, PartialsRescalingScheme.ALWAYS, false,
                new PreOrderSettings(false, false, false, true, true, false));
        return new TreeDataLikelihood(delegate, treeModel, new StrictClockBranchRates(clockRate));
    }

    public void testRejectedSubstitutionMoveIsRestored() {
        if (!actionAvailable()) {
            System.out.println("ActionRestoreTest skipped: the BEAGLE library does not support action-based likelihoods");
            return;
        }

        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 2.0);
        Parameter clockRate = new Parameter.Default(1.0);
        TreeDataLikelihood likelihood = createLikelihood(kappa, clockRate);

        final double before = likelihood.getLogLikelihood();

        // A kappa move that is rejected.
        likelihood.storeModelState();
        kappa.setParameterValue(0, 6.0);
        final double proposed = likelihood.getLogLikelihood();
        assertTrue("the kappa move should change the likelihood", Math.abs(proposed - before) > 1e-3);
        likelihood.restoreModelState();
        assertEquals("restored likelihood", before, likelihood.getLogLikelihood(), 0.0);

        // The next move changes only the clock rate, so the substitution model is not re-sent.
        likelihood.storeModelState();
        clockRate.setParameterValue(0, 1.3);
        final double incremental = likelihood.getLogLikelihood();
        likelihood.makeDirty();
        final double full = likelihood.getLogLikelihood();

        // The same parameter values, from scratch.
        final double expected = createLikelihood(new Parameter.Default(HKYParser.KAPPA, 2.0), new Parameter.Default(1.3)).getLogLikelihood();

        System.out.println("action restore: before " + before + "  proposed " + proposed + "  next move: incremental "
                + incremental + "  makeDirty " + full + "  from scratch " + expected);
        assertEquals("incremental likelihood after a rejected kappa move", expected, incremental, 1e-8 * Math.abs(expected));
        assertEquals("full likelihood after a rejected kappa move", expected, full, 1e-8 * Math.abs(expected));
    }
}
