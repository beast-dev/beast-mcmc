package test.dr.evomodel.treedatalikelihood.action;

import junit.framework.TestCase;

/**
 * The action-likelihood CPU/GPU dispatch rule (BeagleDataLikelihoodDelegate.actionPrefersCpu) against the
 * measured crossovers it encodes (A100 vs 16 CPU threads, November 2026).
 */
public class ActionDispatchTest extends TestCase {

    public ActionDispatchTest(String name) {
        super(name);
    }

    private static boolean cpu(int taxa, int states, int patterns, int categories) {
        return dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate.actionPrefersCpu(states, patterns, categories, taxa);
    }

    public void testTinyProblemsUseTheCpu() {
        assertTrue(cpu(10, 32, 1, 1));    // CPU 3.5x faster
        assertTrue(cpu(10, 61, 1, 1));    // 3.3x
        assertTrue(cpu(10, 32, 16, 1));   // 1.5x
        assertTrue(cpu(10, 128, 4, 1));   // 1.1x
        assertTrue(cpu(50, 32, 1, 1));    // 1.6x
    }

    public void testEverythingElseUsesTheGpu() {
        assertFalse(cpu(10, 128, 1, 4));  // GPU 2.0x faster
        assertFalse(cpu(10, 61, 4, 4));   // 1.5x
        assertFalse(cpu(50, 128, 1, 1));  // 2.4x
        assertFalse(cpu(10, 256, 1, 1));  // 1.3x
        assertFalse(cpu(10, 488, 1, 1));  // 2.2x
        assertFalse(cpu(50, 3721, 1, 4));
    }
}
