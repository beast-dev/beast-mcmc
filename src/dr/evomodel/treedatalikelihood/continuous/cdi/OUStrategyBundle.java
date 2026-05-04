package dr.evomodel.treedatalikelihood.continuous.cdi;

public final class OUStrategyBundle {

    public enum Kind {
        DIAGONAL,
        GENERAL,
        BLOCK
    }

    private final Kind kind;
    private final OUGradientStrategy gradientStrategy;

    OUStrategyBundle(final Kind kind, final OUGradientStrategy gradientStrategy) {
        this.kind = kind;
        this.gradientStrategy = gradientStrategy;
    }

    public Kind getKind() {
        return kind;
    }

    public OUGradientStrategy getGradientStrategy() {
        return gradientStrategy;
    }

    public void applyTo(final SafeMultivariateActualizedWithDriftIntegrator integrator) {
        switch (kind) {
            case DIAGONAL:
                integrator.useDiagonalOUActualizationStrategy();
                return;
            case BLOCK:
                integrator.useBlockOUActualizationStrategy();
                return;
            case GENERAL:
            default:
                integrator.useGeneralOUActualizationStrategy();
        }
    }
}
