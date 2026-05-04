package dr.evomodel.treedatalikelihood.continuous.canonical;

/**
 * Utility for capability-based dispatch on {@link CanonicalBranchTransitionProvider} instances.
 *
 * <p>Prefer {@link #requireCapability} over the older {@link #requireOUProvider}: it names
 * both the provider class and the missing capability in its error message, making failures
 * easier to diagnose.
 */
public final class CanonicalOUProviderSupport {

    private CanonicalOUProviderSupport() { }

    /**
     * Casts {@code provider} to {@code capabilityType}, or throws a descriptive
     * {@link UnsupportedOperationException} naming the provider class and the required
     * capability if the cast is not possible.
     *
     * @param provider       the object to inspect
     * @param capabilityType the capability interface that must be implemented
     * @param <T>            the capability type
     * @return {@code provider} cast to {@code T}
     * @throws UnsupportedOperationException if {@code provider} does not implement
     *                                        {@code capabilityType}
     */
    public static <T> T requireCapability(final Object provider, final Class<T> capabilityType) {
        if (!capabilityType.isInstance(provider)) {
            throw new UnsupportedOperationException(
                    provider.getClass().getSimpleName()
                    + " does not implement " + capabilityType.getSimpleName()
                    + "; this capability is required for this operation.");
        }
        return capabilityType.cast(provider);
    }

    /**
     * Casts {@code transitionProvider} to {@link CanonicalOUTransitionProvider}.
     *
     * @deprecated Use {@link #requireCapability(Object, Class)} instead, which produces a
     *     more descriptive error message including the provider's class name.
     */
    @Deprecated
    public static CanonicalOUTransitionProvider requireOUProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        return requireCapability(transitionProvider, CanonicalOUTransitionProvider.class);
    }
}
