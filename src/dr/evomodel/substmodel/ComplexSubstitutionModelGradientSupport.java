/*
 * ComplexSubstitutionModelGradientSupport.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.substmodel;

import dr.inference.model.Parameter;

/**
 * Shared compatibility-check and rates-parameter extraction for gradients
 * wrt a {@link SubstitutionModel} used as a migration-rate process, used by
 * both MASCOT's {@code MascotDynamics} and BASTA's {@code
 * StructuredCoalescentLikelihoodGradient}. Both engines need the same three
 * things from the model before they can map an internal rate-matrix gradient
 * back to an XML-facing {@link Parameter}: it must be a {@link
 * ComplexSubstitutionModel} (rate-ordering assumption), it must not be
 * normalized (the normalization derivative isn't implemented on either side),
 * and it must expose a {@link GeneralSubstitutionModel#getRatesParameter()}.
 */
public final class ComplexSubstitutionModelGradientSupport {

    private ComplexSubstitutionModelGradientSupport() {
    }

    /**
     * @return {@code null} if {@code model} is gradient-compatible, otherwise
     * a human-readable reason it is not, prefixed by {@code context}.
     */
    public static String getCompatibilityError(SubstitutionModel model, String context) {
        if (!(model instanceof ComplexSubstitutionModel)) {
            return context + " currently requires ComplexSubstitutionModel-compatible rate ordering";
        }
        if (((ComplexSubstitutionModel) model).getNormalization()) {
            return context + " requires normalized=\"false\" until the normalization derivative is implemented";
        }
        if (getRatesParameter(model) == null) {
            return context + ": model does not expose a rates parameter";
        }
        return null;
    }

    /**
     * @return the model's rates {@link Parameter}, or {@code null} if it
     * isn't a {@link GeneralSubstitutionModel} (and so has no such parameter
     * to return).
     */
    public static Parameter getRatesParameter(SubstitutionModel model) {
        if (model instanceof GeneralSubstitutionModel) {
            return ((GeneralSubstitutionModel) model).getRatesParameter();
        }
        return null;
    }
}
