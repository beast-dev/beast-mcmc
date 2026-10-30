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
 * normalized unless the caller implements the normalization derivative (BASTA
 * and MASCOT do, for models normalized with their frequency model), and it must
 * expose a {@link GeneralSubstitutionModel#getRatesParameter()}.
 */
public final class ComplexSubstitutionModelGradientSupport {

    private ComplexSubstitutionModelGradientSupport() {
    }

    /**
     * @return {@code null} if {@code model} is gradient-compatible, otherwise
     * a human-readable reason it is not, prefixed by {@code context}.
     */
    public static String getCompatibilityError(SubstitutionModel model, String context) {
        return getCompatibilityError(model, context, false);
    }

    /**
     * @param normalizationSupported whether the caller applies the derivative of
     *                               the normalization constant, which it can only
     *                               do for {@link #normalizesWithFrequencyModel} models
     * @return {@code null} if {@code model} is gradient-compatible, otherwise
     * a human-readable reason it is not, prefixed by {@code context}.
     */
    public static String getCompatibilityError(SubstitutionModel model, String context,
                                               boolean normalizationSupported) {
        if (!(model instanceof ComplexSubstitutionModel)) {
            return context + " currently requires ComplexSubstitutionModel-compatible rate ordering";
        }
        if (((ComplexSubstitutionModel) model).getNormalization()) {
            if (!normalizationSupported) {
                return context + " requires normalized=\"false\" until the normalization derivative is implemented";
            }
            if (!normalizesWithFrequencyModel(model)) {
                return context + " requires normalized=\"false\" for " + namedClass(model).getSimpleName() +
                        ", whose normalization is not supported by the gradient";
            }
        }
        if (getRatesParameter(model) == null) {
            return context + ": model does not expose a rates parameter";
        }
        return null;
    }

    /**
     * @return whether Q is rate * (frequency or 1) * (indicator or 1) per entry and is
     * normalized by sum_i pi_i (-Q_ii) with pi from the frequency model, so that the
     * normalization constant is linear in the rates. Subclasses that change the
     * frequencies used (e.g. the stationary distribution of Q) or the rate
     * parameterization are excluded. Anonymous subclasses, such as the XML
     * parser's, which only choose the eigen system, count as their superclass.
     */
    public static boolean normalizesWithFrequencyModel(SubstitutionModel model) {
        Class<?> type = namedClass(model);
        return type == ComplexSubstitutionModel.class || type == SVSComplexSubstitutionModel.class;
    }

    private static Class<?> namedClass(SubstitutionModel model) {
        Class<?> type = model.getClass();
        while (type.isAnonymousClass()) {
            type = type.getSuperclass();
        }
        return type;
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
