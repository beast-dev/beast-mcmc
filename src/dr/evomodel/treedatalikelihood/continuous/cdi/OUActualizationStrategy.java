package dr.evomodel.treedatalikelihood.continuous.cdi;

interface OUActualizationStrategy {

    void setDiffusionStationaryVariance(SafeMultivariateActualizedWithDriftIntegrator integrator,
                                        int precisionIndex,
                                        double[] basisValues,
                                        double[] basisRotations);

    void updateOrnsteinUhlenbeckDiffusionMatrices(SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                  int precisionIndex,
                                                  int[] probabilityIndices,
                                                  double[] edgeLengths,
                                                  double[] optimalRates,
                                                  double[] basisValues,
                                                  double[] basisRotations,
                                                  int updateCount);

    void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                            int precisionIndex,
                                                            int[] probabilityIndices,
                                                            double[] edgeLengths,
                                                            double[] optimalRates,
                                                            double[] basisValues,
                                                            double[] basisRotations,
                                                            int updateCount);
}
