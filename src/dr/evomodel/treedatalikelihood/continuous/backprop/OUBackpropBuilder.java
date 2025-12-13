package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;

public final class OUBackpropBuilder {

    private final MultivariateElasticModel elasticModel;
    // TODO add also precision parametrization and possibly also mean parametrization

    public OUBackpropBuilder(OUDiffusionModelDelegate delegate) {
        this.elasticModel = delegate.getElasticModel();
    }

    public PrimitiveParameterBackpropStrategy buildPrimitiveStrategy() {

        // Strategy choice based ONLY on forward parametrization information
        if (elasticModel.isBlockDiag()) {
//            return new GeneralBackpropStrategy();
//            throw new UnsupportedOperationException(
//                    "BlockDiagonal S parameterization backprop not yet implemented."
//            );
            return new BlockDiagonalizableBackpropStrategy(elasticModel);
//            return new BlockDiagBackpropStrategy(
////                    elasticModel.getBasisR(),
////                    elasticModel.getBasisRinv(),
////                    elasticModel.getBlockDiagGenerator(),
//                    elasticModel.getBlockStarts(),
//                    elasticModel.getBlockSizes(),
//                    true
//            );
        }

        if (elasticModel.isDiagonal()) {
            throw new UnsupportedOperationException(
                    "Diagonal S parameterization backprop not yet implemented."
            );
//            return new DiagonalBackpropStrategy(
//                    elasticModel.getBasisEigenValues()
//            );
        }

        return new GeneralBackpropStrategy();

//        switch (elasticModel.getParametrization()) {
//            case AS_DECOMPOSED:
//                return new DecomposedEigenBackpropStrategy(
//                        elasticModel.getBasisR(),
//                        elasticModel.getBasisRinv(),
//                        elasticModel.getBasisEigenValues()
//                );
//
//            case GENERAL:
//                return new FullMatrixBackpropStrategy();
//        }

//        throw new IllegalStateException("Unknown OU S parametrization.");
    }
}

