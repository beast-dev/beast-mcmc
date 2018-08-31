package dr.evomodel.treedatalikelihood;

import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;

/**
 * @author Marc A. Suchard
 */
public interface PostOrderStatistics {

    void executePostOrder(TreeDataLikelihood treeDataLikelihood);

    enum Discrete implements PostOrderStatistics {

        FULL {
            @Override
            public void executePostOrder(TreeDataLikelihood treeDataLikelihood) {
                treeDataLikelihood.getLogLikelihood();
            }
        }
    }

    enum Continuous implements PostOrderStatistics {

        NONE {

            @Override
            public void executePostOrder(TreeDataLikelihood treeDataLikelihood) {
               // Do nothing
            }
        },

        MINIMAL {

            @Override
            public void executePostOrder(TreeDataLikelihood treeDataLikelihood) {

                ContinuousDataLikelihoodDelegate delegate = (ContinuousDataLikelihoodDelegate)
                        treeDataLikelihood.getDataLikelihoodDelegate();

                delegate.setComputeWishartStatistics(false);
                delegate.setComputeRemainders(false);

                try {
                    treeDataLikelihood.executePostOrderComputation();
                } catch (DataLikelihoodDelegate.LikelihoodException e) {
                    throw new RuntimeException(e.getMessage());
                }

                delegate.setComputeWishartStatistics(true);
                delegate.setComputeRemainders(true);

            }
        },

        REMAINDERS {

            @Override
            public void executePostOrder(TreeDataLikelihood treeDataLikelihood) {
                treeDataLikelihood.getLogLikelihood();
            }
        },

        REMAINDERS_AND_WISHART {

            @Override
            public void executePostOrder(TreeDataLikelihood treeDataLikelihood) {
                treeDataLikelihood.getLogLikelihood();
            }
        }
    }
}
