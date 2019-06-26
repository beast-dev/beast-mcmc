package dr.evomodel.substmodel;

import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import optimization.functionImplementation.ObjectiveFunctionNonLinear;
import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import solvers.NonlinearEquationSolver;


public class CorrectedF3x4 extends FrequencyModel {

    private final FrequencyModel nucleotideFrequencyModel;
    private final Codons dataType;
    private final Parameter codonFrequencies;

    public CorrectedF3x4(Codons dataType,
                         FrequencyModel nucleotideFrequencyModel,
                         Parameter codonFrequencies) {
        super(dataType, codonFrequencies); // TODO To make updateFrequencyParameter() lazy, codonFrequencies should be a within-class proxy

        if (nucleotideFrequencyModel.getDataType() != Nucleotides.INSTANCE) {
            throw new IllegalArgumentException("Must provide a nucleotide frequency model");
        }

        this.codonFrequencies = codonFrequencies;
        this.dataType = dataType;
        this.nucleotideFrequencyModel = nucleotideFrequencyModel;
        updateFrequencyParameter();
        addModel(nucleotideFrequencyModel);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        assert (model == nucleotideFrequencyModel);

        updateFrequencyParameter();

        fireModelChanged(model);
    }

    private DMatrixRMaj nonlinearcf3x4() {

        double[] freq3x4 = nucleotideFrequencyModel.getFrequencies();
        double pi1_A=freq3x4[0]; double pi1_C=freq3x4[1]; double pi1_G=freq3x4[2]; double pi1_T=freq3x4[3];
        double pi2_A=freq3x4[4]; double pi2_C=freq3x4[5]; double pi2_G=freq3x4[6]; double pi2_T=freq3x4[7];
        double pi3_A=freq3x4[8]; double pi3_C=freq3x4[9]; double pi3_G=freq3x4[10]; double pi3_T=freq3x4[11];


        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(9, 1);
                double x0 = x.get(0); double x1 = x.get(1); double x2 = x.get(2); double x3 = x.get(3);
                double x4 = x.get(4); double x5 = x.get(5); double x6 = x.get(6); double x7 = x.get(7);
                double x8 = x.get(8);
                double phi_9 = (1-(x0+x3+x6));
                double phi_den = (1-(phi_9*x1*x2+phi_9*x1*x8+phi_9*x7*x2));

                f.set(0,0,(x0)/phi_den - pi1_A);
                f.set(1,0,(x3)/phi_den - pi1_C);
                f.set(2,0,(x6)/phi_den - pi1_G);
                f.set(3,0,(x1*(1-(phi_9*x2+phi_9*x8)))/phi_den - pi2_A);
                f.set(4,0,(x4)/phi_den - pi2_C);
                f.set(5,0,(x7*(1-(phi_9*x2)))/phi_den - pi2_G);
                f.set(6,0,(x2*(1-(phi_9*x1+phi_9*x7)))/phi_den - pi3_A);
                f.set(7,0,(x5)/phi_den - pi3_C);
                f.set(8,0,(x8*(1-(phi_9*x1)))/phi_den - pi3_G);
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                return null;
            }

        };
        //initial guess
        DMatrixRMaj initialGuess = new DMatrixRMaj(9, 1);
        for (int i = 0; i < 9; i++) {
            initialGuess.set(i, 0.25);
        }
        //options
        Options options = new Options(9);
        options.setAnalyticalJacobian(false);
        options.setAlgorithm(Options.LINE_SEARCH);
        options.setSaveIterationDetails(false);
        options.setAllTolerances(1e-12);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        //solve and print output
        nonlinearSolver.solve(new DMatrixRMaj(initialGuess));

        DMatrixRMaj nonlinearSolverx = nonlinearSolver.getX();
        return nonlinearSolverx;
    }


    private void updateFrequencyParameter() {
        DMatrixRMaj result = nonlinearcf3x4();
        double phi_0 = result.get( 0,0); //a1
        double phi_1 = result.get( 1,0); //a2
        double phi_2 = result.get( 2,0); //a3
        double phi_3 = result.get( 3,0); //c1
        double phi_4 = result.get( 4,0); //c2
        double phi_5 = result.get( 5,0); //c3
        double phi_6 = result.get( 6,0); //g1
        double phi_7 = result.get( 7,0); //g2
        double phi_8 = result.get( 8,0); //g3
        double phi_9 = 1 - (phi_0 + phi_3 + phi_6); //t1
        double phi_10 = 1 - (phi_1 + phi_4 + phi_7); //t1
        double phi_11 = 1 - (phi_2 + phi_5 + phi_8); //t1

        // TODO Maybe make this lazy if necessary

        for (int nuc1 = 0; nuc1 < 4; ++nuc1) {
            double freq1;
            if (nuc1 == 0){
                freq1 = phi_0;
            } else if (nuc1 == 1){
                freq1 = phi_3;
            } else if (nuc1 == 2) {
                freq1 = phi_6;
            } else {
                freq1 = phi_9;
            }
            for (int nuc2 = 0; nuc2 < 4; ++nuc2) {
                double freq2;
                if (nuc2 == 0){
                    freq2 = phi_1;
                } else if (nuc2 == 1){
                    freq2 = phi_4;
                } else if (nuc2 == 2) {
                    freq2 = phi_7;
                } else {
                    freq2 = phi_10;
                }
                for (int nuc3 = 0; nuc3 < 4; ++nuc3) {
                    double freq3;
                    if (nuc3 == 0){
                        freq3 = phi_2;
                    } else if (nuc3 == 1){
                        freq3 = phi_5;
                    } else if (nuc3 == 2) {
                        freq3 = phi_8;
                    } else {
                        freq3 = phi_11;
                    }

                    int state = dataType.getState(nuc1, nuc2, nuc3);

                    if (!(dataType.isStopCodon(state))) {
                        frequencyParameter.setParameterValue( state,
                                freq1 * freq2 * freq3);
                    }
                }
            }
        }

        final double sum = getSumOfFrequencies(frequencyParameter);

        for (int i = 0; i < frequencyParameter.getSize(); i++) {

            frequencyParameter.setParameterValue(i, frequencyParameter.getParameterValue(i) / sum);

        }
    }

    private String getDimensionName(int dim) {
        return codonFrequencies.getParameterName() + "." + dataType.getTriplet(dim);
    }

}

