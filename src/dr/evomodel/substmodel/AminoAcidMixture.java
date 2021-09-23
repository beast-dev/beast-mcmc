package dr.evomodel.substmodel;

import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;

import java.util.List;

public class AminoAcidMixture extends DesignMatrix {

    public AminoAcidMixture(List<EmpiricalAminoAcidModel> rateMatrix) {
        super("name", getParameters(rateMatrix), false);
    }

    private static Parameter[] getParameters(List<EmpiricalAminoAcidModel> rateMatrix) {
        Parameter[] p = new Parameter[rateMatrix.size()];

        for (int i = 0; i < p.length; ++i) {

            EmpiricalAminoAcidModel model = rateMatrix.get(i);
            double[] rates = model.getEmpiricalRateMatrix().getEmpiricalRates(); // this is flat and includes diaganoal elements
            // we want flattened upper triangular ordered, then flatten lower triangular order
            p[i] = new Parameter.Default(rates); // this is wrong, because we need to permute the elements in rates


            // TODO
        }
        return p;
    }


}
