package dr.app.beauti.components.continuous;

public enum ContinuousModelExtensionType {
    RESIDUAL("Residual variance extension"),
    LATENT_FACTORS("Latent factor model");

    ContinuousModelExtensionType(String name) {
        this.name = name;
    }

    public String toString() {
        return name;
    }

    private final String name;

}
