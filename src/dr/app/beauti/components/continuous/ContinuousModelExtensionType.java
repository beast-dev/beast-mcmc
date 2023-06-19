package dr.app.beauti.components.continuous;

public enum ContinuousModelExtensionType {
    NONE("None"),
    RESIDUAL("Residual variance"),
    LATENT_FACTORS("Latent factor model");

    ContinuousModelExtensionType(String name) {
        this.name = name;
    }

    public String toString() {
        return name;
    }

    private final String name;

}
