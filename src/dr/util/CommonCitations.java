package dr.util;

/**
 * @author Marc Suchard
 *         <p/>
 *         A class to hold common BEAST model and inference machinery citations.
 *         Avoids repeated definitions in classes that share citations
 */
public class CommonCitations {

    public static Citation LEMEY_2010 = new Citation(
            new Author[]{
                    new Author("P", "Lemey"),
                    new Author("A", "Rambaut"),
                    new Author("JJ", "Welch"),
                    new Author("MA", "Suchard")
            },
            "Phylogeography takes a relaxed random walk in continous space and time",
            "Molecular Biology and Evolution",
            Citation.Status.ACCEPTED
    );

    public static Citation OBRIEN_2009 = new Citation(
            new Author[]{
                    new Author("JB", "O'Brien"),
                    new Author("VN", "Minin"),
                    new Author("MA", "Suchard")

            },
            "Learning to count: robust estimates for labeled distances between molecular sequences",
            2009,
            "Molecular Biology and Evolution",
            26,
            801, 814,
            Citation.Status.PUBLISHED
    );
}
