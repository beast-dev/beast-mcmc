package dr.app.beauti.generator;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.XMLWriter;

/**
 * This is an abstract base class for component generators. Currently it simply
 * wraps each insertion with a pair of comments.
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class BaseComponentGenerator extends Generator implements ComponentGenerator {
    protected BaseComponentGenerator(final BeautiOptions options) {
        super(options, null);
    }

    public void generateAtInsertionPoint(final InsertionPoint point, final Object item, final XMLWriter writer) {
        writer.writeComment("START " + getCommentLabel());
        generate(point, item, writer);
        writer.writeComment("END " + getCommentLabel());
        writer.writeText("");
    }

    protected abstract void generate(final InsertionPoint point, final Object item, final XMLWriter writer);

    protected abstract String getCommentLabel();
}
