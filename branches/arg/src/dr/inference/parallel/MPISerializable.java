package dr.inference.parallel;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 23, 2007
 * Time: 9:12:11 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MPISerializable {

	public void sendState(int toRank);

	public void receiveState(int fromRank);

}
