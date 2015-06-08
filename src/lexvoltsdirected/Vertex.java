package lexvoltsdirected;

/**
 * Created by sushant on 1/31/15.
 */
public class Vertex {
    // An in-edge is one that has the vertex as its dest
    // An out-edge is one that has the vertex as its source

    public int deg; // Denotes the current degree -- counts edges from all pressure levels, but does not count terminal-terminal edges
    public int actDeg; // Denotes the current active degree. Does not count terminal-terminal edges

    public boolean isTerminal;
    public boolean isOrigTerminal;
    public double voltage;

    public int actNbrsEnd; // Gives the end of the active neighbors of each vertex
    // NOTE: This is the index of the last active entry

    double vlow;
    double vhigh;

    double highTermVolt;
    double lowTermVolt;
    double highTermDist;
    double lowTermDist;
}