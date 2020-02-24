import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTree;
import javax.swing.tree.*;
import java.awt.Color;

/*
 * Created on 9 janv. 2007
 */
/** This is a class to put all sequences together and score alignments, 
 * then help to build the phylogenetic tree with in-results. 
 * The alignment computes distances between sequences and gives a score, 
 * all scores are put in a matrix (2-dimentionnal table). 
 * @author Gabriel Chandesris
 */
public class PhyloTools {
    /** Lists all sequences with their name. */
    private LifeSequence ls_tab[];
    /**  Working attribute to built the tree, initializes with alignment matrix scores. */
    private int scorematrix[][];
    /** The root of the tree obtained with UPGMA algorithm. */
    private DefaultMutableTreeNode rootUPGMA;
    /** The root of the tree obtained with Neighbour-Joining algorithm. */
    private DefaultMutableTreeNode rootNJ;
    /** The interface to put following traces of execution.  */
    Phylogenetik item;
    
    /** Default constructor could be used, do not forget to add 
     *  LifeSequence's in this instance of class...*/
    public PhyloTools () { 
        this.ls_tab = null;
        this.scorematrix = null; 
        this.rootUPGMA = null;
        this.rootNJ = null;
        this.item = null;
    }
    /** Constructor starting with a ready-made table of LifeSequence's.  
     * @param tab The table of LifeSequence's which is compute before to instancuate class. */
    public PhyloTools (LifeSequence tab[]) {
        this.ls_tab = tab;
        this.scorematrix = null;
        this.rootUPGMA = null;
        this.rootNJ = null;
        this.item = null;
    }
    /** Constructor starting with a ready-made table of LifeSequence's and the interface. 
     * @param tab The table of LifeSequence's which is compute before to instancuate class. 
     * @param obj The Phylogenetik interface used. */
    public PhyloTools (LifeSequence tab[],Phylogenetik obj) {
        this.ls_tab = tab;
        this.scorematrix = null;
        this.rootUPGMA = null;
        this.rootNJ = null;
        this.item = obj;
    }
    
    /** Constructor starting with the interface alone.  
     * @param obj The Phylogenetik interface used. */
    public PhyloTools (Phylogenetik obj) {
        this.ls_tab = null;
        this.scorematrix = null;
        this.rootUPGMA = null;
        this.rootNJ = null;
        this.item = obj;
    }
    
    /** Application fo the Neighbour-Joining algorithm. */
    public void buildtreeNJ() {
        // the cluster data structure to build tree, 
        // first element when finish is root for java built
        DefaultMutableTreeNode treenodes[] = 
            new DefaultMutableTreeNode[this.ls_tab.length];
        // initializing this table of treenodes
        for (int i = 0 ; i < treenodes.length ; i++) {
            treenodes[i] = new DefaultMutableTreeNode(this.ls_tab[i]);
        }
        // initializing the matrix use to built tree
        int treematrix[][] = new int[this.scorematrix.length][this.scorematrix.length];
        for (int i = 0 ; i < this.scorematrix.length ; i++) {
            for (int j = 0 ; j < this.scorematrix[i].length ; j++) {
                treematrix[i][j] = this.scorematrix[i][j];
        } }
        
        this.setMatrixInterface(treematrix,"NJ before div");
        
        // NJ matrix built here !!
        // Compute divergency
        int divergency[] = new int[treematrix.length];
        for (int i = 0 ; i < treematrix.length ; i++) {
            int sum = 0;
            for (int j = 0 ; j < treematrix[i].length ; j++) {
                sum = sum + treematrix[i][j];
            }
            divergency[i] = sum;
        }
        
        // re-compute the treematrix
        int temp = 0;
        for (int i = 0 ; i < divergency.length ; i++) {
            for (int j = 0 ; j < divergency.length ; j++) {
                treematrix[i][j] = // temp;
                    (treematrix[i][j] - 
                    		((divergency[i] + divergency[j])
                    		        /(divergency.length-2))); 
            }
        }
        
        this.setMatrixInterface(treematrix,"NJ after div");
        
        // Making the tree (NJ) : a while iterative loop
        int test = 0;
        while (treenodes.length > 1) {
            // find the minus element in the table
            int imin = 0;int jmin = 1;
            for (int i = 0 ; i < treematrix.length ; i++) {
                for (int j = 0 ; j < treematrix[i].length-1 ; j++) { //*
                    if (treematrix[i][j] < treematrix[imin][jmin]) {
                        // here + (*) = use only half the treematrix
                        if (j > i) { imin = i;jmin = j; } } } }
            // making those two elements together
            // can't be the same element / LifeSequence's / node
            // Modify the cluster : bothing the choosen elements
            DefaultMutableTreeNode newtreenodes[] = 
                new DefaultMutableTreeNode[treenodes.length-1];
            int stopvalue =  0;int stopvalu2 = 0;
            if (imin < jmin) { stopvalue = imin;stopvalu2 = jmin; }
            		else { stopvalue = jmin;stopvalu2 = imin; }
            // // copying first part
            for (int i = 0 ; i < stopvalue ; i++) {
                newtreenodes[i] = treenodes[i]; }
            // // new node = adding the two child-nodes
            newtreenodes[stopvalue] = 
                new DefaultMutableTreeNode("Dist = "+
                        this.scorematrix[imin][jmin]/2);
            newtreenodes[stopvalue].add(treenodes[imin]);
            newtreenodes[stopvalue].add(treenodes[jmin]);
            // // copying second part
            for (int i = stopvalue+1 ; i < newtreenodes.length ; i++) {
                int imodif = 0;
                if (i >= stopvalu2) { imodif=1; }
                newtreenodes[i] = treenodes[i+imodif];
            }
            treenodes = newtreenodes;
            
            //compute the new treematrix
            int newtreematrix[][] = 
                new int[treematrix.length-1][treematrix.length-1];
            // // the new line/column
            int newlinecolumn[] = new int[treematrix.length-1];
            boolean addtoi = false;boolean addtoj = false;int modif = 0;
            newlinecolumn[0] = treematrix.length-1;
            for (int i = 0 ; i < newlinecolumn.length ; i++) {
                if ((i == imin) || (i == jmin)) { modif++; }
                newlinecolumn[i] = 
                    (treematrix[imin][i]+treematrix[i][jmin])/2; }
            // // newtreematrix is here
            int modifi = 0;int modifj = 0;
            for (int i = 0 ; i < newtreematrix.length ; i++) {
                for (int j = 0 ; j < newtreematrix[i].length ; j++) {
                    if ((i == imin) || (i == jmin)) {
                        if (!addtoi) { 
                            newtreematrix[i] = newlinecolumn;
                            addtoi = true; }
                        else {                             
                            if (addtoi) { modifi = 1; } else { modifi = 0; }
	                        if (addtoj) { modifj = 1; } else { modifj = 0; }
	                        newtreematrix[i][j] = 
	                            treematrix[i+modifi][j+modifj]; }
                    } else {
                        if ((j == imin) || (j == jmin)) {
                            if (!addtoj) { 
                                for (int k = 0 ; k < newlinecolumn.length ; k++) {
                                    newtreematrix[k][j] = newlinecolumn[k]; }
                                addtoj = true; }
                            else { 
                                if (addtoi) { modifi = 1; } else { modifi = 0; }
                                if (addtoj) { modifj = 1; } else { modifj = 0; }
                                newtreematrix[i][j] = 
                                    treematrix[i+modifi][j+modifj]; }
                        } else {
                            if (addtoi) { modifi = 1; } else { modifi = 0; }
                            if (addtoj) { modifj = 1; } else { modifj = 0; }
                            newtreematrix[i][j] = 
                                treematrix[i+modifi][j+modifj]; }
                    } 
            } }
            treematrix = newtreematrix;
            this.setMatrixInterface(treematrix,"NJ step");
        }
        this.rootNJ = treenodes[0];
    }
    
    /** Application of the UPGMA Algorithm. */
    public void buildtreeUPGMA() {
        // the cluster data structure to build tree, 
        // first element when finish is root for java built
        DefaultMutableTreeNode treenodes[] = 
            new DefaultMutableTreeNode[this.ls_tab.length];
        for (int i = 0 ; i < treenodes.length ; i++) {
            treenodes[i] = new DefaultMutableTreeNode(this.ls_tab[i]);
        }
        //      May be part following to put into UPGMA constructor
        int treematrix[][] = new int[this.scorematrix.length][this.scorematrix.length];
        for (int i = 0 ; i < this.scorematrix.length ; i++) {
            for (int j = 0 ; j < this.scorematrix[i].length ; j++) {
                treematrix[i][j] = this.scorematrix[i][j];
        } }
        
        this.setMatrixInterface(treematrix,"UPGMA start");
        
        // Making the tree (UPGMA)
        int test = 0;
        while (treenodes.length > 1) {
            // find the minus element in the table
            // unless first line and column
            int imin = 0;int jmin = 1;
            for (int i = 0 ; i < treematrix.length ; i++) {
                for (int j = 0 ; j < treematrix[i].length-1 ; j++) { //*
                    if (treematrix[i][j] < treematrix[imin][jmin]) {
                        // here + (*) = use only half the treematrix
                        if (j > i) {imin = i;jmin = j; } } } }
            // making those two elements together
            DefaultMutableTreeNode newtreenodes[] = 
                new DefaultMutableTreeNode[treenodes.length-1];
            int stopvalue =  0;int stopvalu2 = 0;
            if (imin < jmin) { stopvalue = imin;stopvalu2 = jmin; }
            		else { stopvalue = jmin;stopvalu2 = imin; }
            for (int i = 0 ; i < stopvalue ; i++) {
                newtreenodes[i] = treenodes[i]; }
            newtreenodes[stopvalue] = 
                new DefaultMutableTreeNode("Dist = "+
                        this.scorematrix[imin][jmin]/2);
            newtreenodes[stopvalue].add(treenodes[imin]);
            newtreenodes[stopvalue].add(treenodes[jmin]);
            for (int i = stopvalue+1 ; i < newtreenodes.length ; i++) {
                int imodif = 0;
                if (i >= stopvalu2) { imodif=1; }
                newtreenodes[i] = treenodes[i+imodif];
            }
            treenodes = newtreenodes;
            
            //compute the new treematrix
            int newtreematrix[][] = 
                new int[treematrix.length-1][treematrix.length-1];
            // // the new line/column
            int newlinecolumn[] = new int[treematrix.length-1];
            boolean addtoi = false;boolean addtoj = false;int modif = 0;
            newlinecolumn[0] = treematrix.length-1;
            for (int i = 0 ; i < newlinecolumn.length ; i++) {
                if ((i == imin) || (i == jmin)) { modif++; }
                newlinecolumn[i] = 
                    (treematrix[imin][i]+treematrix[i][jmin])/2; }
            // // newtreematrix is here
            int modifi = 0;int modifj = 0;
            for (int i = 0 ; i < newtreematrix.length ; i++) {
                for (int j = 0 ; j < newtreematrix[i].length ; j++) {
                    if ((i == imin) || (i == jmin)) {
                        if (!addtoi) { 
                            newtreematrix[i] = newlinecolumn;
                            addtoi = true; }
                        else {                             
                            if (addtoi) { modifi = 1; } else { modifi = 0; }
	                        if (addtoj) { modifj = 1; } else { modifj = 0; }
	                        newtreematrix[i][j] = 
	                            treematrix[i+modifi][j+modifj]; }
                    } else {
                        if ((j == imin) || (j == jmin)) {
                            if (!addtoj) { 
                                for (int k = 0 ; k < newlinecolumn.length ; k++) {
                                    newtreematrix[k][j] = newlinecolumn[k]; }
                                addtoj = true; }
                            else { 
                                if (addtoi) { modifi = 1; } else { modifi = 0; }
                                if (addtoj) { modifj = 1; } else { modifj = 0; }
                                newtreematrix[i][j] = 
                                    treematrix[i+modifi][j+modifj]; }
                        } else {
                            if (addtoi) { modifi = 1; } else { modifi = 0; }
                            if (addtoj) { modifj = 1; } else { modifj = 0; }
                            newtreematrix[i][j] = 
                                treematrix[i+modifi][j+modifj]; }
                    } 
            } }
            treematrix = newtreematrix;
            this.setMatrixInterface(treematrix,"UPGMA step");
            // }//end if (imin==jmin)
        }
        // root
        this.rootUPGMA = treenodes[0];
    }
    
    
   
    
    /** This method to add a new sequence into instance of PhyloTools. 
     * @param lifeseq A LifeSequence instance already defined. */
    public void addLifeSequence(LifeSequence lifeseq) {
        if (this.ls_tab == null) { this.ls_tab = new LifeSequence[1]; } 
        else { 
            LifeSequence newtable[] = new LifeSequence[this.ls_tab.length+1];
            for (int i = 0 ; i < this.ls_tab.length ; i++) {
                newtable[i] = this.ls_tab[i];
            }
            this.ls_tab = newtable;
        }
        this.ls_tab[this.ls_tab.length-1] = lifeseq;
    }
    
    /** This method to add a new sequence into...
     * @param name Name of sequence. 
     * @param sequence The sequence itself.  */
    public void addLifeSequence(String name, String sequence) {
        this.addLifeSequence(new LifeSequence(name,sequence));
    }
    
    /** This method to get LifeSequence (name and sequence) at a specific position. */
    public LifeSequence getSequence (int i) { return ls_tab[i]; }
    
    /** This method to return the number of LifeSequence's. */
    public int getSeqNumber() { return this.ls_tab.length; }
    
    /** This method to get the last sequence added to the Phylotool */
    public LifeSequence getLastSeq () { return this.ls_tab[this.ls_tab.length-1]; }

    /** Return the matrix of scores. */
    public int[][] getMatrixScore() { return scorematrix; }
    
    /** If the interface exist, put text in it. 
     * @param txt Strinf txt to put. */
    public void setTextInterface(String txt) {
        if (this.item != null) {
            this.item.addtoaff(txt);
        }
    }
    
    /** If the interface exist, put a matrix of int in it. */
    public void setMatrixInterface(int[][] matrix,String title) {
        if (this.item != null) {
            this.item.addtoaff(title);
            this.item.addtoaff(Phylogenetik.getAffMatrix(matrix));
        }
    }
    
    /** Set the interface (use for transmitting text for the user's see). */
    public void setInterface(Phylogenetik obj) { this.item = obj; }
    
    /** Launch the methods makeAlign(), makeTrees() and showTrees(). */
    public void makeAll() {
        this.makeAlign();
        this.makeTrees();
        this.showTrees();
    }
    
    /** This method is separate from other because it will be call when
     * ALL sequences / Lifequence's were added to Phylotools. 
     * Only because Alignment cost -some- time-spend on each call 
     * (calling this method on each addLifeSequence would be loss of time).
     *  */
    public void makeAlign() {
        // DONE only half the matrix ??
        // not really entirely done but NJ and UPGMA use the inferior part
        int max = this.ls_tab.length;
        this.scorematrix = new int[max][max]; 
        for (int i = 0 ; i < max  ; i++) {
            for (int j = 0 ; j < max ; j++) {
                // TODO choice of other alignment ??
                AlignNW nw = new AlignNW(
                        this.ls_tab[i],this.ls_tab[j]);
                this.scorematrix[i][j] = nw.getScore();
            }
        }
    }
    
    /** This method build trees.  */
    public void makeTrees () {
        this.buildtreeUPGMA();
        this.buildtreeNJ();
    }
    
    /** This methid to show trees. */
    public void showTrees() {
        PhyloTools.showTree(this.rootUPGMA,"UPGMA");
        PhyloTools.showTree(this.rootNJ,"NJ");
    }
    
    /** This methode built a new JFram containing the Tree in Javax.swing.JTree format.  */
    public static void showTree(DefaultMutableTreeNode root,String title) {
        JFrame phylotreeframe = new JFrame(title+" Tree result");
        JPanel phylotreepanel = new JPanel();
        DefaultTreeModel ptmodel = new DefaultTreeModel(root);
        JTree phylotree = new JTree(ptmodel);
        TreeCellRenderer cellRenderer = phylotree.getCellRenderer();
        if (cellRenderer instanceof DefaultTreeCellRenderer) {
	          DefaultTreeCellRenderer renderer = (DefaultTreeCellRenderer)cellRenderer;
	          renderer.setBackgroundNonSelectionColor(Color.gray); 
	          renderer.setBackgroundSelectionColor(Color.black);
	          renderer.setTextSelectionColor(Color.white);
	          renderer.setTextNonSelectionColor(Color.black); 
	          phylotree.setBackground(Color.gray);
        }
        phylotree.putClientProperty("JTree.lineStyle","Horizontal");
        phylotree.setShowsRootHandles(true);
        phylotreepanel.add(phylotree);
        phylotree.setBounds(0,0,400,400);
        phylotreeframe.getContentPane().add(phylotreepanel);
        phylotreepanel.setBounds(0,0,400,400);
        phylotreeframe.setSize(400,400);
        phylotreeframe.show();
    }
    
    /** To get the number of LifeSequence's into an instance of this class. */
    public int numOfSeq() { 
        try { return ls_tab.length; }
        catch (NullPointerException e) { return 0; }
    }
    
    /** Return String overview of the object instance. 
     * @see Phylogenetik#getAffMatrix */
    public String toString() { 
		String aff_txt = "";
        for (int i = 0 ; i < this.numOfSeq() ; i++) {
            aff_txt = aff_txt+ls_tab[i].toString();
        }
        aff_txt = aff_txt + Phylogenetik.getAffMatrix(scorematrix);
        return aff_txt;
    }
}
