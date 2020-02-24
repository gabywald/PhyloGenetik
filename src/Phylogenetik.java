/*
 * Created on 11 janv. 2007
 */

/** This is the main class of the application, method main 
 * and graphical interface with users. 
 * @author Gabriel Chandesris
 */
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
public class Phylogenetik extends JFrame {
    /** This attribute contain sequences, table of alignment and score matrix */
    private PhyloTools buildingtool;
    //Graphical attributes
    /** Define the font appearance of graphical interface */
    private Font font = new Font("monaco", Font.BOLD, 20);
    /** This button is for running the program on the selected files/sequences.  */
	private JButton runbutton = new JButton("Run building tree");
	/** This button to add a new sequence to analyse.  */
	private JButton addbutton = new JButton("Add a sequence");
	/** This button to add new sequenceS to analyse.  */
	private JButton addmorebutton = new JButton("Add several sequences");
	/** This button to add new sequenceS contained by a directory.  */
	private JButton adddirbutton = new JButton("Add a group (directory)");
	/** On this text area to follow the work in progress... */
	private JTextArea aff_action = new JTextArea("",33,65);
	/** Making scrolling on the main text area of the interface. */
	private JScrollPane aff_scroll;
	/** List of the sequences. */
	private JTextArea list_sequences = new JTextArea("",5,30);

	/** The constructor of the graphic interface of this program. 
	 * Construction of the staic instance and graphical view. */
    public Phylogenetik () {
        // Most important thing here...
        // Construct the static instance
        this.buildingtool = new PhyloTools(this);
        
        // drawing the interface... two panels / JPanels
        JPanel buttonspanel = new JPanel() {
		    	public void paint(Graphics g) {
		    		super.paint(g);
		    		g.setColor(Color.black);
		    		g.setFont(new Font ("serif", Font.ITALIC+Font.BOLD, 20));
		    		g.drawString("Phylogenetic's",10,20);
		    		g.drawString("Build Tree",10,40);
		    		g.drawString("List of sequences : ",5,140);
		    	}
	    	};
        JPanel aff_panel = new JPanel();
        
        // the left part of the frame
        buttonspanel.setLayout(null);
        buttonspanel.add(this.addbutton);
        this.addbutton.setBounds(0,50,200,25);
        this.addbutton.addActionListener(new ChooseActionListen(this,1));
        buttonspanel.add(this.addmorebutton);
        this.addmorebutton.setBounds(0,75,200,25);
        this.addmorebutton.addActionListener(new ChooseActionListen(this,20));
        buttonspanel.add(this.adddirbutton);
        this.adddirbutton.setBounds(0,100,200,25);
        this.adddirbutton.addActionListener(new ChooseDirActionListen(this));
        buttonspanel.add(this.list_sequences);
        this.list_sequences.setBounds(15,150,170,400);
        this.list_sequences.setFont(new Font("times",Font.BOLD,12));
        this.list_sequences.setEditable(false);
        
        // the right part of the frame (textearea, scroll and run button)
        aff_panel.setLayout(null);
        Dimension pref = new Dimension(580,530);
        this.aff_action.setPreferredSize(pref);
        this.aff_action.getScrollableBlockIncrement(new Rectangle(pref),
                SwingConstants.VERTICAL,
                SwingConstants.HORIZONTAL);
        this.aff_scroll = new JScrollPane(this.aff_action
                ,JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED
                ,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        aff_panel.add(this.aff_scroll);
        // //aff_action.setBounds(5,5,590,545);
        this.aff_scroll.setBounds(5,5,590,545);
        this.aff_action.setFont(new Font("courier",Font.BOLD,12));
        this.aff_action.setEditable(false);
        aff_panel.add(runbutton);
        this.runbutton.setBounds(0,550,600,25);
        this.runbutton.addActionListener(new RunActionListen(this));
        // all is put in place
        Container global = this.getContentPane();
	    global.setLayout(null);
	    global.add(buttonspanel);
	    buttonspanel.setBounds(0,0,200,600);
	    global.add(aff_panel);
	    aff_panel.setBounds(200,0,800,600);
	    this.setSize(800,600);
	    this.setLocationRelativeTo(null);
	    this.setResizable(false);
	    this.setTitle("Phylogenetik");
	    this.show();
	    
	    /* First text to appear : a little help ! */
	    this.aff_action.setText(""+
	    "Welcome to this application -- Phylogenetik"+"\n"+
	    "\tThis software works with fasta sequences (nucleic or proteic)."+"\n"+
	    "\tUse of this software is easy : "+"\n"+
	    "\tClick on thoses buttons on left side of the window,"+"\n"+
	    "\n<------ This button to choose one sequence at a time,"+"\n"+
	    "\n<------ This button to choose a set of several sequences (less than 20),"+"\n"+
	    "\n<------ This button to select sequences in a complete directory."+"\n"+
	    "\t(Just select a sequence in the directory you want completely add. )"+"\n"+
	    "\n\n Name of sequences, identifiers and sequences will appear under this text."+"\n"+
	    "\n\n Identifiers will be added to the area on the left (list). "+"\n"+
	    "<--------- identifiers listed in this area"+"\n"+
	    "\n\n\nWhen you have finish to add sequences, just click on the button on the \n" +
	    		"bottom of this window. (Run building tree)"+"\n"+
	    	"\n=>Files of sequences must be in FASTA format (and filename ends with .fasta)."+"\n"+
	    	"\n=>Take care of your sequences, add of the same type : only nucleic OR only proteic"+"\n"+
	    "\n--------------------------------------------------------------------------\n");
	    
	    /* For / When testing, put code here ! */
    }
    
    /** The method to add a file to the PhyloTools object. */
    public void addFile(File toadd) {
        //choosenfile = chooser.getSelectedFile().getName();
        //File newsequence = chooser.getSelectedFile();
        BufferedReader sequencereader = readFile(toadd);
        String line = "test",name = "",sequence = "";
        // FASTA format
        /* >name_of_sequence and other references
         * sequence after line-break or end-of-line
         */
        while (line != null) {
            try { 
                line = sequencereader.readLine(); 
                if (line != null) { 
                    // Here is sequence and its name
                    line = line.replaceAll("\n","");
                    line = line.replaceAll("\r","");
                    if (line.startsWith(">")) { name = line; }
                    else { sequence = sequence + line; }
                }
            }
            catch (IOException e) { 
                this.addtoaff(toadd.getName()+" : Problem to get the content of the file..."); }
        	}
        
        if ((!name.equals("")) && (!sequence.equals(""))) {
            // add into Phylotools.buildingtool
            this.addLifeSequence(name,sequence);
        }
        
        try { sequencereader.close(); }
        catch (IOException e) { this.addtoaff(toadd.getName()+" : Problem to close the file..."); }
	    if (!toadd.getName().equals("")) {
	        this.addtoaff(toadd.getName());
	        this.addtoseq(this.getTool().getLastSeq().getName());
	        this.addtoaff("Name : ");
	        this.addtoaff(this.getTool().getLastSeq().getName());
	        this.addtoaff("Sequence : ");
	        this.addtoaff(this.getTool().getLastSeq().getSequence());
	        this.addtoaff("--------------------\n\n");
	        this.addtoaff("Number of sequences : "+this.getTool().numOfSeq());
	    }
    }
    
    /** The method to run when user want to add more than one sequence (i.e. several). */
    public int chooseFileNumber(int a) {
        int filenumber = 1;
        DialogFileNumber dialog = new DialogFileNumber(this,a);
        filenumber = dialog.getInt();
        return filenumber;
    }
    
    /** The method to add a group of sequences : select a directory. */
    public void chooseDir() {
        JFileChooser chooser = new JFileChooser();
        int returnVal = chooser.showOpenDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            File directory = chooser.getCurrentDirectory();
            File listfiles[] = directory.listFiles();
            for (int i = 0 ; i < listfiles.length ; i++) {
                if ((listfiles[i].isFile()) && 
                        (listfiles[i].getName().endsWith(".fasta"))) {
                    this.addFile(listfiles[i]); }
            }
        }
    }
    
    /** The method run when user wants to add a new sequence */
    public String chooseFile() {
        String choosenfile = "";
        JFileChooser chooser = new JFileChooser();
        int returnVal = chooser.showOpenDialog(this);
        if ( (returnVal == JFileChooser.APPROVE_OPTION) &&
                (chooser.getSelectedFile().getName().endsWith(".fasta")) ){
            choosenfile = chooser.getSelectedFile().getName();
            File newsequence = chooser.getSelectedFile();
            this.addFile(newsequence);
        }
        return choosenfile;
    }
    
    /** Specific method to read external files. */
    private BufferedReader readFile (File data) {
        BufferedReader tampon;
        try {
            tampon = new BufferedReader (new FileReader(data));
        } catch (IOException e) {
            this.addtoaff("Problem to find the file... ");
            tampon = null;
        }
        return tampon;
    }
    
    /** The method to add a new LifeSequence, with an instance of LifeSequence.
     * @see PhyloTools#addLifeSequence(LifeSequence)
     * @param lifeseq An instanciation of LifeSequence.  */
    public void addLifeSequence(LifeSequence lifeseq) {
        this.buildingtool.addLifeSequence(lifeseq);
    }
    
    /** The method to add a new LifeSequence, with name and sequence. 
     * @see PhyloTools#addLifeSequence
     * @param name Name of the sequence.
     * @param seq The sequence.
     */
    public void addLifeSequence(String name,String seq) {
        this.buildingtool.addLifeSequence(name,seq);
    }
    
    /** The method to make the alignment. 
     * @see PhyloTools#makeAlign */
    public void makeAlign() { this.buildingtool.makeAlign(); }
  
    /** The method to build the trees. 
     * @see PhyloTools#makeTrees */
    public void makeTree() { this.buildingtool.makeTrees(); }
    
    /** The method to show the trees. 
     * @see PhyloTools#showTrees */
    public void showTrees() { this.buildingtool.showTrees(); }
    
    /** The method to make all things we need in this application (alignment, trees and showing the trees). 
     * @see PhyloTools#makeAll */
    public void makeAll() { this.buildingtool.makeAll(); }
    
    /** Return the matrix score resulting from alignment. 
     * @see PhyloTools#getMatrixScore */
    public int[][] getMatrixScore() { return this.buildingtool.getMatrixScore(); }
    
    /** Return the form String of this object.  */
    public String toString() { // == getAff() 
        return "Building Tool : \n"+this.buildingtool.toString();
    }
    
    /** Return the matrix of score to be put into the textarea of interface.  */
    public static String getAffMatrix(int[][] matrix) { 
        String aff_txt = "";
        String linebreak = " -";
        int k = 0;
        while (k < matrix.length) { 
            linebreak = linebreak+"-------";
            k++; 
        }
        aff_txt = aff_txt+"\n"+linebreak+"\n |";
        for (int i = 0 ; i < matrix.length ; i++) {
            for (int j = 0 ; j < matrix[i].length ; j++) {
                String temp = "";
                String insert = "|";
                // if ( (i == 0) || (j == 0) ) {
                temp = matrix[i][j]+temp; 
                // } else {
                //    temp = scorematrix[i][j].getScore()+temp;
                // }
                if (temp.length() >= 6) { insert = temp+insert; }
                if (temp.length() == 5) { insert = " "+temp+insert; }
                if (temp.length() == 4) { insert = " "+temp+" "+insert; }
                if (temp.length() == 3) { insert = "  "+temp+" "+insert; }
                if (temp.length() == 2) { insert = "  "+temp+"  "+insert; }
                if (temp.length() == 1) { insert = "   "+temp+"  "+insert; }
                if (temp.length() == 0) { insert = "   "+temp+"   "+insert; }
                aff_txt = aff_txt+insert;
            }
            aff_txt = aff_txt+"\n"+linebreak+"\n |";
        }
        return aff_txt; }
    
    /** This method has to add arguments to the text area to follow actions 
     * made by the program and see some parts of the results. */
    public void addtoaff (String txt) {
        int ancientnblines = this.aff_action.getLineCount();
        this.aff_action.append(txt+"\n");
        int nblines = this.aff_action.getLineCount();
        if (nblines>40) {
	        Dimension d = this.aff_action.getPreferredSize();
	        // here count of nblines is useful how to grow the textarea
	        d.setSize(d.getWidth(),d.getHeight()+13.0*(nblines-ancientnblines));
	        // (one line + one interline)*(number of new lines) is added
	        this.aff_action.setPreferredSize(d);
        }
        this.aff_action.setCaretPosition(this.aff_action.getDocument().getLength());
    }
    
    /** The method to add name of sequences in the viewing interface (left part). */
    public void addtoseq (String txt) { 
        txt.replaceAll("\n","");
        this.list_sequences.append(txt+"\n"); }
    
    /** The method to run the alignments and building trees.
     * @see Phylogenetik#makeAlign
     * @see Phylogenetik#makeTree 
     * @see Phylogenetik#showTrees*/
    public void runprogram() {
        if (this.buildingtool.getSeqNumber() >= 3) {
	        this.addtoaff("-----");
	        this.addtoaff("Alignment Running... Please wait...");
	        this.buildingtool.makeAlign();
	        this.addtoaff(this.buildingtool.toString());
	        this.addtoaff("\n-----");
	        this.addtoaff("MakingTrees.... Please wait...");
	        this.buildingtool.makeTrees();
	        this.addtoaff("-----\n-----");
	        this.buildingtool.showTrees();
        } else {
            this.addtoaff("3 sequences are required to run...");
        }
    }
    
    /** This method to write a string argument into a file called "export-test.txt".
     * @see Phylogenetik#addtoaff
     * @param txt The String you want to export. 
     */
    public void export(String txt) {
        try { 
            PrintWriter exittxt = new PrintWriter 
        			(new FileWriter("export-test.txt"));
            exittxt.println(txt);
            exittxt.close();
        } catch (IOException e) { this.addtoaff("Export failed"); }
    }
    
    /** This method to get the PhyloTools from this instance. 
     * Other methods are more useful and returns values of this object 
     * or runs from an instance of Phylogenetik. 
     * @see Phylogenetik#getMatrixScore 
     * @see Phylogenetik#makeAlign 
     * @see Phylogenetik#addLifeSequence */
    public PhyloTools getTool() { return buildingtool; }
    
    /** This class to get action of the user to add a sequence.  */
	private class ChooseActionListen implements ActionListener {
	    /** To remember in which area to show the results. */
	    private Phylogenetik remainder;
	    /** How many time it will be asked. */
	    int ask;
	    /** The constructor which needs the object which event is attached. 
	     * @param o Phylogenetik interface. 
	     * @param i The number of times user will be asked for a sequence. */
	    public ChooseActionListen(Phylogenetik o,int i) { 
	        this.remainder = o;
	        this.ask = i;
	    }
	    /** Event capture. 
	     * @see Phylogenetik#chooseFile*/
	    public void actionPerformed (ActionEvent e) {
        		// if choose of multiple
	        if (this.ask != 1) {
	            this.ask = chooseFileNumber(this.ask);
	        }
	        if (this.ask == 0) { remainder.addtoaff("Choose number is canceled."); }
	        else {
		        for (int i = 0 ; i < this.ask ; i++) {
			        // DONE chooseFile
			        String file_chosen = remainder.chooseFile();
			        if (file_chosen.equals("")) { remainder.addtoaff("Add is canceled."); }
		        }
	        }
	    }
	}
    /** This class to get action of the user to add a group of sequences. */
	private class ChooseDirActionListen implements ActionListener {
	    /** To remember in which area to show the results. */
	    private Phylogenetik remainder;
	    /** The constructor which needs the object which event is attached. */
	    public ChooseDirActionListen (Phylogenetik o) { this.remainder = o; }
	    /** Event capture. 
	     * @see Phylogenetik#chooseDir */
	    public void actionPerformed (ActionEvent e) { remainder.chooseDir(); }
	}
	/** When user wants to run the program on the sequences.  */
	private class RunActionListen implements ActionListener {
	    /** To remember in which area to show the results. */
	    private Phylogenetik remainder;
	    /** The constructor which needs the object where event is attached */
	    public RunActionListen(Phylogenetik o) { this.remainder = o; }
	    /** Event capture.  */
	    public void actionPerformed (ActionEvent e) {
	        // DONE Running the program / algorithm itself !!
            this.remainder.addtoaff("Running... please wait !");
            this.remainder.runprogram();
	    }
	}
    /** A Dialog box if the user want to ask more than one sequence. */
	private class DialogFileNumber extends JDialog implements ActionListener {
	    /** The interface object Phylogenetik. */
	    private Phylogenetik item;
	    /** value choosen by user. */
	    private int selectvalue;
	    /** List of can-be values. */
	    private Integer numbers[];
	    /** The box which appears.  */
	    private JComboBox numberlist;
	    /** Button of the JComboBox (numberlist). */
	    private JButton ok_button,cancel_button;
	    /** The text which appears on the box.  */
	    private JLabel seq_txt;
	    /** Constructor of this class
	     * @param ancestor The phylogenetik application interface.  
	     * @param a The number of elements in the list where to choose. */
	    public DialogFileNumber (Phylogenetik ancestor,int a) {
	        super(ancestor,"Choose a number",true);
	        this.item = ancestor;
	        this.selectvalue = 0;//default value if closed
	        this.setSize(250,100);
	        this.setResizable(false);
	        this.setLocation(350,300);
	        this.numbers = new Integer[a];
	        for (int i = 0 ; i < this.numbers.length ; i++) {
	            this.numbers[i] = new Integer(i+1);
	        }
	        this.numberlist = new JComboBox(numbers);
	        this.ok_button = new JButton("OK");
	        this.cancel_button = new JButton("Cancel");
	        this.seq_txt = new JLabel("sequence(s) to add");
	        Container global = this.getContentPane();
	        global.add(this.numberlist);
	        global.add(this.seq_txt);
	        global.add(this.ok_button);
	        global.add(this.cancel_button);
	        global.setLayout(null);
	        this.numberlist.setBounds(10,10,100,20);
	        this.seq_txt.setBounds(115,15,135,20);
	        this.ok_button.setBounds(20,50,50,25);
	        this.cancel_button.setBounds(75,50,100,25);
	        this.ok_button.addActionListener(this);
	        this.cancel_button.addActionListener(this);
	        this.setVisible(true);
	    }
	    /** Event capture.  */
	    public void actionPerformed(ActionEvent e) {
	        if (e.getSource() == this.ok_button) {
	            String value = numberlist.getSelectedItem().toString();
	            selectvalue = Integer.parseInt(value);
	            this.setVisible(false);
	            this.item.addtoaff("OK--"+selectvalue);
	        }
	        if (e.getSource() == this.cancel_button) {
	            selectvalue = 0;
	            this.setVisible(false);
	        }
	    }
	    /** Getting the choosen element.  */
	    public int getInt() { return selectvalue; }
	}
	
	
	/** The main method : program starts here */
    public static void main(String[] args) { new Phylogenetik(); }
}
