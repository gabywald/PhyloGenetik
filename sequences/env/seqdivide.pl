#! usr/bin/perl -w

## ce petit programme pour fabriquer les tronçons 1-100 ; 101-200 ; 201-300 et 301-400

opendir(DIRHANDLE,".")||die "ERROR: impossible de lire le répertoire courant\n"; 
foreach $file (readdir(DIRHANDLE)){ 
	if ($file =~ /.fasta/) {
		$compteur = 1;
		$tampon = "";
		$firstline = "";
		open SEQFILE,$file;
		@lines = <SEQFILE>;
		chomp(@lines);
		close SEQFILE;
		foreach $line (@lines) {
			#print $line."\n";
			$line =~ /\n/;
			$line =~ /\r/;
			if ($line =~ />/) {
				$firstline = substr($line,1,length($line)-1);
			} else {
				for ($i = 0;$i<length($line);$i++) {
					$char = substr($line,$i,1);
					$tampon = $tampon.$char;
					#print $char.$compteur."\n";
					$compteur++;
					if ($compteur == 101) {
						open TESTFILE,">../env1-100/$file";
						$count = ($compteur-100+1);
						print TESTFILE ">$count-$compteur :: $firstline";
						print TESTFILE $tampon;
						close TESTFILE;
						$tampon="";
					}
					if ($compteur == 201) {
						open TESTFILE,">../env101-200/$file";
						$count = ($compteur-100+1);
						print TESTFILE ">$count-$compteur :: $firstline";
						print TESTFILE $tampon;
						close TESTFILE;
						$tampon="";
					}
					if ($compteur == 301) {
						open TESTFILE,">../env201-300/$file";
						$count = ($compteur-100+1);
						print TESTFILE ">$count-$compteur :: $firstline";
						print TESTFILE $tampon;
						close TESTFILE;
						$tampon="";
					}
					if ($compteur == 401) {
						open TESTFILE,">../env301-400/$file";
						$count = ($compteur-100+1);
						print TESTFILE ">$count-$compteur :: $firstline";
						print TESTFILE $tampon;
						close TESTFILE;
						$tampon="";
					}
				}
			}
		}
	}
} 
closedir DIRHANDLE;


