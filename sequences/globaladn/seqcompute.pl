#! usr/bin/perl -w

## ce petit programme pour fabriquer les tronçons 300-600 et 600-900

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
					if ($compteur == 299) {
						$tampon="";
					}
					if ($compteur == 600) {
						open TESTFILE,">../globaladn300-600/$file";
						$count = ($compteur-300);
						print TESTFILE ">$count-$compteur :: $firstline";
						print TESTFILE $tampon;
						close TESTFILE;
						$tampon=substr($line,length($tampon)-1,1);
					}
					if ($compteur == 900) {
						open TESTFILE,">../globaladn600-900/$file";
						$count = ($compteur-300);
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


