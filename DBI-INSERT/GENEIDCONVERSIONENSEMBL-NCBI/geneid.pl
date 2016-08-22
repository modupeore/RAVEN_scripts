#!/usr/bin/perl -w
# importing ensembl and ncbi gene id conversion table
 use DBI;
 use DBD::mysql;
 open(DATA, "<en_geneid.txt") or die "Can't open en_geneid.txt";

###attributes
 $dsn = 'dbi:mysql:PENGUIN';
 $user = 'modupe';
 $passwd = 'penguin123'; 

###connect
 $dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";
 $i = 0;

#create table
 $syntax = "Create table GENEID(ensembl_geneid varchar(50), entrez_geneid varchar(50),primary key(ensembl_geneid,entrez_geneid))";
#"create table GENE_ID(no_id int not null auto_increment,ensembl_geneid varchar(50), entrez_geneid varchar(50),primary key(no_id))";
 $sth = $dbh->prepare($syntax);
 $sth->execute or die "SQL Error: $DBI::errstr\n";

###INSERT
 $sth = $dbh->prepare("insert into GENEID(ensembl_geneid, entrez_geneid) values (?,?)");
 while (<DATA>){
 	chomp;
 	if(/^ENSGAL.*/){
 		$i++;
 		my ($ensembl, $entrez ) = split /\t/;
 		if (length($entrez)<1){$entrez = "-";}
 		$sth ->execute( $ensembl, $entrez );
	}
}
close DATA;


###statement
 print "\nNumber of lines = $i\n\n";
 $syntax = "select count(*) from GENEID";

 $sth = $dbh->prepare($syntax);

 $sth->execute or die "SQL Error: $DBI::errstr\n";

print "\tTable Information\n";
print "\t======================\n";
 while (@row = $sth->fetchrow_array() ) {print "Row: @row\n";}


##disconnect
$dbh-> disconnect();
exit;
