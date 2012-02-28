CREATE TABLE nodes_to_index
(
	ind	int primary key,
	id	varchar(125) not null unique
)engine=myisam;

CREATE TABLE sets_for_nodes
(
	i	int primary key,
	nodes_ind	int not null,
	start	int not null,
	end	int  not null,
	index idx1 (nodes_ind)
)engine=myisam;

CREATE TABLE parsed
(
	i1	int not null,
	i2	int not null,
	opposite_strand	tinyint not null,
	unique key (i1,i2),
	index idx1 (i2)
)engine=myisam;

CREATE TABLE seeds
(
	i int primary key
)engine=myisam;

# this is for Rfam
CREATE TABLE Rfam_fasta
(
	id	varchar(125) primary key,
	rfam_id  varchar(125) not null,
	rfam_fam varchar(125) not null,
	seq longtext
)engine=myisam;
