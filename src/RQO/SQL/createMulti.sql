drop table client cascade constraints;
create table client (
   clientNo char(5) primary key,
	fName    char(8),
	lName    char(10),
	telNo    char(15),
	prefType char(8),
	maxrent  integer
);

insert into client values ('CR76', 'John', 'Kay',     '0207-774-5632','Flat', '425');
insert into client values ('CR56', 'Aline','Stewart', '0141-848-1825','Flat', '350');
insert into client values ('CR74', 'Mike', 'Richie',  '01475-392178', 'House','750' );
insert into client values ('CR62', 'Mary', 'Tregear', '01224-196720', 'Flat', '600' );

drop table viewing cascade constraints;
create table viewing (
   clientNo    char(5) references client,
	propertyNo  char(5),
	viewDate    date,
	komment     char(20),
	primary key (clientNo, propertyNo)
);

insert into viewing values ('CR56', 'PA14', '24-MAY-04', 'too small');
insert into viewing values ('CR76', 'PG4',  '20-APR-04', 'too remote');
insert into viewing values ('CR56', 'PG4',  '26-MAY-04', null );
insert into viewing values ('CR62', 'PA14', '14-MAY-04', 'no dining room');
insert into viewing values ('CR56', 'PG36', '28-APR-04', null );


