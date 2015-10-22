---------------------------------------------------------------
-- Mail Order Database; Create Tables Script
-- Chapter 2; Oracle 9i Programming -- A Primer
--            by R. Sunderraman
---------------------------------------------------------------
drop table zipcodes cascade constraints;

-- "cascade constraints" means drop any constraints that other
-- tables may have on this table when it is deleted

create table zipcodes (
  zip      number(5),
  city     varchar2(30),
  primary key (zip));

drop table employees cascade constraints;
create table employees (
  eno      number(4) not null primary key, 
  ename    varchar2(30),
  zip      number(5) references zipcodes,
  hdate    date);

drop table parts cascade constraints;
create table parts(
  pno      number(5) not null primary key,
  pname    varchar2(30),
  qoh      integer check(qoh >= 0),
  price    number(6,2) check(price >= 0.0),
  olevel   integer);

drop table customers cascade constraints;
create table customers (
  cno      number(5) not null primary key,
  cname    varchar2(30),
  street   varchar2(30),
  zip      number(5) references zipcodes,
  phone    char(12));
 
drop table orders cascade constraints;
create table orders (
  ono      number(5) not null primary key,
  cno      number(5) references customers,
  eno      number(4) references employees,
  received date,
  shipped  date);

drop table odetails cascade constraints;
create table odetails (
  ono      number(5) not null references orders,
  pno      number(5) not null references parts,
  qty      integer check(qty > 0),
  primary key (ono,pno));

------------- now populate relations ---------

insert into  zipcodes values
  (67226,'Wichita');
insert into  zipcodes values
  (60606,'Fort Dodge');
insert into  zipcodes values
  (50302,'Kansas City');
insert into  zipcodes values
  (54444,'Columbia');
insert into  zipcodes values
  (66002,'Liberal');
insert into  zipcodes values
  (61111,'Fort Hays');

insert into employees values
  (1000,'Jones',67226,'12-DEC-95');
insert into employees values
  (1001,'Smith',60606,'01-JAN-92');
insert into employees values
  (1002,'Brown',50302,'01-SEP-94');

insert into parts values
  (10506,'Land Before Time I',200,19.99,20);
insert into parts values
  (10507,'Land Before Time II',156,19.99,20);
insert into parts values
  (10508,'Land Before Time III',190,19.99,20); 
insert into parts values
  (10509,'Land Before Time IV',60,19.99,20);
insert into parts values
  (10601,'Sleeping Beauty',300,24.99,20);
insert into parts values
  (10701,'When Harry Met Sally',120,19.99,30);
insert into parts values
  (10800,'Dirty Harry',140,14.99,30);
insert into parts values
  (10900,'Dr. Zhivago',100,24.99,30);

insert into customers values
  (1111,'Charles','123 Main St.',67226,'316-636-5555');
insert into customers values
  (2222,'Bertram','237 Ash Avenue',67226,'316-689-5555');
insert into customers values
  (3333,'Barbara','111 Inwood St.',60606,'316-111-1234');

insert into orders values
  (1020,1111,1000,'10-DEC-94','12-DEC-94');
insert into orders values
  (1021,1111,1000,'12-JAN-95','15-JAN-95');
insert into orders values
  (1022,2222,1001,'13-FEB-95','20-FEB-95');
insert into orders values
  (1023,3333,1000,'20-JUN-97',null);

insert into odetails values
  (1020,10506,1);
insert into odetails values
  (1020,10507,1);
insert into odetails values
  (1020,10508,2);
insert into odetails values
  (1020,10509,3);
insert into odetails values
  (1021,10601,4);
insert into odetails values
  (1022,10601,1);
insert into odetails values
  (1022,10701,1);
insert into odetails values
  (1023,10800,1);
insert into odetails values
  (1023,10900,1);
