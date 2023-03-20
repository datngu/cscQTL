use strict;
use Test::More;
use base qw(Test::Class);

$| = 1;

my $test_count = 0;

sub initial:Test(setup){
	$test_count ++;
	print "-"x80, "\n";
	print "Begin Test $test_count ...\n";
}

sub end:Test(teardown){
	print "End Test $test_count;\n";
	print "-"x80, "\n";
}

sub init:Test(startup){
	print "#"x80, "\n";
}

sub shutdwn:Test(shutdown){
	print "Finished $test_count testcases.\n";
	print "#"x80, "\n";
}

sub test_scanning1:Test(2){
	ok(&existed_list('sample.sam.list') == 1, 'Test for 1st scanning ready');
	my ($expected, $received);
	$expected = &read_list('tmp.list');
	$received = &read_list('sample.sam.list');
	is_deeply($expected, $received, 'Test for 1st scanning') or diag("Test for 1st scanning failed.");
	system 'rm sample.sam.list';
}

sub test_scanning2:Test(2){
	ok(&existed_list('sample.sam.list2') == 1, 'Test for 2nd scanning ready');
	my ($expected, $received);
	$expected = &read_list('tmp.list2');
	$received = &read_list('sample.sam.list2');
	is_deeply($expected, $received, 'Test for 2nd scanning') or diag("Test for 2nd scanning failed.");
	system 'rm sample.sam.list2';
}

sub test_integration:Test(2){
	ok(&existed_list('test.ciri') == 1, 'Integration test ready');
	my ($expected, $received);
	$expected = &read_list('sample.ciri');
	$received = &read_list('test.ciri');
	is_deeply($expected, $received, 'Integration test') or diag("Test for final output failed.");
	system 'rm test.ciri';
	system 'rm test.ciri.log';
}

sub existed_list{
	if(-e $_[0] and -r $_[0]){
		1;
	}else{
		0;
	}
}

sub read_list{
	my %contents;
	open FILE, $_[0] or die "cannot open $_[0]: $!";
	while(<FILE>){
		chomp;
		$contents{$_} ++;
	}
	\%contents;
}

my $dir;
if (rindex($0, "/") >= 0) {
	$dir = substr($0, 0, rindex($0, "/")+1);
} else {
	$dir = "./";
}

print "This script is used to test CIRI2.\n";
unless (-w $dir){
	die "Please make sure the directory $dir is writable.";
}

chdir $dir;

print "Preparing to test ...";
system 'perl ../CIRI2.pl -I sample.sam -O test.ciri -F chr1.fa -D -Q';
print " Done\n";
my $final = Test::Class->runtests();
if($final == 1){
	print "CIRI2 was tested successfully.\n";
}else{
	print "CIRI2 was tested unsuccessfully.\n";
}

