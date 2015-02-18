#!/bin/bash

check_exit_code() {
	# $1 is a message about what the last command was
	# e.g. "run yotools deoligo command on empty folder"
	if [ "$?" -ne 0 ]
	then
		echo ""
		echo "	Fail. Was trying to $1"
		echo ""
		exit -1
	fi
	echo ":) Successfully managed to $1"
}

check_directory() {
	# $1 is the directory to check
	# $2 is the command that should have created it
	if [ ! -d $1 ]
	then
		echo "	Failed to find directory $1."
		echo "	Should have been created by $2."
		exit -1
	fi
	echo ":) Verified creation of directory $1"
}

check_file() {
	# Verifies that file exists and is nonempty
	# $1 is the file to check
	# $2 is the command that should have created it
	if [ ! -s $1 ]
	then
		echo "	Failed to find file $1."
		echo "	Should have been created by $2."
		exit -1
	fi
	echo ":) Verified creation of file $1"
}

check_files_match() {
	diff $1 $2 > /dev/null
	check_exit_code "diff $1 $2"
}
