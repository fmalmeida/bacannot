#!/bin/bash

# Simple script to start the Bacannot server
# that concatenates and brings together all
# of its outputs.
#
# Author: Felipe M. Almeida (almeidafmarques@outlook.com)

# Help
Help()
{

	# Display Help
	echo
	echo "Simple script for starting the bacannot server for its output interrogation"
	echo
	echo "Syntax: run_jbrowse.sh [-h|s|d] [-p <number>]"
	echo "options:"
	echo
	echo "h					Print this help"
	echo "p					Set out a custom PORT for server listening. [ Default: 3838 ]"
	echo "s					Start bacannot server"
	echo "d					Only download the required docker image"
	echo ""
	echo
}

# Download image
Download()
{
	# Run docker pull
	docker pull fmalmeida/bacannot:server
}

# Start server
Start()
{

	# Tells user
	echo "The server has started in: http://localhost:${PORT}/"

	# Start server in current directory
	echo "When finished, run the command:"
	echo -n "	docker rm -f "
	docker run -v $(pwd):/work -d --rm --platform linux/amd64 -p 4567:4567 -p 3838:"$PORT" fmalmeida/bacannot:server

}

# Default
PORT="3838"

# Get the options
if [ "${1}" == "" ] ; then
	Help
	exit
fi
while getopts "hsd:p" option; do
   case $option in
      h) # display Help
         Help
         exit;;
			p) # get custom PORT
				 PORT="$OPTARG"
				 ;;
			s) # Start server
				 Start
				 exit;;
			d) # Only download image
				 Download
				 ;;
   esac
done
