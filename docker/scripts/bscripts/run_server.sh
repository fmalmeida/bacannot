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
	echo "k                 Kill bacannot server"
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
	echo -n "	docker rm -f ServerBacannot or ./run_server.sh -k"
	docker run -v $(pwd):/work -d --rm --platform linux/amd64 -p "$PORT":3838 -p 4567:4567 --name ServerBacannot fmalmeida/bacannot:server &> /dev/null
	echo

}

# Kill server
KillServer()
{
	docker rm -f ServerBacannot
}

# Default
PORT="3838"
EXEC='false'

# Get the options
if [ "${1}" == "" ] ; then
	Help
	exit
fi
while getopts "hskdp:" option; do
   case $option in
      h) # display Help
         Help
         exit
		 ;;
	  k) # kill server
	  	 KillServer
		 exit
		 ;;
	  p) # get custom PORT
		 PORT="$OPTARG"
		 ;;
	  s) # Start server
		 EXEC='true'
		 ;;
	  d) # Only download image
		 Download
		 ;;
   esac
done

if [ "$EXEC" == 'true' ]; then
	Start
fi