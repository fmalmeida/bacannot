#!/bin/bash

# log message
cat << EOF
# Bacannot shell script for updating the docker image which contains the database files
#
# It is useful to maintain the databases up-to-date. By default, with github actions, the
# image is updated in the first day of each month, however, this script enables that you
# update the database image any time.
#
# When building (locally or in github actions), the image will always update the databases.
#
# Author: Felipe M. Almeida (almeidafmarques@outlook.com)
#
# The script will now begin the image update!
#
# Remember: docker must be available in \$PATH


EOF

REPLY="N"
read -p "Do you really want to build it locally? (y/N) " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Starting docker local build!"
    wget --quiet -O my_dockerfile https://github.com/fmalmeida/bacannot/raw/master/docker/Dockerfile_bacannot
    docker build -t fmalmeida/bacannot:latest -f my_dockerfile .
    docker rmi -f fmalmeida/bacannot_main_tools
    echo
    echo "Image built!"
    echo
else
    echo "You've chosen not to build! Giving up"
    echo
fi
