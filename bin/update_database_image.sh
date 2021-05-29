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

# default tag
get_latest_release() {
  curl --silent "https://api.github.com/repos/$1/releases/latest" | # Get latest release from GitHub api
    grep '"tag_name":' |                                            # Get tag line
    sed -E 's/.*"([^"]+)".*/\1/'                                    # Pluck JSON value
}
DEF_TAG=$(get_latest_release fmalmeida/bacannot)

# default reply
REPLY="N"

# actual code
read -p "Do you really want to build it locally? (y/N) " -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then

    echo "Set a different release? Default is ${DEF_TAG}"
    echo "It will create the image fmalmeida/bacannot:${DEF_TAG}"
    read -p "Your answer? (y/N) " -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        read -p "Which tag? " -r
        echo    # (optional) move to a new line
        DEF_TAG=$REPLY
    else
      echo
    fi
    echo
    echo "# Starting docker image fmalmeida/bacannot:${DEF_TAG} local build!"
    wget --quiet -O my_dockerfile https://github.com/fmalmeida/bacannot/raw/master/docker/Dockerfile_bacannot
    docker build -t fmalmeida/bacannot:${DEF_TAG} -f my_dockerfile .
    rm my_dockerfile
    echo
    echo "Image built!"
    echo

    echo "# Environment cleaning"
    echo
    echo "Do you want to remove the base image fmalmeida/bacannot_main_tools?"
    echo "When updating the databases again, if available, docker will automatically"
    echo "skips the base image download."
    read -p "Your answer? (y/N) " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
      docker rmi -f fmalmeida/bacannot_main_tools
    else
      echo
    fi
else
    echo "You've chosen not to build! Giving up"
    echo
fi
