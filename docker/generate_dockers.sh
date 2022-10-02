source set_version.sh
for dir in $(ls -d */) ; do ( cd $dir && bash build.sh $NEW_VERSION ) ; done