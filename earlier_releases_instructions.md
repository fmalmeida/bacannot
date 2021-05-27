# Building docker images of earlier (or specific) releases

Since the tools wrapped in the docker images may change as well as the image file structure, the latest docker images available in Dockerhub may not be compatible with earlier releases. For using earlier releases, users must build the specific images locally as shown below.

> Remember that this will build images named the same way as those on Dockerhub. So if you build and earlier release image and after that want to run the latest pipeline version, will must remember to pull the latest image from dockerhub. The images for the latest pipeline versions will always be in Dockerhub.

```bash
# clone the repo
git clone https://github.com/fmalmeida/bacannot.git

# checkout to specific release tag
git checkout tags/<tag_name>

# enter docker dir
cd docker

# build images
docker build -t fmalmeida/bacannot:latest    -f Dockerfile_bacannot  .
docker build -t fmalmeida/bacannot:kofamscan -f Dockerfile_kofamscan .
docker build -t fmalmeida/bacannot:jbrowse   -f Dockerfile_jbrowse   .
docker build -t fmalmeida/bacannot:renv      -f Dockerfile_renv      .
docker build -t fmalmeida/bacannot:server    -f Dockerfile_server    .
```
