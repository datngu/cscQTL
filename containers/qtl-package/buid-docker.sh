## BUILD DOCKER AUTOMATICLY
docker build -t qtl-package:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag qtl-package:v0.0.0 ndatth/qtl-package:v0.0.0
docker push ndatth/qtl-package:v0.0.0
echo DONE


### test docker

docker run -it -v /sigma4:/sigma4 --name qtl ndatth/qtl-package:v0.0.0
docker start qtl
docker attach qtl

# /sigma4/projects/rnaQTL_hg38/containers/qtl-package/bin/QTLtools /usr/bin/

# docker commit --change "LABEL Docker image containing all requirements for running fastQTL" qtl ndatth/qtl-package:v0.0.0
# docker push