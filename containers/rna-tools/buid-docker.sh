## BUILD DOCKER AUTOMATICLY
docker build -t rna-tools:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag rna-tools:v0.0.0 ndatth/rna-tools:v0.0.0
docker push ndatth/rna-tools:v0.0.0
echo DONE


### test docker
docker run -it --rm -v /sigma4:/sigma4 --name rna rna-tools:v0.0.0
docker run -it --rm -v /sigma4:/sigma4 --name rna ndatth/rna-tools:v0.0.0
docker start rna
docker attach rna