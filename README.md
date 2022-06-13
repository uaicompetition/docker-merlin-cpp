# Dockerized Container

## Preparation in local environment
This version copies statically complied exeuctable inside image.

```
FROM ubuntu:18.04

RUN apt-get -y update && apt-get -y --no-install-recommends install build-essential

ENV MODEL model.uai
ENV EVID evid.uai
ENV QUERY query.uai
ENV RESULT result.uai
ENV LD_LIBRARY_PATH /merlin/boost_1_53_0/stage/lib:$LD_LIBRARY_PATH

COPY . /merlin
WORKDIR /merlin

RUN tar xzf ./boost_1_53_0.tar.gz && \
cd ./boost_1_53_0 && \
./bootstrap.sh && \
./b2 --with-program_options

RUN ln -s /merlin/boost_1_53_0/stage/lib/libboost_program_options.so.1.53.0 /usr/local/lib/libboost_program_options.so
RUN ln -s /merlin/boost_1_53_0/stage/lib/libboost_program_options.a /usr/local/lib/libboost_program_options.a
RUN ln -s /merlin/boost_1_53_0/boost /usr/local/include/boost

RUN make clean && make


CMD ["bash", "-c", "echo $INPUT; /merlin/bin/merlin --input-file $MODEL --evidence-file $EVID --task MAR --algorithm wmb --ibound 4 --output-file $RESULT"] 

```

This step has to be done by each participants and submitted solvers should work without error.

Once we have a working docker script we have two options
(1) push to docker hub (public but easier access)
(2) keep local Docker script and source to build locally

To build a Docker image and running container,

`$ docker build -t junkyu/merlinbin .`
```
# -t is the tag of image, junkyul:merlin follows convention that junkyul is ID of docker hub, 
merlin is the name of Docker image repo. To use Dockerhub, make an ID and create repo, 
install Docker hub Desktop and login (follow instructions https://docs.docker.com/docker-hub/)
```

After building it, we can run container by

`$ docker run -v /home/junkyul/samples:/merlins/problems -e INPUT=/merlins/problems/Alchemy_11.uai -e EVID=/merlins/problems/Alchemy_11.uai.evid junkyul/merlin`

```
# -v is for mounting a volum at host (local hard drive) to Docker container so that we can use Docker as an executable. 
Here we are mounting  /home/junkyul/samples, which contains Alchemy_11 as sample problem to /merlins/problems
# we use environment variable to pass file names to Docker
# ENV MODEL model.uai
# ENV EVID evid.uai
# ENV QUERY query.uai
# ENV RESULT result.uai
```

After running the container we can confirm that the output file is also written in the same folder by
`CMD ["bash", "-c", "echo $MODEL; /merlin/bin/merlin --input-file $MODEL --evidence-file $EVID --task MAR --algorithm wmb --ibound 4 --output-file $RESULT"] `

This step checks whether the Docker image is valid in local environment.
Next, upload image to Docker hub (preferred) or send ZIP file containing all the necessary files


## Preparation in local environment
Once Docker image is working, then we can upload it to dockerhub.
If you don't have dockerhub account, make an account and install docker desktop 
following Step 1 to Step 3 in Quickstart https://docs.docker.com/docker-hub/
Then, simply run the following commands to upload image
```
$ docker login
$ docker build -t <ID>/<REPONAME> .
$ docker push <ID>/<REPONAME>
```


## running Docker containers using Singularity

Our openlab uses singularity for running Docker containers.

Build a SIF singularity image by
`singularity build --bind /home/junkyul/samples:/merlin/problems merlinbin.sif docker://junkyul/merlinbin`

```
# --bind is like -v in docker build
# merlin.sif is the name of image
# docker://junkyul/merlin will directly pull image from docker hub
# if image was local, docker://localhost:5000/merlin:latest
```

After building an image, we can run it by
`singularity run --bind /home/junkyul/samples:/merlin/problems --env INPUT=/merlin/problems/Alchemy_11.uai,EVID=/merlin/problems/Alchemcy_11.uai.evid merlinbin.sif`

```
# the command is like $ singularity run [OPTIONS] merlin.sif
# we also pass filenames through environment variables
# we stored sample files in /home/junkyul/samples
# the bind was writable --writable was not needed?
# after running it we can confirm that the output file was written in /home/junkyul/samples
```

