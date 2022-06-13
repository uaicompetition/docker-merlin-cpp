FROM ubuntu:18.04

# install OS
RUN apt-get -y update && apt-get -y --no-install-recommends install build-essential

# give dummy names for the environment variables that will be overwritten when running
ENV MODEL model.uai
ENV EVID evid.uai
ENV QUERY query.uai
ENV RESULT result.uai

COPY . /merlin
WORKDIR /merlin		# this can be any name that stores the files 

CMD ["bash", "-c", "echo $MODEL; /merlin/merlin-static --input-file $MODEL --evidence-file $EVID --task MAR --algorithm wmb --ibound 4 --output-file $RESULT"] 


