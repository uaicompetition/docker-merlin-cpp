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


