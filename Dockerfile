FROM ubuntu
SHELL ["/bin/bash"]
ADD . /opt/cpipe
RUN apt-get update && apt install -y curl make build-essential libssl-dev zlib1g-dev libbz2-dev libsqlite3-dev
RUN ./install.sh
ENTRYPOINT ./cpipe
