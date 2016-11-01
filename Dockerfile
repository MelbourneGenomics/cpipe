FROM ubuntu

# Start using bash
SHELL ["/bin/bash", "-c"]

# Install basic linux tools that we need to make python
RUN apt-get update && apt-get install -y curl make build-essential libssl-dev zlib1g-dev libsqlite3-dev libreadline-dev git libbz2-dev libz-dev libssl-dev gfortran xorg-dev libcurl4-openssl-dev ncurses-dev openjdk-8-jdk

# Copy in the repository
ADD . /opt/cpipe

# Move into the cpipe dir
WORKDIR /opt/cpipe

# Run the install script
RUN ./install.sh

# Run the main script
ENTRYPOINT ["./cpipe"]
